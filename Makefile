# Universal macOS build for Apple Silicon + Intel
#
# Usage:
#   make universal            # builds a universal binary at $(OUT)/$(UNIVERSAL)
#   make build                # builds per-arch release bins
#   make fmt clippy test      # hygiene targets
#   make clean
#
# Customize these if your crate/binary name differs
BIN ?= tagbam
OUT ?= dist
UNIVERSAL ?= $(BIN)-universal

# Targets
AARCH := aarch64-apple-darwin
X86   := x86_64-apple-darwin

CARGO := cargo
LIPO  := lipo
FILE  := file

# Wider macOS compatibility when lipo-ing
MACOSX_DEPLOYMENT_TARGET ?= 11.0
export MACOSX_DEPLOYMENT_TARGET

# Prevent any shell RUSTFLAGS leaking into cross builds (config rustflags still apply)
export RUSTFLAGS=

.PHONY: all setup build universal clean fmt clippy test

all: universal

setup:
	rustup target add $(AARCH) $(X86)

build: setup
	$(CARGO) build --release --target $(AARCH)
	$(CARGO) build --release --target $(X86)

universal: build
	mkdir -p $(OUT)
	$(LIPO) -create \
	  -output $(OUT)/$(UNIVERSAL) \
	  target/$(AARCH)/release/$(BIN) \
	  target/$(X86)/release/$(BIN)
	$(FILE) $(OUT)/$(UNIVERSAL)

clean:
	rm -rf $(OUT)
	$(CARGO) clean

fmt:
	$(CARGO) fmt --all

clippy:
	$(CARGO) clippy --all-targets -- -D warnings

test:
	$(CARGO) test
