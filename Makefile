.PHONY: all build release debug check test clean fmt install help

BIN     := target/release/alphabet
INPUT   ?= input.10kb.fasta
OUTDIR  ?= .
THREADS ?= 8

all: release

build: release

release:
	cargo build --release

debug:
	cargo build

check:
	cargo check

test:
	cargo test --release

fmt:
	cargo fmt

clean:
	cargo clean

install: release
	cargo install --path . --force

# End-to-end pipeline: discover → cut → annotate.
# Usage: make run INPUT=path/to/arrays.10kb.fasta OUTDIR=out THREADS=96
run: release
	$(BIN) run $(INPUT) -o $(OUTDIR) -t $(THREADS)

help:
	@echo "Alphabet — make targets"
	@echo "  make                       build release binary (default)"
	@echo "  make debug                 build debug binary"
	@echo "  make check                 cargo check"
	@echo "  make test                  run tests"
	@echo "  make fmt                   cargo fmt"
	@echo "  make clean                 cargo clean"
	@echo "  make install               cargo install --path ."
	@echo "  make run INPUT=... OUTDIR=... THREADS=..."
	@echo ""
	@echo "Binary: $(BIN)"
	@echo "CLI:    $(BIN) --help"
