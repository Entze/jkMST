
m_CPLEX_STUDIO_BINARIES := $(shell julia --color=yes src/locateCPLEX.jl | tail -n 1)

all: install

install: build

build: export CPLEX_STUDIO_BINARIES=${m_CPLEX_STUDIO_BINARIES}
build:
	julia --project=@. --color=yes --eval "import Pkg; Pkg.build();"

update: export CPLEX_STUDIO_BINARIES=${m_CPLEX_STUDIO_BINARIES}
update:
	julia --project=@. --color=yes --eval "import Pkg; Pkg.update();"

clean:
	$(strip $(RM) $(filter-out $(NEVERDELETE),))

NEVERDELETE =$(strip $(wildcard src/*) $(wildcard data/*.dat) jkMST LICENSE LICENCE Makefile Project.toml Manifest.toml)

.PHONY = build install update all