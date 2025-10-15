.RECIPEPREFIX = +
SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c

# Detect system Python
PYTHON_SYSTEM := $(shell command -v python3 || command -v python || echo "")
ifeq ($(PYTHON_SYSTEM),)
$(error "No Python found. Please install Python 3.")
endif

# Path to virtualenv Python
PYTHON := venv/bin/python3
PY = src tests

.PHONY: venv
venv: .venv

.venv:
+ $(PYTHON_SYSTEM) -m venv venv
+ touch $@

.PHONY: install
install: .install

.install: requirements.txt .venv
+ $(PYTHON) -m pip install -r requirements.txt
+ $(PYTHON) -m pip install -r requirements-dev.txt
+ touch $@

.PHONY: ruff
ruff: .install
+ $(PYTHON) -m ruff check $(PY)

.PHONY: lint
lint: .install
+ $(PYTHON) -m pylint $(PY)

.PHONY: check
check: .install
+ $(PYTHON) -m mypy $(PY)

.PHONY: test
test: .install
+ $(PYTHON) -m pytest -s tests
