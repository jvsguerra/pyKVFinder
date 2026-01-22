.PHONY: install
install: ## Install the virtual environment and install the pre-commit hooks
	@echo "🚀 Creating virtual environment using uv"
	@uv sync
	@uv pip install .[dev,docs]
	@echo "🔧 Installing pre-commit hooks"
	@uv run pre-commit install

.PHONY: check
check: ## Run code quality tools.
	@echo "🚀 Checking lock file consistency with 'pyproject.toml'"
	@uv lock --locked
	@echo "🚀 Linting code: Running pre-commit"
	@uv run pre-commit run -a
	@echo "🚀 Static type checking: Running mypy"
	@uv run mypy

.PHONY: build
build: build ## Build wheel file
	@echo "🚀 Creating wheel file"
	@uvx --from build pyproject-build --installer uv

.PHONY: clean-build
clean-build: ## Clean build artifacts
	@echo "🚀 Removing build artifacts"
	@uv run python -c "import shutil; import os; shutil.rmtree('dist') if os.path.exists('dist') else None"

.PHONY: tests
tests: unit integration ## Test the code
	@echo "🚀 Testing package from tests/"
	@uv run pytest

unit: ## Run unit tests
	@echo "🚀 Running unit tests from tests/unit/"
	@uv run pytest tests/unit/ -v

integration: ## Run integration tests
	@echo "🚀 Running integration tests from tests/integration/"
	@uv run pytest tests/integration/ -v

.PHONY: clean
clean: ## Remove all untracked files (except .venv and data/)
	@echo "🚀 Removing untracked files"
	@git clean -fdx -e .venv -e data/ -e results/ .

.PHONY: docs
docs: ## Build HTML documentation
	@cd docs && make html
