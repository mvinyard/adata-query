# CI/CD Workflows

This directory contains GitHub Actions workflows for continuous integration and deployment.

## Workflows

### 1. `tests.yml` - Main Test Suite

**Triggers:**
- Push to `main` or `dev` branches
- Pull requests to `main` or `dev` branches

**What it does:**
- Runs tests on Python 3.11, 3.12, and 3.13
- Installs dependencies using `uv` for fast setup
- Runs pytest with coverage reporting
- Uploads coverage to Codecov (Python 3.11 only)
- Runs linting checks with ruff

**Status:** This is the main workflow that validates all code changes.

### 2. `pr-check.yml` - Pull Request Validation

**Triggers:**
- When a PR is opened, synchronized, or reopened
- Only runs on non-draft PRs

**What it does:**
- **Code Quality Checks:**
  - Validates code formatting with ruff
  - Runs linter checks
  - Ensures test coverage is at least 85%
  
- **PR Metadata Checks:**
  - Ensures PR title is descriptive (min 10 characters)
  - Checks for PR description
  
- **Test Coverage Check:**
  - Detects if new Python files were added
  - Warns if corresponding tests are missing (non-blocking)
  
- **Coverage Comment:**
  - Posts coverage report in the PR summary

**Status:** This workflow provides additional validation specific to PRs.

## Setting Up CI/CD

### For Repository Owners

1. **Enable GitHub Actions** (if not already enabled)
   - Go to repository Settings → Actions → General
   - Ensure "Allow all actions and reusable workflows" is selected

2. **Optional: Set up Codecov**
   - Sign up at [codecov.io](https://codecov.io)
   - Add your repository
   - Add `CODECOV_TOKEN` to repository secrets (Settings → Secrets and variables → Actions)

3. **Branch Protection Rules** (recommended)
   - Go to Settings → Branches
   - Add rule for `main` branch:
     - ✅ Require a pull request before merging
     - ✅ Require status checks to pass before merging
       - Select: `Test on Python 3.11`, `Test on Python 3.12`, `Test on Python 3.13`
       - Select: `Validate PR`, `Check PR Metadata`
     - ✅ Require conversation resolution before merging

### For Contributors

1. **Local Development Setup**
   ```bash
   # Install dependencies
   uv sync --extra dev
   
   # Install pre-commit hooks (optional but recommended)
   uv run pre-commit install
   ```

2. **Before Creating a PR**
   ```bash
   # Run tests locally
   uv run pytest tests/ -v
   
   # Check code formatting
   uv run ruff format src/ tests/
   
   # Run linter
   uv run ruff check src/ tests/ --fix
   
   # Check coverage
   uv run pytest tests/ --cov=adata_query --cov-report=term-missing
   ```

3. **Using Pre-commit Hooks** (optional)
   ```bash
   # Automatically run checks before each commit
   uv run pre-commit install
   
   # Run manually on all files
   uv run pre-commit run --all-files
   ```

## Workflow Status

You can check the status of workflows:
- In the "Actions" tab of the repository
- On the PR page (checks will appear at the bottom)
- In commit status (green checkmark or red X)

## Adding Badges to README

Add these badges to your `README.md`:

```markdown
[![Tests](https://github.com/YOUR_USERNAME/adata-query/actions/workflows/tests.yml/badge.svg)](https://github.com/YOUR_USERNAME/adata-query/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/YOUR_USERNAME/adata-query/branch/main/graph/badge.svg)](https://codecov.io/gh/YOUR_USERNAME/adata-query)
```

Replace `YOUR_USERNAME` with your GitHub username.

## Troubleshooting

### Tests Fail Locally but Pass in CI (or vice versa)

- Ensure you're using the same Python version
- Check that `uv.lock` is up to date: `uv lock`
- Clear pytest cache: `rm -rf .pytest_cache`

### Coverage Fails

- Run locally: `uv run pytest tests/ --cov=adata_query --cov-fail-under=85`
- Add tests for uncovered code
- Check `htmlcov/index.html` for detailed coverage report

### Linter Errors

- Auto-fix most issues: `uv run ruff check --fix src/ tests/`
- Format code: `uv run ruff format src/ tests/`

## Configuration Files

- `.github/workflows/*.yml` - Workflow definitions
- `.pre-commit-config.yaml` - Pre-commit hook configuration
- `pyproject.toml` - Tool configuration (ruff, pytest, etc.)

## Future Enhancements

Potential additions:
- Release automation workflow
- Documentation building and deployment
- Performance benchmarking
- Security scanning (Dependabot, CodeQL)
- Automated changelog generation

