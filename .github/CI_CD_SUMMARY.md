# CI/CD Setup Summary

This document provides a quick overview of the CI/CD setup for `adata-query`.

## ✅ What Was Added

### GitHub Actions Workflows

1. **`tests.yml`** - Main Test Suite
   - Runs on push to `main`/`dev` and on PRs
   - Tests on Python 3.11, 3.12, and 3.13
   - Generates coverage reports (90% coverage achieved!)
   - Uploads to Codecov (optional, requires token)
   - Runs linting checks (non-blocking)

2. **`pr-check.yml`** - Pull Request Validation
   - Validates PR metadata (title, description)
   - Runs code formatting checks (non-blocking)
   - Runs linting (non-blocking)
   - Ensures test coverage ≥85%
   - Checks if tests were added for new code

### Configuration Files

- **`.pre-commit-config.yaml`** - Pre-commit hooks for local development
- **`.github/CONTRIBUTING.md`** - Contributor guidelines
- **`pyproject.toml`** - Updated with:
  - Build system configuration
  - Ruff linter configuration
  - Dev dependencies
  - Test dependencies

## 🚀 Quick Start for Contributors

```bash
# Clone and setup
git clone https://github.com/YOUR_USERNAME/adata-query.git
cd adata-query

# Install with dev dependencies
uv sync --extra dev

# Optional: Install pre-commit hooks
uv run pre-commit install

# Run tests
uv run pytest tests/ -v

# Check code quality
uv run ruff check --fix src/ tests/
uv run ruff format src/ tests/
```

## 📋 What Checks Run on PRs

When you open a PR, these checks run automatically:

### Required (Must Pass)
- ✅ Tests on Python 3.11, 3.12, 3.13
- ✅ Test coverage ≥85%
- ✅ PR has descriptive title (≥10 characters)

### Recommended (Non-Blocking)
- ⚠️ Code formatting (ruff format)
- ⚠️ Linting (ruff check)
- ⚠️ Tests added for new code
- ⚠️ PR has description

## 🔧 Commands

### Running Tests
```bash
# All tests
uv run pytest tests/ -v

# With coverage
uv run pytest tests/ --cov=adata_query --cov-report=term-missing

# Specific test
uv run pytest tests/test_fetcher.py::TestFetchFunction::test_fetch_simple -v
```

### Code Quality
```bash
# Auto-fix linting issues
uv run ruff check --fix src/ tests/

# Format code
uv run ruff format src/ tests/

# Check without fixing
uv run ruff check src/ tests/
uv run ruff format --check src/ tests/
```

### Pre-commit Hooks
```bash
# Install hooks (runs checks before each commit)
uv run pre-commit install

# Run manually
uv run pre-commit run --all-files

# Skip hooks for a commit (not recommended)
git commit --no-verify
```

## 📊 Current Status

- **Test Suite**: 76 tests (70 passed, 6 skipped)
- **Code Coverage**: 90%
- **Python Versions**: 3.11, 3.12, 3.13
- **Package Manager**: uv (fast, modern)

## 🎯 Best Practices

1. **Before Creating a PR**:
   ```bash
   uv run pytest tests/ -v          # Run tests
   uv run ruff format src/ tests/   # Format code
   uv run ruff check --fix src/ tests/  # Fix linting
   ```

2. **Writing Tests**:
   - Add tests for new features
   - Aim for ≥85% coverage
   - Use existing fixtures in `tests/conftest.py`

3. **Code Style**:
   - Follow PEP 8
   - Use type hints
   - Keep imports sorted alphabetically
   - Add docstrings to public functions

4. **PR Guidelines**:
   - Descriptive title and description
   - Reference related issues
   - Keep changes focused
   - Respond to review feedback

## 🔗 Useful Links

- [GitHub Actions Docs](https://docs.github.com/en/actions)
- [pytest Documentation](https://docs.pytest.org/)
- [Ruff Documentation](https://docs.astral.sh/ruff/)
- [uv Documentation](https://docs.astral.sh/uv/)

## 🆘 Troubleshooting

### Tests fail in CI but pass locally
- Check Python version matches
- Ensure `uv.lock` is up to date: `uv lock`
- Clear cache: `rm -rf .pytest_cache`

### Coverage too low
- Run: `uv run pytest --cov=adata_query --cov-report=html`
- Open: `htmlcov/index.html` in browser
- Add tests for uncovered lines

### Linting errors
- Auto-fix: `uv run ruff check --fix src/ tests/`
- Format: `uv run ruff format src/ tests/`
- Some errors may need manual fixing

### Pre-commit hooks failing
- Update hooks: `uv run pre-commit autoupdate`
- Skip temporarily: `git commit --no-verify` (not recommended)

## 📝 Notes

- Linting checks are **non-blocking** in CI (show warnings but don't fail)
- This allows gradual improvement without blocking development
- Coverage requirement (85%) ensures code quality
- Tests must pass on all Python versions

## 🎉 Next Steps

To fully enable CI/CD:

1. **Push these changes to GitHub**
2. **Enable GitHub Actions** (if not already enabled)
3. **Optional: Set up Codecov**
   - Sign up at [codecov.io](https://codecov.io)
   - Add `CODECOV_TOKEN` to repository secrets
4. **Optional: Add branch protection rules**
   - Require PR reviews
   - Require status checks to pass
5. **Optional: Add badges to README**
   ```markdown
   [![Tests](https://github.com/YOUR_USERNAME/adata-query/actions/workflows/tests.yml/badge.svg)](https://github.com/YOUR_USERNAME/adata-query/actions/workflows/tests.yml)
   ```

## 📬 Questions?

- Open an issue for bugs/features
- Check the workflows README: `.github/workflows/README.md`
- See contributing guidelines: `.github/CONTRIBUTING.md`

