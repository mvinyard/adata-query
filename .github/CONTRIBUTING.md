# Contributing to adata-query

Thank you for considering contributing to adata-query! This document provides guidelines and instructions for contributing.

## Quick Start

1. **Fork and Clone**
   ```bash
   git clone https://github.com/YOUR_USERNAME/adata-query.git
   cd adata-query
   ```

2. **Set Up Development Environment**
   ```bash
   # Install dependencies with development tools
   uv sync --extra dev
   
   # Optional: Install pre-commit hooks
   uv run pre-commit install
   ```

3. **Create a Branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

## Development Workflow

### Running Tests

```bash
# Run all tests
uv run pytest tests/ -v

# Run specific test file
uv run pytest tests/test_locator.py -v

# Run with coverage
uv run pytest tests/ --cov=adata_query --cov-report=term-missing
```

### Code Quality

```bash
# Format code
uv run ruff format src/ tests/

# Check for linting issues
uv run ruff check src/ tests/

# Auto-fix linting issues
uv run ruff check --fix src/ tests/
```

### Pre-commit Checks

If you installed pre-commit hooks, they will run automatically before each commit. To run them manually:

```bash
uv run pre-commit run --all-files
```

## Submitting Changes

### 1. Write Tests

- Add tests for new features in the `tests/` directory
- Ensure all tests pass locally
- Aim for at least 85% code coverage

### 2. Update Documentation

- Add docstrings to new functions/classes
- Update README.md if needed
- Add examples for new features

### 3. Create a Pull Request

1. **Push your branch:**
   ```bash
   git push origin feature/your-feature-name
   ```

2. **Open a PR on GitHub:**
   - Provide a clear, descriptive title (min 10 characters)
   - Include a detailed description of your changes
   - Reference any related issues

3. **Wait for CI checks:**
   - All tests must pass
   - Code coverage should be ≥85%
   - Linting checks must pass

### PR Review Process

- A maintainer will review your PR
- Address any feedback or requested changes
- Once approved, your PR will be merged!

## CI/CD Checks

Your PR will automatically be checked for:

✅ **Tests** - Run on Python 3.11, 3.12, and 3.13
✅ **Code Coverage** - Must be at least 85%
✅ **Linting** - Code must pass ruff checks
✅ **Formatting** - Code must be formatted with ruff
✅ **PR Metadata** - Title and description must be descriptive

## Coding Standards

### Style Guide

- Follow PEP 8
- Use type hints where appropriate
- Keep imports alphabetically sorted
- Use descriptive variable names

### Import Organization

```python
# -- import packages: ---------------------------------------------------------
import package_name

# -- set type hints: ----------------------------------------------------------
from typing import Any, Dict, List

# -- import local dependencies: -----------------------------------------------
from ._module import function
```

### Docstring Format

```python
def function_name(param: str) -> int:
    """Brief description.
    
    Longer description if needed.
    
    Args:
        param: Description of parameter.
        
    Returns:
        Description of return value.
        
    Example:
        >>> function_name("example")
        42
    """
```

## Testing Guidelines

### Test Structure

- One test file per module
- Group related tests in classes
- Use descriptive test names

Example:
```python
class TestMyFeature:
    """Test suite for MyFeature."""
    
    def test_basic_functionality(self, sample_fixture):
        """Test that basic functionality works as expected."""
        result = my_feature(sample_fixture)
        assert result == expected_value
```

### Fixtures

- Use fixtures defined in `tests/conftest.py`
- Add new fixtures there if they'll be reused
- Keep test data minimal but representative

## Getting Help

- Open an issue for bugs or feature requests
- Ask questions in discussions
- Check existing issues and PRs first

## Code of Conduct

- Be respectful and inclusive
- Provide constructive feedback
- Focus on the code, not the person

Thank you for contributing to adata-query! 🎉

