# Test Suite for adata-query

This directory contains the test suite for the `adata-query` package.

## Running Tests

### With `uv` (recommended)

First, sync the test dependencies:

```bash
uv sync --extra test
```

Then run tests:

```bash
# Run all tests
uv run pytest tests/

# Run with verbose output
uv run pytest tests/ -v

# Run with coverage
uv run pytest tests/ --cov=adata_query --cov-report=term-missing --cov-report=html

# Run specific test file
uv run pytest tests/test_locator.py -v

# Run specific test
uv run pytest tests/test_locator.py::TestAnnDataLocator::test_init_default -v
```

### With standard pip/pytest

To run all tests:

```bash
pytest tests/
```

To run tests with verbose output:

```bash
pytest tests/ -v
```

To run tests with coverage (requires pytest-cov):

```bash
pytest tests/ --cov=adata_query --cov-report=term-missing --cov-report=html
```

To run a specific test file:

```bash
pytest tests/test_locator.py -v
```

To run a specific test:

```bash
pytest tests/test_locator.py::TestAnnDataLocator::test_init_default -v
```

## Test Structure

- `conftest.py` - Shared fixtures and test configuration
- `test_locator.py` - Tests for the `_locator` module (AnnDataLocator class and locate function)
- `test_formatter.py` - Tests for the `_formatter` module (DataFormatter class and format_data function)
- `test_fetcher.py` - Tests for the `_fetcher` module (AnnDataFetcher class and fetch function)

## Test Coverage

The test suite includes:

### Locator Tests (18 tests)
- Initialization with default and custom searchable parameters
- Attribute storage and retrieval
- Data intake from AnnData objects
- Cross-referencing keys across attributes
- Error message formatting
- Key location in different AnnData attributes (obsm, layers, etc.)

### Formatter Tests (29 tests)
- Initialization with different data types
- Device type detection (CPU, CUDA, MPS)
- Data type identification (numpy, torch, ArrayView)
- Conversion between numpy and torch
- Device placement for tensors
- Handling of sparse matrices and ArrayViews
- Value preservation through conversions

### Fetcher Tests (29 tests)
- Simple data fetching without grouping
- Grouped data fetching (as dict and as list)
- Fetching different data types (X, obsm, layers)
- Torch tensor conversion
- Device placement
- Sparse matrix handling
- Data consistency checks
- Multiple groupby columns

## Fixtures

The test suite includes the following fixtures (defined in `conftest.py`):

- `sample_adata` - A basic AnnData object with various attributes
- `sample_numpy_array` - A simple numpy array for testing
- `sample_torch_tensor` - A simple torch tensor for testing
- `sparse_adata` - An AnnData object with sparse matrices

## Notes

- Some tests are skipped when CUDA is not available (marked with `@pytest.mark.skipif`)
- Tests for MPS (Apple Silicon GPU) will run if available
- All tests use reproducible random data for consistency

