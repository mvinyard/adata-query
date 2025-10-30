# -- import packages: ---------------------------------------------------------
import numpy as np
import pytest
import torch

# -- import from package: -----------------------------------------------------
from adata_query._core._fetcher import AnnDataFetcher, fetch


# -- test AnnDataFetcher class: -----------------------------------------------
class TestAnnDataFetcher:
    """Test suite for AnnDataFetcher class."""

    def test_init(self):
        """Test AnnDataFetcher initialization."""
        fetcher = AnnDataFetcher()
        assert fetcher is not None

    def test_grouped_property(self, sample_adata):
        """Test _GROUPED property returns pandas GroupBy object."""
        fetcher = AnnDataFetcher()
        fetcher._adata = sample_adata
        fetcher._groupby = 'cell_type'

        grouped = fetcher._GROUPED
        assert grouped is not None
        # Check it's a groupby object
        assert hasattr(grouped, 'groups')

    def test_forward_X(self, sample_adata):
        """Test _forward method fetching X data."""
        fetcher = AnnDataFetcher()
        fetcher._torch = False
        fetcher._device = 'cpu'

        result = fetcher._forward(sample_adata, 'X')
        assert isinstance(result, np.ndarray)
        assert result.shape == sample_adata.X.shape

    def test_forward_obsm(self, sample_adata):
        """Test _forward method fetching obsm data."""
        fetcher = AnnDataFetcher()
        fetcher._torch = False
        fetcher._device = 'cpu'

        result = fetcher._forward(sample_adata, 'X_pca')
        assert isinstance(result, np.ndarray)
        assert result.shape == sample_adata.obsm['X_pca'].shape

    def test_forward_layer(self, sample_adata):
        """Test _forward method fetching layer data."""
        fetcher = AnnDataFetcher()
        fetcher._torch = False
        fetcher._device = 'cpu'

        result = fetcher._forward(sample_adata, 'raw')
        assert isinstance(result, np.ndarray)
        assert result.shape == sample_adata.layers['raw'].shape

    def test_forward_with_torch(self, sample_adata):
        """Test _forward method with torch=True."""
        fetcher = AnnDataFetcher()
        fetcher._torch = True
        fetcher._device = 'cpu'

        result = fetcher._forward(sample_adata, 'X_pca')
        assert isinstance(result, torch.Tensor)

    def test_call_simple(self, sample_adata):
        """Test __call__ method without groupby."""
        fetcher = AnnDataFetcher()
        result = fetcher(sample_adata, 'X_pca', torch=False)

        assert isinstance(result, np.ndarray)
        assert result.shape == sample_adata.obsm['X_pca'].shape

    def test_call_with_groupby_dict(self, sample_adata):
        """Test __call__ method with groupby returning dict."""
        fetcher = AnnDataFetcher()
        result = fetcher(sample_adata, 'X_pca', groupby='cell_type', as_dict=True)

        assert isinstance(result, dict)
        assert 'A' in result
        assert 'B' in result
        assert isinstance(result['A'], np.ndarray)
        assert isinstance(result['B'], np.ndarray)
        # Each group should have 5 observations
        assert result['A'].shape[0] == 5
        assert result['B'].shape[0] == 5

    def test_call_with_groupby_list(self, sample_adata):
        """Test __call__ method with groupby returning list."""
        fetcher = AnnDataFetcher()
        result = fetcher(sample_adata, 'X_pca', groupby='cell_type', as_dict=False)

        assert isinstance(result, list)
        assert len(result) == 2  # Two groups: A and B
        assert all(isinstance(item, np.ndarray) for item in result)

    def test_call_with_torch(self, sample_adata):
        """Test __call__ method with torch=True."""
        fetcher = AnnDataFetcher()
        result = fetcher(sample_adata, 'X_pca', torch=True, device='cpu')

        assert isinstance(result, torch.Tensor)
        assert result.device.type == 'cpu'

    def test_call_with_groupby_and_torch(self, sample_adata):
        """Test __call__ method with both groupby and torch."""
        fetcher = AnnDataFetcher()
        result = fetcher(sample_adata, 'X_pca', groupby='cell_type',
                        torch=True, as_dict=True)

        assert isinstance(result, dict)
        assert 'A' in result
        assert 'B' in result
        assert isinstance(result['A'], torch.Tensor)
        assert isinstance(result['B'], torch.Tensor)

    def test_call_X_data(self, sample_adata):
        """Test __call__ method fetching X matrix."""
        fetcher = AnnDataFetcher()
        result = fetcher(sample_adata, 'X')

        assert isinstance(result, np.ndarray)
        assert result.shape == sample_adata.X.shape

    def test_grouped_subroutine_as_dict(self, sample_adata):
        """Test _grouped_subroutine yields correct format with as_dict=True."""
        fetcher = AnnDataFetcher()
        fetcher._adata = sample_adata
        fetcher._groupby = 'cell_type'
        fetcher._torch = False
        fetcher._device = 'cpu'
        fetcher._as_dict = True

        results = list(fetcher._grouped_subroutine(sample_adata, 'X_pca'))

        assert len(results) == 2
        # Each result should be a tuple of (group_name, data)
        assert all(isinstance(item, tuple) for item in results)
        assert all(len(item) == 2 for item in results)
        assert all(isinstance(item[1], np.ndarray) for item in results)

    def test_grouped_subroutine_as_list(self, sample_adata):
        """Test _grouped_subroutine yields correct format with as_dict=False."""
        fetcher = AnnDataFetcher()
        fetcher._adata = sample_adata
        fetcher._groupby = 'cell_type'
        fetcher._torch = False
        fetcher._device = 'cpu'
        fetcher._as_dict = False

        results = list(fetcher._grouped_subroutine(sample_adata, 'X_pca'))

        assert len(results) == 2
        # Each result should be just the data
        assert all(isinstance(item, np.ndarray) for item in results)


# -- test fetch function: -----------------------------------------------------
class TestFetchFunction:
    """Test suite for fetch function."""

    def test_fetch_simple(self, sample_adata):
        """Test fetch function without groupby."""
        result = fetch(sample_adata, 'X_pca')

        assert isinstance(result, np.ndarray)
        assert result.shape == sample_adata.obsm['X_pca'].shape

    def test_fetch_X(self, sample_adata):
        """Test fetch function for X matrix."""
        result = fetch(sample_adata, 'X')

        assert isinstance(result, np.ndarray)
        assert result.shape == sample_adata.X.shape

    def test_fetch_layer(self, sample_adata):
        """Test fetch function for layer data."""
        result = fetch(sample_adata, 'raw')

        assert isinstance(result, np.ndarray)
        assert result.shape == sample_adata.layers['raw'].shape

    def test_fetch_with_groupby_dict(self, sample_adata):
        """Test fetch function with groupby returning dict."""
        result = fetch(sample_adata, 'X_pca', groupby='cell_type', as_dict=True)

        assert isinstance(result, dict)
        assert 'A' in result
        assert 'B' in result
        assert isinstance(result['A'], np.ndarray)
        assert isinstance(result['B'], np.ndarray)

    def test_fetch_with_groupby_list(self, sample_adata):
        """Test fetch function with groupby returning list."""
        result = fetch(sample_adata, 'X_pca', groupby='cell_type', as_dict=False)

        assert isinstance(result, list)
        assert len(result) == 2
        assert all(isinstance(item, np.ndarray) for item in result)

    def test_fetch_with_torch(self, sample_adata):
        """Test fetch function with torch=True."""
        result = fetch(sample_adata, 'X_pca', torch=True, device='cpu')

        assert isinstance(result, torch.Tensor)
        assert result.device.type == 'cpu'

    @pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
    def test_fetch_with_cuda(self, sample_adata):
        """Test fetch function with CUDA device."""
        result = fetch(sample_adata, 'X_pca', torch=True, device='cuda')

        assert isinstance(result, torch.Tensor)
        assert result.device.type == 'cuda'

    def test_fetch_with_groupby_and_torch(self, sample_adata):
        """Test fetch function with both groupby and torch."""
        result = fetch(sample_adata, 'X_pca', groupby='cell_type',
                      torch=True, as_dict=True)

        assert isinstance(result, dict)
        assert all(isinstance(v, torch.Tensor) for v in result.values())

    def test_fetch_multiple_keys(self, sample_adata):
        """Test fetch function with different keys."""
        result_pca = fetch(sample_adata, 'X_pca')
        result_umap = fetch(sample_adata, 'X_umap')

        assert result_pca.shape[1] == 3  # PCA has 3 components
        assert result_umap.shape[1] == 2  # UMAP has 2 components

    def test_fetch_preserves_grouping(self, sample_adata):
        """Test that fetch preserves group structure correctly."""
        result = fetch(sample_adata, 'X', groupby='cell_type', as_dict=True)

        # Manually verify group sizes
        cell_type_counts = sample_adata.obs['cell_type'].value_counts()
        for group_name, group_data in result.items():
            expected_size = cell_type_counts[group_name]
            assert group_data.shape[0] == expected_size

    def test_fetch_with_batch_grouping(self, sample_adata):
        """Test fetch with different groupby column."""
        result = fetch(sample_adata, 'X_pca', groupby='batch', as_dict=True)

        assert isinstance(result, dict)
        assert 'batch1' in result
        assert 'batch2' in result

    def test_fetch_sparse_data(self, sparse_adata):
        """Test fetch function with sparse matrices."""
        result = fetch(sparse_adata, 'X')

        assert isinstance(result, np.ndarray)
        assert result.shape == sparse_adata.X.shape

    def test_fetch_sparse_with_groupby(self, sparse_adata):
        """Test fetch function with sparse matrices and groupby."""
        result = fetch(sparse_adata, 'X', groupby='cell_type', as_dict=True)

        assert isinstance(result, dict)
        assert all(isinstance(v, np.ndarray) for v in result.values())

    def test_fetch_data_consistency(self, sample_adata):
        """Test that fetch returns consistent data."""
        # Fetch same data twice
        result1 = fetch(sample_adata, 'X_pca')
        result2 = fetch(sample_adata, 'X_pca')

        np.testing.assert_array_equal(result1, result2)

    def test_fetch_grouped_data_consistency(self, sample_adata):
        """Test that grouped fetch returns consistent data."""
        result1 = fetch(sample_adata, 'X_pca', groupby='cell_type', as_dict=True)
        result2 = fetch(sample_adata, 'X_pca', groupby='cell_type', as_dict=True)

        assert result1.keys() == result2.keys()
        for key in result1.keys():
            np.testing.assert_array_equal(result1[key], result2[key])

