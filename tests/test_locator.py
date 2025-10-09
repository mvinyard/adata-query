# -- import packages: ---------------------------------------------------------
import pytest

# -- import from package: -----------------------------------------------------
from adata_query._core._locator import AnnDataLocator, locate


# -- test AnnDataLocator class: -----------------------------------------------
class TestAnnDataLocator:
    """Test suite for AnnDataLocator class."""

    def test_init_default(self):
        """Test AnnDataLocator initialization with default parameters."""
        locator = AnnDataLocator()
        assert hasattr(locator, '_ATTRS')
        assert hasattr(locator, '_searchable')
        assert locator._searchable == ['X']
        assert locator._ATTRS == {}

    def test_init_with_searchable(self):
        """Test AnnDataLocator initialization with custom searchable list."""
        searchable = ['custom_attr']
        locator = AnnDataLocator(searchable=searchable)
        assert 'X' in locator._searchable
        assert 'custom_attr' in locator._searchable
        assert len(locator._searchable) == 2

    def test_stash(self):
        """Test _stash method stores attributes correctly."""
        locator = AnnDataLocator()
        locator._stash('test_attr', 'test_value')
        assert 'test_attr' in locator._ATTRS
        assert locator._ATTRS['test_attr'] == 'test_value'
        assert hasattr(locator, 'test_attr')
        assert locator.test_attr == 'test_value'

    def test_intake(self, sample_adata):
        """Test _intake method processes AnnData object correctly."""
        locator = AnnDataLocator()
        locator._intake(sample_adata)

        # Check that obsm_keys was processed
        assert 'obsm_keys' in locator._ATTRS
        assert isinstance(locator._ATTRS['obsm_keys'], type(sample_adata.obsm_keys()))

        # Check that layers was processed
        assert 'layers' in locator._ATTRS
        assert isinstance(locator._ATTRS['layers'], list)

        # Check that X was processed
        assert 'X' in locator._ATTRS

    def test_cross_reference_found(self, sample_adata):
        """Test _cross_reference finds keys correctly."""
        locator = AnnDataLocator()
        locator._intake(sample_adata)

        # Test finding a key that exists in obsm
        result = locator._cross_reference('X_pca')
        assert len(result) == 1
        assert 'obsm_keys' in result[0]

    def test_cross_reference_not_found(self, sample_adata):
        """Test _cross_reference returns empty list for non-existent key."""
        locator = AnnDataLocator()
        locator._intake(sample_adata)

        result = locator._cross_reference('nonexistent_key')
        assert len(result) == 0

    def test_query_str_vals(self):
        """Test _query_str_vals formats output correctly."""
        locator = AnnDataLocator()
        result = locator._query_str_vals(['attr1', 'attr2', 'attr3'])
        assert result == 'attr1, attr2, attr3'

    def test_format_error_msg_multiple_matches(self):
        """Test _format_error_msg for multiple matches."""
        locator = AnnDataLocator()
        msg = locator._format_error_msg('test_key', ['match1', 'match2'])
        assert 'more than one match' in msg.lower()
        assert 'match1' in msg
        assert 'match2' in msg

    def test_format_error_msg_no_matches(self):
        """Test _format_error_msg for no matches."""
        locator = AnnDataLocator()
        msg = locator._format_error_msg('test_key', [])
        assert 'NOT FOUND' in msg
        assert 'test_key' in msg

    def test_format_output_str(self):
        """Test _format_output_str extracts attribute name correctly."""
        locator = AnnDataLocator()
        result = locator._format_output_str(['obsm_keys'])
        assert result == 'obsm'

        result = locator._format_output_str(['varm_keys'])
        assert result == 'varm'

    def test_forward_success(self, sample_adata):
        """Test _forward method successfully locates a key."""
        locator = AnnDataLocator()
        result = locator._forward(sample_adata, 'X_pca')
        assert result == 'obsm'

        result = locator._forward(sample_adata, 'X_umap')
        assert result == 'obsm'

    def test_forward_key_not_found(self, sample_adata):
        """Test _forward raises KeyError when key is not found."""
        locator = AnnDataLocator()
        with pytest.raises(KeyError) as exc_info:
            locator._forward(sample_adata, 'nonexistent_key')
        assert 'NOT FOUND' in str(exc_info.value)

    def test_call_method(self, sample_adata):
        """Test __call__ method works correctly."""
        locator = AnnDataLocator()
        result = locator(sample_adata, 'X_pca')
        assert result == 'obsm'

    def test_locate_layers(self, sample_adata):
        """Test locating keys in layers."""
        locator = AnnDataLocator()
        result = locator(sample_adata, 'raw')
        assert result == 'layers'

        result = locator(sample_adata, 'normalized')
        assert result == 'layers'


# -- test locate function: ----------------------------------------------------
class TestLocateFunction:
    """Test suite for locate function."""

    def test_locate_obsm_key(self, sample_adata):
        """Test locate function finds keys in obsm."""
        result = locate(sample_adata, 'X_pca')
        assert result == 'obsm'

        result = locate(sample_adata, 'X_umap')
        assert result == 'obsm'

    def test_locate_layer_key(self, sample_adata):
        """Test locate function finds keys in layers."""
        result = locate(sample_adata, 'raw')
        assert result == 'layers'

        result = locate(sample_adata, 'normalized')
        assert result == 'layers'

    def test_locate_not_found(self, sample_adata):
        """Test locate function raises KeyError for non-existent key."""
        with pytest.raises(KeyError) as exc_info:
            locate(sample_adata, 'nonexistent_key')
        assert 'NOT FOUND' in str(exc_info.value)

    def test_locate_with_searchable(self, sample_adata):
        """Test locate function with custom searchable parameter."""
        result = locate(sample_adata, 'X_pca', searchable=None)
        assert result == 'obsm'

    def test_locate_sparse_adata(self, sparse_adata):
        """Test locate function works with sparse matrices."""
        result = locate(sparse_adata, 'X_pca')
        assert result == 'obsm'

