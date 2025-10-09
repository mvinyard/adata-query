# -- import packages: ---------------------------------------------------------
import anndata
import numpy as np
import pytest
import scipy.sparse
import torch


# -- fixtures: ----------------------------------------------------------------
@pytest.fixture
def sample_adata():
    """Create a sample AnnData object for testing.
    
    Returns:
        AnnData object with various attributes populated.
    """
    n_obs = 10
    n_vars = 5

    # Create basic AnnData object
    X = np.random.rand(n_obs, n_vars)
    adata = anndata.AnnData(X=X)

    # Add observations metadata
    adata.obs['cell_type'] = ['A'] * 5 + ['B'] * 5
    adata.obs['batch'] = ['batch1', 'batch2'] * 5

    # Add obsm (multidimensional observations)
    adata.obsm['X_pca'] = np.random.rand(n_obs, 3)
    adata.obsm['X_umap'] = np.random.rand(n_obs, 2)

    # Add layers
    adata.layers['raw'] = np.random.rand(n_obs, n_vars)
    adata.layers['normalized'] = np.random.rand(n_obs, n_vars)

    # Add variables metadata
    adata.var['gene_name'] = [f'gene_{i}' for i in range(n_vars)]

    return adata


@pytest.fixture
def sample_numpy_array():
    """Create a sample numpy array for testing.
    
    Returns:
        numpy array with shape (5, 3).
    """
    return np.array([[1.0, 2.0, 3.0],
                     [4.0, 5.0, 6.0],
                     [7.0, 8.0, 9.0],
                     [10.0, 11.0, 12.0],
                     [13.0, 14.0, 15.0]])


@pytest.fixture
def sample_torch_tensor():
    """Create a sample torch tensor for testing.
    
    Returns:
        torch tensor with shape (5, 3).
    """
    return torch.tensor([[1.0, 2.0, 3.0],
                        [4.0, 5.0, 6.0],
                        [7.0, 8.0, 9.0],
                        [10.0, 11.0, 12.0],
                        [13.0, 14.0, 15.0]])


@pytest.fixture
def sparse_adata():
    """Create an AnnData object with sparse matrices for testing.
    
    Returns:
        AnnData object with sparse X matrix.
    """
    n_obs = 10
    n_vars = 5

    # Create sparse matrix
    X = scipy.sparse.random(n_obs, n_vars, density=0.3, format='csr')
    adata = anndata.AnnData(X=X)

    # Add observations metadata
    adata.obs['cell_type'] = ['A'] * 5 + ['B'] * 5

    # Add obsm
    adata.obsm['X_pca'] = np.random.rand(n_obs, 3)

    return adata

