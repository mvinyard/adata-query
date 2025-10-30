# -- import packages: ---------------------------------------------------------
import numpy as np
import pytest
import torch

# -- import from package: -----------------------------------------------------
from adata_query._core._formatter import DataFormatter, format_data


# -- test DataFormatter class: ------------------------------------------------
class TestDataFormatter:
    """Test suite for DataFormatter class."""

    def test_init_numpy(self, sample_numpy_array):
        """Test DataFormatter initialization with numpy array."""
        formatter = DataFormatter(sample_numpy_array)
        assert hasattr(formatter, '_data')
        assert isinstance(formatter._data, np.ndarray)

    def test_init_torch(self, sample_torch_tensor):
        """Test DataFormatter initialization with torch tensor."""
        formatter = DataFormatter(sample_torch_tensor)
        assert hasattr(formatter, '_data')
        assert isinstance(formatter._data, torch.Tensor)

    def test_device_type_cpu_numpy(self, sample_numpy_array):
        """Test device_type property returns 'cpu' for numpy arrays."""
        formatter = DataFormatter(sample_numpy_array)
        assert formatter.device_type == 'cpu'

    def test_device_type_cpu_torch(self):
        """Test device_type property returns 'cpu' for CPU tensors."""
        tensor = torch.tensor([1, 2, 3])
        formatter = DataFormatter(tensor)
        assert formatter.device_type == 'cpu'

    @pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
    def test_device_type_cuda(self):
        """Test device_type property returns 'cuda' for CUDA tensors."""
        tensor = torch.tensor([1, 2, 3]).cuda()
        formatter = DataFormatter(tensor)
        assert formatter.device_type == 'cuda'

    @pytest.mark.skipif(not torch.backends.mps.is_available(), reason="MPS not available")
    def test_device_type_mps(self):
        """Test device_type property returns 'mps' for MPS tensors."""
        tensor = torch.tensor([1, 2, 3]).to('mps')
        formatter = DataFormatter(tensor)
        assert formatter.device_type == 'mps'

    def test_is_numpy_array(self, sample_numpy_array):
        """Test is_numpy_array property correctly identifies numpy arrays."""
        formatter = DataFormatter(sample_numpy_array)
        assert formatter.is_numpy_array is True
        assert formatter.is_torch_Tensor is False

    def test_is_torch_tensor(self, sample_torch_tensor):
        """Test is_torch_Tensor property correctly identifies torch tensors."""
        formatter = DataFormatter(sample_torch_tensor)
        assert formatter.is_torch_Tensor is True
        assert formatter.is_numpy_array is False

    def test_is_ArrayView(self, sample_adata):
        """Test is_ArrayView property correctly identifies ArrayViews."""
        array_view = sample_adata.obsm['X_pca']
        formatter = DataFormatter(array_view)
        # ArrayView should be detected
        if formatter.is_ArrayView:
            assert formatter.is_ArrayView is True

    def test_on_cpu_numpy(self, sample_numpy_array):
        """Test on_cpu property returns True for numpy arrays."""
        formatter = DataFormatter(sample_numpy_array)
        assert formatter.on_cpu is True
        assert formatter.on_gpu is False

    def test_on_cpu_torch(self):
        """Test on_cpu property returns True for CPU tensors."""
        tensor = torch.tensor([1, 2, 3])
        formatter = DataFormatter(tensor)
        assert formatter.on_cpu is True
        assert formatter.on_gpu is False

    @pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
    def test_on_gpu_cuda(self):
        """Test on_gpu property returns True for CUDA tensors."""
        tensor = torch.tensor([1, 2, 3]).cuda()
        formatter = DataFormatter(tensor)
        assert formatter.on_gpu is True
        assert formatter.on_cpu is False

    def test_to_numpy_from_numpy(self, sample_numpy_array):
        """Test to_numpy with numpy array input."""
        formatter = DataFormatter(sample_numpy_array)
        result = formatter.to_numpy()
        assert isinstance(result, np.ndarray)
        np.testing.assert_array_equal(result, sample_numpy_array)

    def test_to_numpy_from_torch_cpu(self, sample_torch_tensor):
        """Test to_numpy with CPU torch tensor input."""
        formatter = DataFormatter(sample_torch_tensor)
        result = formatter.to_numpy()
        assert isinstance(result, np.ndarray)
        np.testing.assert_array_almost_equal(result, sample_torch_tensor.numpy())

    @pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
    def test_to_numpy_from_torch_gpu(self):
        """Test to_numpy with GPU torch tensor input."""
        tensor = torch.tensor([[1.0, 2.0], [3.0, 4.0]]).cuda()
        formatter = DataFormatter(tensor)
        result = formatter.to_numpy()
        assert isinstance(result, np.ndarray)
        np.testing.assert_array_almost_equal(result, tensor.cpu().numpy())

    def test_to_torch_from_numpy(self, sample_numpy_array):
        """Test to_torch with numpy array input."""
        formatter = DataFormatter(sample_numpy_array)
        result = formatter.to_torch(device='cpu')
        assert isinstance(result, torch.Tensor)
        np.testing.assert_array_almost_equal(result.numpy(), sample_numpy_array)

    def test_to_torch_from_torch(self, sample_torch_tensor):
        """Test to_torch with torch tensor input."""
        formatter = DataFormatter(sample_torch_tensor)
        result = formatter.to_torch(device='cpu')
        assert isinstance(result, torch.Tensor)
        torch.testing.assert_close(result, sample_torch_tensor)

    def test_to_torch_default_device(self, sample_numpy_array):
        """Test to_torch with default device (AutoDevice)."""
        formatter = DataFormatter(sample_numpy_array)
        result = formatter.to_torch(device=None)
        assert isinstance(result, torch.Tensor)

    @pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
    def test_to_torch_cuda_device(self, sample_numpy_array):
        """Test to_torch with CUDA device."""
        formatter = DataFormatter(sample_numpy_array)
        result = formatter.to_torch(device='cuda')
        assert isinstance(result, torch.Tensor)
        assert result.device.type == 'cuda'

    def test_to_torch_from_array_view(self, sample_adata):
        """Test to_torch with ArrayView input."""
        array_view = sample_adata.obsm['X_pca']
        formatter = DataFormatter(array_view)
        result = formatter.to_torch(device='cpu')
        assert isinstance(result, torch.Tensor)


# -- test format_data function: -----------------------------------------------
class TestFormatDataFunction:
    """Test suite for format_data function."""

    def test_format_data_numpy_to_numpy(self, sample_numpy_array):
        """Test format_data returning numpy array from numpy input."""
        result = format_data(sample_numpy_array, torch=False)
        assert isinstance(result, np.ndarray)
        np.testing.assert_array_equal(result, sample_numpy_array)

    def test_format_data_numpy_to_torch(self, sample_numpy_array):
        """Test format_data returning torch tensor from numpy input."""
        result = format_data(sample_numpy_array, torch=True, device='cpu')
        assert isinstance(result, torch.Tensor)
        np.testing.assert_array_almost_equal(result.numpy(), sample_numpy_array)

    def test_format_data_torch_to_numpy(self, sample_torch_tensor):
        """Test format_data returning numpy array from torch input."""
        result = format_data(sample_torch_tensor, torch=False)
        assert isinstance(result, np.ndarray)
        np.testing.assert_array_almost_equal(result, sample_torch_tensor.numpy())

    def test_format_data_torch_to_torch(self, sample_torch_tensor):
        """Test format_data returning torch tensor from torch input."""
        result = format_data(sample_torch_tensor, torch=True, device='cpu')
        assert isinstance(result, torch.Tensor)
        torch.testing.assert_close(result, sample_torch_tensor)

    def test_format_data_default_device(self, sample_numpy_array):
        """Test format_data with default device parameter."""
        result = format_data(sample_numpy_array, torch=True)
        assert isinstance(result, torch.Tensor)

    @pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
    def test_format_data_cuda_device(self, sample_numpy_array):
        """Test format_data with CUDA device."""
        result = format_data(sample_numpy_array, torch=True, device='cuda')
        assert isinstance(result, torch.Tensor)
        assert result.device.type == 'cuda'

    def test_format_data_preserves_values(self):
        """Test that format_data preserves data values through conversions."""
        original = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])

        # numpy -> torch -> numpy
        tensor = format_data(original, torch=True, device='cpu')
        back_to_numpy = format_data(tensor, torch=False)

        np.testing.assert_array_almost_equal(original, back_to_numpy)

    def test_format_data_with_array_view(self, sample_adata):
        """Test format_data with AnnData ArrayView."""
        array_view = sample_adata.obsm['X_pca']

        # Convert to numpy
        result_numpy = format_data(array_view, torch=False)
        assert isinstance(result_numpy, np.ndarray)

        # Convert to torch
        result_torch = format_data(array_view, torch=True, device='cpu')
        assert isinstance(result_torch, torch.Tensor)

