# -- import packages: ----------------------------------------------------------
import ABCParse
import autodevice
import anndata
import logging
import numpy as np
import torch as _torch


# -- set typing: ---------------------------------------------------------------
from typing import Union

# -- configure logger: ---------------------------------------------------------
logger = logging.getLogger(__name__)


# -- operational class: --------------------------------------------------------
class DataFormatter(ABCParse.ABCParse):
    """Format data to interface with numpy or torch, on a specified device."""
    def __init__(self, data: Union[_torch.Tensor, np.ndarray], *args, **kwargs):
        self.__parse__(locals())
        logger.debug(f"Initialized DataFormatter with data type: {type(data)}")

    @property
    def device_type(self) -> str:
        """Returns device type"""
        if hasattr(self._data, "device"):
            return self._data.device.type
        return "cpu"

    @property
    def is_ArrayView(self) -> bool:
        """Checks if device is of type ArrayView"""
        return isinstance(self._data, anndata._core.views.ArrayView)

    @property
    def is_numpy_array(self) -> bool:
        """Checks if device is of type np.ndarray"""
        return isinstance(self._data, np.ndarray)

    @property
    def is_torch_Tensor(self) -> bool:
        """Checks if device is of type torch.Tensor"""
        return isinstance(self._data, _torch.Tensor)

    @property
    def on_cpu(self) -> bool:
        """Checks if device is on cuda or mps"""
        return self.device_type == "cpu"

    @property
    def on_gpu(self) -> bool:
        """Checks if device is on cuda or mps"""
        return self.device_type in ["cuda", "mps"]

    def to_numpy(self) -> np.ndarray:
        """Sends data to np.ndarray"""
        logger.debug("Converting data to numpy array")
        if self.is_torch_Tensor:
            if self.on_gpu:
                logger.debug("Converting GPU tensor to numpy array")
                return self._data.detach().cpu().numpy()
            logger.debug("Converting CPU tensor to numpy array")
            return self._data.numpy()
        elif self.is_ArrayView:
            logger.debug("Converting ArrayView to numpy array")
            return self._data.toarray()
        logger.debug("Data already in numpy format")
        return self._data

    def to_torch(self, device=autodevice.AutoDevice()) -> _torch.Tensor:
        """Convert data to torch tensor on specified device"""
        logger.debug(f"Converting data to torch tensor on device: {device}")
        self.__update__(locals())

        if self.is_torch_Tensor:
            logger.debug(f"Moving existing tensor to device: {device}")
            return self._data.to(self._device)
        elif self.is_ArrayView:
            logger.debug("Converting ArrayView to numpy before torch conversion")
            self._data = self._data.toarray()
        logger.debug("Converting numpy array to torch tensor")
        return _torch.Tensor(self._data).to(self._device)


# -- functional wrap: ----------------------------------------------------------
def format_data(
    data: Union[np.ndarray, _torch.Tensor], 
    torch: bool = False, 
    device: _torch.device = autodevice.AutoDevice(),
):
    logger.info(f"Formatting data as {'torch tensor' if torch else 'numpy array'}" + 
                (f" on device: {device}" if torch else ""))
    formatter = DataFormatter(data=data)
    if torch:
        return formatter.to_torch(device=device)
    return formatter.to_numpy()
