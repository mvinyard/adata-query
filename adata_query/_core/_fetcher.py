# -- import packages: ----------------------------------------------------------
import ABCParse
import autodevice
import anndata
import logging
import torch as _torch
import numpy as np

# -- configure logger: ---------------------------------------------------------
logger = logging.getLogger(__name__)


# -- import local dependencies: ------------------------------------------------
from ._locator import locate
from ._formatter import format_data


# -- set typing: ---------------------------------------------------------------
from typing import Dict, List, Optional, Union


# -- operational class: --------------------------------------------------------
class AnnDataFetcher(ABCParse.ABCParse):
    """Operational class powering the fetch function."""
    def __init__(self, *args, **kwargs) -> None:
        self.__parse__(locals(), public=[None])
        logger.debug("Initialized AnnDataFetcher")

    @property
    def _GROUPED(self):
        logger.debug(f"Grouping data by: {self._groupby}")
        return self._adata.obs.groupby(self._groupby)

    def _forward(self, adata, key):
        logger.debug(f"Fetching data for key: {key}")
        if key == "X":
            data = getattr(adata, "X")
            logger.debug("Retrieved data from adata.X")
        else:
            attr = locate(adata, key)
            data = getattr(adata, attr)[key]
            logger.debug(f"Retrieved data from adata.{attr}['{key}']")
        return format_data(data=data, torch=self._torch, device=self._device)

    def _grouped_subroutine(self, adata, key):
        logger.info(f"Processing grouped data for key: {key}")
        if self._as_dict:
            for group, group_df in self._GROUPED:
                logger.debug(f"Processing group: {group}")
                yield group, self._forward(adata[group_df.index], key)
        else:
            for group, group_df in self._GROUPED:
                logger.debug(f"Processing group: {group}")
                yield self._forward(adata[group_df.index], key)

    def __call__(
        self,
        adata: anndata.AnnData,
        key: str,
        groupby: Optional[str] = None,
        torch: bool = False,
        device: _torch.device = autodevice.AutoDevice(),
        as_dict: bool = True,
    ):
        """
        adata: anndata.AnnData [ required ]
            Annotated single-cell data object.
        
        key: str [ required ]
            Key to access a matrix in adata. For example, if you wanted to access
            adata.obsm['X_pca'], you would pass: "X_pca".
        
        groupby: Optional[str], default = None
            Optionally, one may choose to group data according to a cell-specific
            annotation in adata.obs. This would invoke returning data as List
            
        torch: bool, default = False
            Boolean indicator of whether data should be formatted as torch.Tensor. If
            False (default), data is formatted as np.ndarray.device (torch.device) =
            autodevice.AutoDevice(). Should torch=True, the device ("cpu", "cuda:N", 
            "mps:N") may be set. The default value, autodevice.AutoDevice() will 
            indicate the use of GPU, if available.

        device: torch.device, default = autodevice.AutoDevice()
            
    
        as_dict: bool, default = True
            Only relevant when `groupby` is not None. Boolean indicator to return
            data in a Dict where the key for each value corresponds to the respective
            `groupby` value. If False, returns List.
        """
        logger.info(f"Fetch called for key: {key}" + (f" with groupby: {groupby}" if groupby else ""))
        self.__update__(locals(), public=[None])

        if hasattr(self, "_groupby"):
            logger.info(f"Returning grouped data as {'dictionary' if self._as_dict else 'list'}")
            if self._as_dict:
                return dict(self._grouped_subroutine(adata, key))
            return list(self._grouped_subroutine(adata, key))
        return self._forward(adata, key)

def fetch(
    adata: anndata.AnnData,
    key: str,
    groupby: Optional[str] = None,
    torch: bool = False,
    device: _torch.device = autodevice.AutoDevice(),
    as_dict: bool = True,
    *args,
    **kwargs,
) -> Union[
    _torch.Tensor,
    np.ndarray,
    List[Union[_torch.Tensor, np.ndarray]],
    Dict[Union[str, int], Union[_torch.Tensor, np.ndarray]],
]:
    """
    Given, adata and a key that points to a specific matrix stored in adata,
    return the data, formatted either as np.ndarray or torch.Tensor. If formatted
    as torch.Tensor, device may be specified based on available devices.

    Parameters
    ----------
    adata: anndata.AnnData [ required ]
        Annotated single-cell data object.

    key: str [ required ]
        Key to access a matrix in adata. For example, if you wanted to access
        adata.obsm['X_pca'], you would pass: "X_pca".

    groupby: Optional[str], default = None
        Optionally, one may choose to group data according to a cell-specific
        annotation in adata.obs. This would invoke returning data as List

    torch: bool, default = False
        Boolean indicator of whether data should be formatted as torch.Tensor. If
        False (default), data is formatted as np.ndarray.device (torch.device) =
        autodevice.AutoDevice(). Should torch=True, the device ("cpu", "cuda:N",
        "mps:N") may be set. The default value, autodevice.AutoDevice() will
        indicate the use of GPU, if available.

    as_dict: bool, default = True
        Only relevant when `groupby` is not None. Boolean indicator to return
        data in a Dict where the key for each value corresponds to the respective
        `groupby` value. If False, returns List.

    Returns
    -------
    data: Union[torch.Tensor, np.ndarray, List[Union[torch.Tensor, np.ndarray]], Dict[Union[str, int], Union[torch.Tensor, np.ndarray]]
        Formatted data as np.ndarray or torch.Tensor. If torch=True the torch.Tensor
        is allocated to the device indicated by the device argument. If `groupby` is passed,
        returned as Dict[np.ndarray] or Dict[torch.Tensor]. If groupby is passed and `as_dict`
        = False, returns List[np.ndarray] or List[torch.Tensor].
    """
    logger.info(f"Fetch function called for key: {key}" + (f" with groupby: {groupby}" if groupby else ""))
    fetcher = AnnDataFetcher()
    return fetcher(
        adata=adata,
        key=key,
        groupby=groupby,
        torch=torch,
        device=device,
        as_dict=as_dict,
        *args,
        **kwargs,
    )
