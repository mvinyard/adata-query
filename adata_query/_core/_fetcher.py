
# -- import packages: ----------------------------------------------------------
import ABCParse
import autodevice
import anndata
import torch as _torch


# -- import local dependencies: ------------------------------------------------
from ._locator import locate
from ._formatter import format_data


# -- set typing: ---------------------------------------------------------------
from typing import Optional, Dict, List


# -- operational class: --------------------------------------------------------
class AnnDataFetcher(ABCParse.ABCParse):
    """Operational class powering the fetch function."""
    def __init__(self, *args, **kwargs):

        self.__parse__(locals(), public=[None])

    @property
    def _GROUPED(self):
        return self._adata.obs.groupby(self._groupby)

    def _forward(self, adata, key):
        data = getattr(adata, locate(adata, key))[key]
        return format_data(data=data, torch = self._torch, device = self._device)

    def _grouped_subroutine(self, adata, key):
        for group, group_df in self._GROUPED:
            if self._as_dict:
                yield key, self._forward(adata[group_df.index], key)
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

        self.__update__(locals(), public=[None])

        if hasattr(self, "_groupby"):
            if self._as_dict:
                return dict(self._grouped_subgroutine(adata, key))
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
) -> Union[torch.Tensor, np.ndarray, List[torch.Tensor], List[np.ndarray], Dict[torch.Tensor], Dict[np.ndarray]]:
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
    data: Union[torch.Tensor, np.ndarray, List[torch.Tensor], List[np.ndarray], Dict[torch.Tensor], Dict[np.ndarray]]
        Formatted data as np.ndarray or torch.Tensor. If torch=True the torch.Tensor
        is allocated to the device indicated by the device argument. If groupby is passed,
        returned as a List[np.ndarray] or List[torch.Tensor]
    """
    
    
    fetcher = AnnDataFetcher()
    
    return fetcher(
        adata = adata,
        key = key,
        groupby = groupby,
        torch = torch,
        device = device,
        as_dict = as_dict,
        *args,
        **kwargs,
    )
