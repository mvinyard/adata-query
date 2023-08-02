
# -- import packages: ----------------------------------------------------------
import ABCParse
import autodevice
import anndata
import torch as _torch


# -- import local dependencies: ------------------------------------------------
from ._locator import locate
from ._formatter import format_data


# -- set typing: ---------------------------------------------------------------
from typing import Optional


# -- operational class: --------------------------------------------------------
class AnnDataFetcher(ABCParse.ABCParse):
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
            yield self._forward(adata[group_df.index], key)

    def __call__(
        self,
        adata: anndata.AnnData,
        key: str,
        groupby: Optional[str] = None,
        torch: bool = False,
        device: _torch.device = autodevice.AutoDevice(),
    ):
        """
        adata: anndata.AnnData [required]
        
        key: str [required]
        
        groupby: Optional[str], default = None
        
        torch: bool, default = False
        
        device: torch.device, default = autodevice.AutoDevice()
        """

        self.__update__(locals(), public=[None])

        if hasattr(self, "_groupby"):
            return list(self._grouped_subroutine(adata, key))
        return self._forward(adata, key)

def fetch(
    adata: anndata.AnnData,
    key: str,
    groupby: Optional[str] = None,
    torch: bool = False,
    device: _torch.device = autodevice.AutoDevice(),
    *args,
    **kwargs,
    ):
    
    fetcher = AnnDataFetcher()
    
    return fetcher(
        adata = adata,
        key = key,
        groupby = groupby,
        torch = torch,
        device = device,
        *args,
        **kwargs,
    )