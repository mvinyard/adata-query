
# -- set typing: --------------------------------------------------------------
from typing import Dict, List, Optional, Union


# -- import packages: ---------------------------------------------------------
import ABCParse
import autodevice
import anndata
import torch as _torch
import pandas as pd
import numpy as np


# -- import local dependencies: -----------------------------------------------
from ._locator import locate
from ._formatter import format_data


# -- operational class: -------------------------------------------------------
class AnnDataFetcher(ABCParse.ABCParse):
    def __init__(self, *args, **kwargs):
        self.__parse__(locals(), public=[None])

    @property
    def _GROUPED(self) -> pd.core.groupby.DataFrameGroupBy:
        return self._adata.obs.groupby(self._groupby)

    def _forward(self, adata: anndata.AnnData, key: str) -> np.ndarray:
        if key == "X":
            data = getattr(adata, "X")
        else:
            data = getattr(adata, locate(adata, key))[key]
        return format_data(data=data, torch=self._torch, device=self._device)

    def _grouped_subroutine(
        self,
        adata: anndata.AnnData,
        key: str,
    ) -> Union[List, Dict[str, np.ndarray]]:
        if self._as_dict:
            for group, group_df in self._GROUPED:
                yield group, self._forward(adata[group_df.index], key)
        else:
            for group, group_df in self._GROUPED:
                yield self._forward(adata[group_df.index], key)

    def __call__(
        self,
        adata: anndata.AnnData,
        key: str,
        groupby: Optional[str] = None,
        torch: bool = False,
        device: _torch.device = autodevice.AutoDevice(),
        as_dict: bool = True,
    ) -> Union[List, Dict[str, np.ndarray]]:
        self.__update__(locals(), public=[None])

        if hasattr(self, "_groupby"):
            if self._as_dict:
                return dict(self._grouped_subroutine(adata, key))
            return list(self._grouped_subroutine(adata, key))
        return self._forward(adata, key)


# -- API-facing function: -----------------------------------------------------
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
