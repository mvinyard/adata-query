
# -- import packages: ----------------------------------------------------------
import ABCParse
import anndata
import numpy as np


# -- set typing: ---------------------------------------------------------------
from typing import List


# -- operational class: --------------------------------------------------------
class AnnDataLocator(ABCParse.ABCParse):
    """Query available key values of AnnData. Operational class powering the `locate` function."""
    def __init__(self, *args, **kwargs) -> None:
        """
        Parameters
        ----------
        
        Returns
        -------
        None, initializes class object.
        """
        
        self._ATTRS = {}

    def _stash(self, attr: str, attr_val: np.ndarray) -> None:
        """
        Parameters
        ----------
        attr: str
        
        attr_val: np.ndarray
        
        Returns
        -------
        None, updates `self._ATTRS` and sets the (attr, attr_val) key, value pair.
        """
        self._ATTRS[attr] = attr_val
        setattr(self, attr, attr_val)

    def _intake(self, adata: anndata.AnnData) -> None:
        """
        Parameters
        ----------
        adata
        
        Returns
        -------
        
        """
        for attr in adata.__dir__():
            if "key" in attr:
                attr_val = getattr(adata, attr)()
                self._stash(attr, attr_val)
            if attr == "layers":
                attr_val = list(getattr(adata, attr))
                self._stash(attr, attr_val)

    def _cross_reference(self, passed_key: str) -> List[str]:
        """
        Parameters
        ----------
        
        Returns
        -------
        
        """
        return [key for key, val in self._ATTRS.items() if passed_key in val]

    def _query_str_vals(self, query_result: List[str]) -> str:
        """
        Parameters
        ----------
        
        Returns
        -------
        
        """
        return ", ".join(query_result)

    def _format_error_msg(self, key: str, query_result: List[str]) -> str:
        """
        Parameters
        ----------
        
        Returns
        -------
        
        """
        if len(query_result) > 1:
            return f"Found more than one match: [{self._query_str_vals(query_result)}]"
        return f"{key} NOT FOUND"

    def _format_output_str(self, query_result: List[str]):
        """
        Parameters
        ----------
        
        Returns
        -------
        
        """
        return query_result[0].split("_keys")[0]

    def _forward(self, adata: anndata.AnnData, key: str) -> str:

        """
        Parameters
        ----------
        
        Returns
        -------
        
        """
        self._intake(adata)
        query_result = self._cross_reference(passed_key=key)

        if len(query_result) != 1:
            raise KeyError(self._format_error_msg(key, query_result))

        return self._format_output_str(query_result)

    def __call__(self, adata: anndata.AnnData, key: str) -> str:
        
        """
        Parameters
        ----------
        adata: anndata.AnnData
        
        key: str
        
        Returns
        -------
        attr: str
        """

        return self._forward(adata, key)

    
def locate(adata: anndata.AnnData, key: str) -> str:
    """
    Given, adata and a key that points to a specific matrix stored in adata,
    return the data, formatted either as np.ndarray or torch.Tensor. If formatted
    as torch.Tensor, device may be specified based on available devices.
    
    Parameters
    ----------
    adata: anndata.AnnData, [ required ]
    
    key: str [ required ]
        Key to access a matrix in adata. For example, if you wanted to access
        adata.obsm['X_pca'], you would pass: "X_pca".
    
    Returns
    -------
    attr_key: str
        Attribute of adata containing the passed key
    """
    locator = AnnDataLocator()
    return locator(adata = adata, key = key)
