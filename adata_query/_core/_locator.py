
# -- set typing: ---------------------------------------------------------------
from typing import List, Optional


# -- import packages: ----------------------------------------------------------
import ABCParse
import anndata
import numpy as np


# -- operational class: --------------------------------------------------------
class AnnDataLocator(ABCParse.ABCParse):
    """Query available key values of AnnData. Operational class powering the `locate` function."""
    def __init__(self, searchable: Optional[List[str]] = None, *args, **kwargs) -> None:
        """Description.

        Args:
            searchable (Optional[List[str]])
            Description.
            **Default**: None

        Returns:
            None
                Initializes class object.
        """

        self._ATTRS = {}
        self._searchable = ["X"]
        if not searchable is None:
            self._searchable += searchable

    def _stash(self, attr: str, attr_val: np.ndarray) -> None:
        """Description.

        Args:
            attr (str)

            attr_val (np.ndarray)

        Returns:
            None, updates `self._ATTRS` and sets the (attr, attr_val) key, value pair.
        """
        self._ATTRS[attr] = attr_val
        setattr(self, attr, attr_val)

    def _intake(self, adata: anndata.AnnData) -> None:
        """Description.

        Args:
            adata (anndata.AnnData)

        Returns:
            None, sets class attributes.

        """
        for attr in adata.__dir__():
            if "key" in attr:
                attr_val = getattr(adata, attr)()
                self._stash(attr, attr_val)
            if attr == "layers":
                attr_val = list(getattr(adata, attr))
                self._stash(attr, attr_val)
            if attr in self._searchable:
                self._stash(attr, attr)

    def _cross_reference(self, passed_key: str) -> List[str]:
        """Description.

        Args:
            passed_key (str)

        Returns:
            (List[str])
        """
        return [key for key, val in self._ATTRS.items() if passed_key in val]

    def _query_str_vals(self, query_result: List[str]) -> str:
        """Description.

        Args:
            query_result (List[str])

        Returns:
            (str)
        """
        return ", ".join(query_result)

    def _format_error_msg(self, key: str, query_result: List[str]) -> str:
        """Description.

        Args
            key (str)

            query_result (List[str])

        Returns:
            (str)
        """
        if len(query_result) > 1:
            return f"Found more than one match: [{self._query_str_vals(query_result)}]"
        return f"{key} NOT FOUND"

    def _format_output_str(self, query_result: List[str]):
        """Description.

        Args:
            query_result (List[str])

        Returns:

        """
        return query_result[0].split("_keys")[0]

    def _forward(self, adata: anndata.AnnData, key: str) -> str:
        """Description.

        Args:
            adata (anndata.AnnData)

            key (str)

        Returns:
            (str)

        """
        self._intake(adata)
        query_result = self._cross_reference(passed_key=key)

        if len(query_result) != 1:
            raise KeyError(self._format_error_msg(key, query_result))

        return self._format_output_str(query_result)

    def __call__(self, adata: anndata.AnnData, key: str) -> str:
        """Description of __call__ func.

        Args:
            adata (anndata.AnnData)

            key (str)

        Returns:
            attr (str)
        """

        return self._forward(adata, key)


# -- API-facing function: -----------------------------------------------------
def locate(adata: anndata.AnnData, key: str, test_kw: List) -> str:
    
    """Given adata and a key that points to a specific matrix stored in adata,
    return the data, formatted either as np.ndarray or torch.Tensor. If
    formatted as torch.Tensor, device may be specified based on available
    devices.

    Args:
        adata (anndata.AnnData). The (annotated) single-cell data matrix of shape
        n_obs × n_vars. Rows correspond to cells and columns to genes.
        For more: https://anndata.readthedocs.io/en/latest/.

        key (str). Key to access a matrix in adata. For example, if you wanted to
            access adata.obsm['X_pca'], you would pass: "X_pca".
            
        test_kw (List). This is just a test. Please find it.

    Returns:
        attr_key (str). Attribute of adata containing the passed key.
    """

    locator = AnnDataLocator()

    return locator(adata=adata, key=key)
