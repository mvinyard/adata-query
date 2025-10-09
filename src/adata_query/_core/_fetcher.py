# -- import packages: ---------------------------------------------------------
import logging
from collections.abc import Generator

# -- set typing: --------------------------------------------------------------
import anndata
import autodevice
import numpy as np
import torch as _torch

from ._formatter import format_data

# -- import local dependencies: -----------------------------------------------
from ._locator import locate

# -- configure logger: --------------------------------------------------------
logger = logging.getLogger(__name__)


# -- operational class: -------------------------------------------------------
class AnnDataFetcher:
    """Fetches and formats data from AnnData objects.
    
    This class provides functionality to retrieve data from AnnData objects, with options
    for grouping, tensor conversion, and device placement. It handles both direct data
    access and grouped data retrieval.
    
    Attributes:
        _GROUPED: Property that returns a pandas GroupBy object when groupby is specified.
    """

    def __init__(self, *args, **kwargs) -> None:
        """Initialize the AnnDataFetcher."""

        logger.debug("Initialized AnnDataFetcher")

    @property
    def _GROUPED(self):
        """Get the grouped data from AnnData.obs.
        
        Returns:
            pandas GroupBy object for the specified groupby column.
        """
        logger.debug(f"Grouping data by: {self._groupby}")
        return self._adata.obs.groupby(self._groupby)

    def _forward(self, adata: anndata.AnnData, key: str) -> _torch.Tensor | np.ndarray:
        """Retrieve and format data for a single key.
        
        Args:
            adata: AnnData object to fetch data from.
            key: Key to fetch data for.
            
        Returns:
            Formatted data as either torch.Tensor or np.ndarray.
        """
        logger.debug(f"Fetching data for key: {key}")
        if key == "X":
            data = adata.X
            logger.debug("Retrieved data from adata.X")
        else:
            attr = locate(adata, key)
            data = getattr(adata, attr)[key]
            logger.debug(f"Retrieved data from adata.{attr}['{key}']")
        return format_data(data=data, torch=self._torch, device=self._device)

    def _grouped_subroutine(
        self,
        adata: anndata.AnnData,
        key: str
    ) -> Generator[tuple[str, _torch.Tensor | np.ndarray] | _torch.Tensor | np.ndarray, None, None]:
        """Process data for each group when groupby is specified.
        
        Args:
            adata: AnnData object to fetch data from.
            key: Key to fetch data for.
            
        Yields:
            If as_dict is True: Tuples of (group_name, formatted_data)
            If as_dict is False: Formatted data for each group
        """
        logger.debug(f"Processing grouped data for key: {key}")
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
        groupby: str | None = None,
        torch: bool = False,
        device: _torch.device = autodevice.AutoDevice(),
        as_dict: bool = True,
    ) -> _torch.Tensor | np.ndarray | list[_torch.Tensor | np.ndarray] | dict[str, _torch.Tensor | np.ndarray]:
        """Fetch and format data from an AnnData object.
        
        Args:
            adata: AnnData object to fetch data from.
            key: Key to fetch data for (e.g., "X_pca" for adata.obsm['X_pca']).
            groupby: Optional column name in adata.obs to group data by.
            torch: Whether to return data as torch.Tensor (True) or np.ndarray (False).
            device: Device to place tensor on if torch=True.
            as_dict: When groupby is specified, whether to return data as a dictionary
                    with group names as keys (True) or as a list (False).
                    
        Returns:
            If groupby is None:
                Single array/tensor for the specified key
            If groupby is specified and as_dict is True:
                Dictionary mapping group names to arrays/tensors
            If groupby is specified and as_dict is False:
                List of arrays/tensors for each group
                
        Example:
            >>> import anndata
            >>> adata = anndata.AnnData(X=[[1, 2], [3, 4]])
            >>> adata.obs['cell_type'] = ['A', 'B']
            >>> adata.obsm['X_pca'] = [[0.1, 0.2], [0.3, 0.4]]
            >>> fetcher = AnnDataFetcher()
            >>> # Get single array
            >>> data = fetcher(adata, "X_pca")
            >>> # Get grouped data as dictionary
            >>> grouped_data = fetcher(adata, "X_pca", groupby="cell_type")
        """
        logger.debug(
            f"Fetch called for key: {key}"
            + (f" with groupby: {groupby}" if groupby else "")
        )

        self._adata = adata
        self._key = key
        self._groupby = groupby
        self._torch = torch
        self._device = device
        self._as_dict = as_dict

        if self._groupby is not None:
            logger.debug(
                f"Returning grouped data as {'dictionary' if self._as_dict else 'list'}"
            )
            if self._as_dict:
                return dict(self._grouped_subroutine(self._adata, self._key))
            return list(self._grouped_subroutine(self._adata, self._key))
        return self._forward(self._adata, self._key)


def fetch(
    adata: anndata.AnnData,
    key: str,
    groupby: str | None = None,
    torch: bool = False,
    device: _torch.device = autodevice.AutoDevice(),
    as_dict: bool = True,
    *args,
    **kwargs,
) -> _torch.Tensor | np.ndarray | list[_torch.Tensor | np.ndarray] | dict[str, _torch.Tensor | np.ndarray]:
    """Fetch and format data from an AnnData object.
    
    This function provides a convenient interface to retrieve and format data from AnnData
    objects. It supports both direct data access and grouped data retrieval, with options
    for tensor conversion and device placement.
    
    Args:
        adata: AnnData object to fetch data from.
        key: Key to fetch data for (e.g., "X_pca" for adata.obsm['X_pca']).
        groupby: Optional column name in adata.obs to group data by.
        torch: Whether to return data as torch.Tensor (True) or np.ndarray (False).
        device: Device to place tensor on if torch=True.
        as_dict: When groupby is specified, whether to return data as a dictionary
                with group names as keys (True) or as a list (False).
                
    Returns:
        If groupby is None:
            Single array/tensor for the specified key
        If groupby is specified and as_dict is True:
            Dictionary mapping group names to arrays/tensors
        If groupby is specified and as_dict is False:
            List of arrays/tensors for each group
            
    Example:
        >>> import anndata
        >>> adata = anndata.AnnData(X=[[1, 2], [3, 4]])
        >>> adata.obs['cell_type'] = ['A', 'B']
        >>> adata.obsm['X_pca'] = [[0.1, 0.2], [0.3, 0.4]]
        >>> # Get single array
        >>> data = fetch(adata, "X_pca")
        >>> # Get grouped data as dictionary
        >>> grouped_data = fetch(adata, "X_pca", groupby="cell_type")
        >>> # Get tensor on GPU
        >>> tensor_data = fetch(adata, "X_pca", torch=True)
    """
    logger.debug(
        f"Fetch function called for key: {key}"
        + (f" with groupby: {groupby}" if groupby else "")
    )
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
