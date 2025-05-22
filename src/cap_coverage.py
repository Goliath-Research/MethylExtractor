import numpy as np
import pandas as pd

def cap_coverage(df: pd.DataFrame, target_cov: int = 500, method: str = "mean",
                 q: float = 0.5, verbose: bool = True) -> pd.DataFrame:
    """
    Downsample coverage in DataFrame to target level, adjusting mC and uC counts in place.

    Args:
        df: DataFrame with 'mC' and 'uC' columns.
        target_cov: Target coverage level.
        method: Method to set target_cov if None ("median", "mean", "quantile").
        q: Quantile probability for "quantile" method.
        verbose: Print summary if True.

    Returns:
        DataFrame with capped mC and uC values (modified in place).

    Raises:
        ValueError: If method is invalid or required columns are missing.
    """
    # Extract arrays
    mC = df["mC"].to_numpy(dtype=np.float64)
    uC = df["uC"].to_numpy(dtype=np.float64)
    cov = mC + uC

    # Set target coverage if None
    target_cov = (
        np.nanmedian(cov) if method == "median" else
        np.nanmean(cov) if method == "mean" else
        np.nanquantile(cov, q, method="linear")
    )

    # Adjust high-coverage sites in place
    mask = cov > target_cov
    if np.any(mask):
        prop = mC[mask] / cov[mask]
        df.loc[mask, "mC"] = np.round(target_cov * prop).astype(int)
        df.loc[mask, "uC"] = target_cov - df.loc[mask, "mC"]

    return df