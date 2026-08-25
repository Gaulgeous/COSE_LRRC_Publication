import pandas as pd
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import kendalltau, spearmanr, multivariate_normal, norm, chi2, kruskal, ranksums, bootstrap, permutation_test
from factor_analyzer.factor_analyzer import calculate_kmo, calculate_bartlett_sphericity
from factor_analyzer import FactorAnalyzer
from factor_analyzer.rotator import Rotator
from sklearn.decomposition import PCA
from itertools import combinations
from statsmodels.stats.correlation_tools import corr_nearest
import matplotlib.pyplot as plt
import warnings
from IPython.display import display
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
import pingouin as pg


# Suppress minor convergence warnings during bootstrap iterations
warnings.filterwarnings("ignore")

def evaluate_kmo_bartlett_factorability(
    df: pd.DataFrame = None,
    corr_matrix: pd.DataFrame = None,
    n_bootstraps: int = 1000,
    ci_level: float = 0.95,
    random_state: int = 42,
) -> dict:
    """Computes overall and per-item KMO and Bartlett's Test of Sphericity

    using factor_analyzer, statsmodels, and scipy.stats.
    """
    rng = np.random.default_rng(random_state)

    if df is not None:
        clean_df = df.apply(pd.to_numeric, errors="coerce").dropna().astype(float)
        cols = list(clean_df.columns)
        n_obs, n_vars = clean_df.shape
        data_arr = clean_df.to_numpy()
        input_data = clean_df if corr_matrix is None else corr_matrix.copy()
        is_corr = corr_matrix is not None
    elif corr_matrix is not None:
        input_data = corr_matrix.copy()
        cols = list(corr_matrix.columns)
        n_obs = 27
        n_vars = len(cols)
        data_arr = None
        is_corr = True
    else:
        raise ValueError("Must provide either `df` or `corr_matrix`.")

    # 1. Ensure Positive Semi-Definiteness if correlation matrix is provided
    if is_corr:
        evals = np.linalg.eigvalsh(input_data.values)
        if np.any(evals <= 1e-6):
            input_data = pd.DataFrame(
                corr_nearest(input_data.values, threshold=1e-6),
                index=cols,
                columns=cols,
            )

    # 2. Point Estimates via factor_analyzer
    kmo_per_item, kmo_total = calculate_kmo(input_data)
    chi_sq, p_val = calculate_bartlett_sphericity(input_data)

    # 3. Bootstrap Confidence Intervals (when raw data is available)
    if data_arr is not None:
        def _boot_kmo_total_stat(sample_rows, axis=0):
            try:
                _, tot = calculate_kmo(sample_rows)
                return np.atleast_1d(tot)
            except Exception:
                return np.nan

        def _boot_kmo_items_stat(sample_rows, axis=0):
            try:
                items, _ = calculate_kmo(sample_rows)
                return np.atleast_1d(items)
            except Exception:
                return np.full(n_vars, np.nan)

        def _boot_bartlett_chi(sample_rows, axis=0):
            try:
                chi, _ = calculate_bartlett_sphericity(sample_rows)
                return np.atleast_1d(chi)
            except Exception:
                return np.nan

        def _boot_bartlett_pval(sample_rows, axis=0):
            try:
                _, p = calculate_bartlett_sphericity(sample_rows)
                return np.atleast_1d(p)
            except Exception:
                return np.nan

        # SciPy Vectorized Bootstrapping
        kmo_tot_boot = bootstrap(
            (data_arr,),
            statistic=_boot_kmo_total_stat,
            confidence_level=ci_level,
            n_resamples=n_bootstraps,
            method="percentile",
            random_state=rng,
        )

        kmo_items_boot = bootstrap(
            (data_arr,),
            statistic=_boot_kmo_items_stat,
            confidence_level=ci_level,
            n_resamples=n_bootstraps,
            method="percentile",
            random_state=rng,
        )

        chi_boot = bootstrap(
            (data_arr,),
            statistic=_boot_bartlett_chi,
            confidence_level=ci_level,
            n_resamples=n_bootstraps,
            method="percentile",
            random_state=rng,
        )

        pval_boot = bootstrap(
            (data_arr,),
            statistic=_boot_bartlett_pval,
            confidence_level=ci_level,
            n_resamples=n_bootstraps,
            method="percentile",
            random_state=rng,
        )

        kmo_tot_ci = (float(kmo_tot_boot.confidence_interval.low), float(kmo_tot_boot.confidence_interval.high))
        chi_sq_ci = (float(chi_boot.confidence_interval.low), float(chi_boot.confidence_interval.high))
        p_val_ci = (float(pval_boot.confidence_interval.low), float(pval_boot.confidence_interval.high))

        kmo_items_ci_df = pd.DataFrame(
            {
                "KMO_Estimate": kmo_per_item,
                "CI_Lower": kmo_items_boot.confidence_interval.low,
                "CI_Upper": kmo_items_boot.confidence_interval.high,
            },
            index=cols,
        )
    else:
        kmo_tot_ci = (np.nan, np.nan)
        chi_sq_ci = (np.nan, np.nan)
        p_val_ci = (np.nan, np.nan)
        kmo_items_ci_df = pd.DataFrame(
            {
                "KMO_Estimate": kmo_per_item,
                "CI_Lower": np.nan,
                "CI_Upper": np.nan,
            },
            index=cols,
        )

    return {
        "kmo_total": float(kmo_total),
        "kmo_total_ci": kmo_tot_ci,
        "kmo_per_item": pd.Series(kmo_per_item, index=cols, name="KMO"),
        "kmo_per_item_cis": kmo_items_ci_df,
        "bartlett_chi_sq": float(chi_sq),
        "bartlett_chi_sq_ci": chi_sq_ci,
        "bartlett_p_val": float(p_val),
        "bartlett_p_val_ci": p_val_ci,
    }


def _kendall_stat(x, y):
    return kendalltau(x, y, variant="b").statistic


def _spearman_stat(x, y):
    return spearmanr(x, y).statistic


def compute_ordinal_correlations(
    df: pd.DataFrame,
    n_bootstraps: int = 2000,
    confidence_level: float = 0.95,
    random_state: int = 42,
) -> dict:
    """Computes Kendall's tau-b and Spearman's rho correlation matrices, permutation

    p-values, and bootstrapped confidence intervals using scipy.stats.
    """
    cols = list(df.columns)
    n_cols = len(cols)
    rng = np.random.default_rng(random_state)

    # Initialize symmetric DataFrames
    results = {
        "kendall_tau": pd.DataFrame(np.eye(n_cols), index=cols, columns=cols),
        "kendall_pvalues": pd.DataFrame(
            np.zeros((n_cols, n_cols)), index=cols, columns=cols
        ),
        "kendall_ci_lower": pd.DataFrame(
            np.eye(n_cols), index=cols, columns=cols
        ),
        "kendall_ci_upper": pd.DataFrame(
            np.eye(n_cols), index=cols, columns=cols
        ),
        "spearman_rho": pd.DataFrame(np.eye(n_cols), index=cols, columns=cols),
        "spearman_pvalues": pd.DataFrame(
            np.zeros((n_cols, n_cols)), index=cols, columns=cols
        ),
        "spearman_ci_lower": pd.DataFrame(
            np.eye(n_cols), index=cols, columns=cols
        ),
        "spearman_ci_upper": pd.DataFrame(
            np.eye(n_cols), index=cols, columns=cols
        ),
    }

    for col_a_name, col_b_name in combinations(cols, 2):
        col_a = df[col_a_name].to_numpy(dtype=float)
        col_b = df[col_b_name].to_numpy(dtype=float)

        # 1. Point Estimates
        tau_val = _kendall_stat(col_a, col_b)
        rho_val = _spearman_stat(col_a, col_b)

        # 2. Permutation P-Values (Exact/Monte Carlo approach for small N)
        tau_perm = permutation_test(
            (col_a, col_b),
            statistic=_kendall_stat,
            permutation_type="pairings",
            n_resamples=n_bootstraps,
            random_state=rng,
        )
        rho_perm = permutation_test(
            (col_a, col_b),
            statistic=_spearman_stat,
            permutation_type="pairings",
            n_resamples=n_bootstraps,
            random_state=rng,
        )

        # 3. Standardized Bootstrap CIs (method='percentile' or 'BCa')
        tau_boot = bootstrap(
            (col_a, col_b),
            statistic=_kendall_stat,
            paired=True,
            confidence_level=confidence_level,
            n_resamples=n_bootstraps,
            method="percentile",
            random_state=rng,
        )
        rho_boot = bootstrap(
            (col_a, col_b),
            statistic=_spearman_stat,
            paired=True,
            confidence_level=confidence_level,
            n_resamples=n_bootstraps,
            method="percentile",
            random_state=rng,
        )

        # Fill matrices symmetrically
        for metric, val, pval, ci in [
            (
                "kendall",
                tau_val,
                tau_perm.pvalue,
                tau_boot.confidence_interval,
            ),
            ("spearman", rho_val, rho_perm.pvalue, rho_boot.confidence_interval),
        ]:
            prefix = "kendall_tau" if metric == "kendall" else "spearman_rho"
            p_prefix = f"{metric}_pvalues"
            low_prefix = f"{metric}_ci_lower"
            high_prefix = f"{metric}_ci_upper"

            for u, v in [(col_a_name, col_b_name), (col_b_name, col_a_name)]:
                results[prefix].loc[u, v] = float(val)
                results[p_prefix].loc[u, v] = float(pval)
                results[low_prefix].loc[u, v] = float(ci.low)
                results[high_prefix].loc[u, v] = float(ci.high)

    return results


def _estimate_bivariate_polychoric(
    x: np.ndarray,
    y: np.ndarray,
    scale_categories: list[int] = None,
    smooth_zeros: bool = True,
) -> float:
    """Estimates the polychoric correlation between two ordinal series using

    Olsson's Two-Stage Maximum Likelihood Estimation (2S-MLE).
    """
    # 1. Zero-variance guard
    if np.std(x) == 0 or np.std(y) == 0:
        return 0.0

    # 2. Build contingency table
    if scale_categories is not None:
        cat_x = pd.Categorical(x, categories=scale_categories)
        cat_y = pd.Categorical(y, categories=scale_categories)
    else:
        cat_x = pd.Categorical(x, categories=np.sort(np.unique(x)))
        cat_y = pd.Categorical(y, categories=np.sort(np.unique(y)))

    freq_table = pd.crosstab(cat_x, cat_y, dropna=False).to_numpy(dtype=float)

    # 3. Small-sample Laplace smoothing for empty cells
    if smooth_zeros and np.any(freq_table == 0):
        freq_table = freq_table.copy()
        freq_table += 0.5

    n_total = np.sum(freq_table)
    if n_total == 0:
        return 0.0

    p_table = freq_table / n_total

    # 4. Latent marginal thresholds via inverse standard normal CDF
    p_x_cum = np.clip(np.cumsum(np.sum(p_table, axis=1))[:-1], 1e-6, 1.0 - 1e-6)
    p_y_cum = np.clip(np.cumsum(np.sum(p_table, axis=0))[:-1], 1e-6, 1.0 - 1e-6)

    tau_x = np.concatenate(([-np.inf], norm.ppf(p_x_cum), [np.inf]))
    tau_y = np.concatenate(([-np.inf], norm.ppf(p_y_cum), [np.inf]))

    r_dim, c_dim = freq_table.shape

    # 5. Bivariate Normal Negative Log-Likelihood
    def _neg_log_lik(rho: float) -> float:
        rho = np.clip(rho, -0.9999, 0.9999)
        cov = np.array([[1.0, rho], [rho, 1.0]])
        mvn = multivariate_normal(
            mean=[0.0, 0.0], cov=cov, allow_singular=True
        )

        log_lik = 0.0
        for r in range(r_dim):
            for c in range(c_dim):
                if freq_table[r, c] > 0:
                    # Joint probability inside the bin rectangle
                    p_cell = (
                        mvn.cdf([tau_x[r + 1], tau_y[c + 1]])
                        - mvn.cdf([tau_x[r], tau_y[c + 1]])
                        - mvn.cdf([tau_x[r + 1], tau_y[c]])
                        + mvn.cdf([tau_x[r], tau_y[c]])
                    )
                    log_lik += freq_table[r, c] * np.log(max(p_cell, 1e-12))
        return -log_lik

    # 6. Brent bounded scalar optimization
    res = minimize_scalar(
        _neg_log_lik, bounds=(-0.9999, 0.9999), method="bounded"
    )
    return float(res.x)


def compute_polychoric_matrix(
    df: pd.DataFrame,
    scale_categories: list[int] = None,
    smooth_zeros: bool = True,
    ensure_positive_definite: bool = True,
) -> pd.DataFrame:
    """Computes a symmetric, positive-definite Polychoric Correlation Matrix

    using Olsson's 2S-MLE and statsmodels.
    """
    clean_df = df.apply(pd.to_numeric, errors="coerce").dropna().astype(float)
    cols = list(clean_df.columns)
    n_cols = len(cols)
    poly_matrix = pd.DataFrame(
        np.eye(n_cols), index=cols, columns=cols, dtype=float
    )

    # Compute pairwise combinations
    for (i, col_a), (j, col_b) in combinations(enumerate(cols), 2):
        x_vals = clean_df[col_a].to_numpy()
        y_vals = clean_df[col_b].to_numpy()

        rho = _estimate_bivariate_polychoric(
            x=x_vals,
            y=y_vals,
            scale_categories=scale_categories,
            smooth_zeros=smooth_zeros,
        )
        poly_matrix.iloc[i, j] = poly_matrix.iloc[j, i] = rho

    # Project to nearest positive definite matrix if eigenvalues <= 0
    if ensure_positive_definite:
        evals = np.linalg.eigvalsh(poly_matrix.values)
        if np.any(evals <= 1e-7):
            adjusted = corr_nearest(poly_matrix.values, threshold=1e-7)
            poly_matrix = pd.DataFrame(adjusted, index=cols, columns=cols)

    return poly_matrix


def compute_pca_scree_analysis_fa(
    corr_matrix: pd.DataFrame,
    kaiser_cutoff: float = 1.0,
) -> dict:
    """Computes eigenvalues and Kaiser retention using factor_analyzer."""
    fa = FactorAnalyzer(
        n_factors=corr_matrix.shape[0], method="principal", is_corr_matrix=True
    )
    fa.fit(corr_matrix)

    ev, _ = fa.get_eigenvalues()
    pct_var, cum_var, _ = fa.get_factor_variance()

    n_kaiser = int(np.sum(ev >= kaiser_cutoff))

    summary_df = pd.DataFrame(
        {
            "Eigenvalue": ev,
            "% of Variance": pct_var * 100.0,
            "Cumulative %": cum_var * 100.0,
        },
        index=[f"Component {i + 1}" for i in range(len(ev))],
    )

    return {
        "eigenvalues": ev,
        "n_factors_kaiser": n_kaiser,
        "variance_table": summary_df,
        "total_variance_kaiser": (
            summary_df["Cumulative %"].iloc[n_kaiser - 1]
            if n_kaiser > 0
            else 0.0
        ),
    }


def perform_horn_parallel_analysis(
    df: pd.DataFrame = None,
    corr_matrix: pd.DataFrame = None,
    n_observations: int = 27,
    n_simulations: int = 1000,
    percentile: float = 95.0,
    random_state: int = 42,
    plot_comparison: bool = True,
    save_path: str = None,
) -> dict:
    """Performs Horn's Parallel Analysis against Monte Carlo simulated uncorrelated noise."""
    rng = np.random.default_rng(random_state)

    if df is not None:
        clean_df = (
            df.apply(pd.to_numeric, errors="coerce").dropna().astype(float)
        )
        n_obs, n_vars = clean_df.shape
        emp_corr = (
            clean_df.corr().to_numpy()
            if corr_matrix is None
            else corr_matrix.to_numpy()
        )
    elif corr_matrix is not None:
        n_obs = n_observations
        n_vars = corr_matrix.shape[0]
        emp_corr = corr_matrix.to_numpy()
    else:
        raise ValueError("Must provide either `df` or `corr_matrix`.")

    # 1. Empirical Eigenvalues
    emp_eigenvalues = np.sort(np.linalg.eigvalsh(emp_corr))[::-1]

    # 2. Monte Carlo Simulation of Uncorrelated Gaussian Datasets (N x P)
    synthetic_evals = np.zeros((n_simulations, n_vars))
    for i in range(n_simulations):
        noise = rng.standard_normal(size=(n_obs, n_vars))
        noise_corr = np.corrcoef(noise, rowvar=False)
        synthetic_evals[i, :] = np.sort(np.linalg.eigvalsh(noise_corr))[::-1]

    simulated_mean = np.mean(synthetic_evals, axis=0)
    simulated_percentile = np.percentile(synthetic_evals, percentile, axis=0)

    # 3. Decision
    retained_mask = emp_eigenvalues > simulated_percentile
    n_factors_retained = int(np.sum(retained_mask))

    comparison_df = pd.DataFrame(
        {
            "Empirical Eigenvalue": emp_eigenvalues,
            f"Synthetic ({int(percentile)}th %ile)": simulated_percentile,
            "Synthetic (Mean)": simulated_mean,
            "Retain": ["Yes" if r else "No" for r in retained_mask],
        },
        index=[f"Component {i + 1}" for i in range(n_vars)],
    )

    if plot_comparison:
        plt.figure(figsize=(8, 5))
        comp_indices = np.arange(1, n_vars + 1)
        plt.plot(
            comp_indices,
            emp_eigenvalues,
            "o-",
            color="#1f77b4",
            linewidth=2,
            label="Empirical Eigenvalues",
        )
        plt.plot(
            comp_indices,
            simulated_percentile,
            "s--",
            color="#d62728",
            linewidth=1.75,
            label=f"Simulated {int(percentile)}th %ile",
        )
        plt.plot(
            comp_indices,
            simulated_mean,
            "^:",
            color="#7f7f7f",
            linewidth=1.25,
            label="Simulated Mean",
        )
        plt.axhline(
            y=1.0,
            color="black",
            linestyle=":",
            alpha=0.5,
            label=r"Kaiser ($\lambda = 1.0$)",
        )
        plt.title("Horn's Parallel Analysis vs. Scree Trace")
        plt.xlabel("Component Number")
        plt.ylabel("Eigenvalue")
        plt.xticks(comp_indices)
        plt.legend()
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300)
        plt.show()

    return {
        "n_factors_parallel": n_factors_retained,
        "empirical_eigenvalues": emp_eigenvalues,
        "synthetic_percentile_eigenvalues": simulated_percentile,
        "synthetic_mean_eigenvalues": simulated_mean,
        "comparison_table": comparison_df,
    }


def perform_abels_pca_oblimin(
    df: pd.DataFrame = None,
    corr_matrix: pd.DataFrame = None,
    n_factors: int = None,
    kaiser_cutoff: float = 1.0,
    rotation: str = "oblimin",
) -> dict:
    """Performs PCA with Direct Oblimin rotation using factor_analyzer.

    Accepts raw DataFrame or precomputed correlation matrix.
    """
    if corr_matrix is not None:
        input_data = corr_matrix.copy()
        is_corr = True
        cols = list(corr_matrix.columns)
    elif df is not None:
        input_data = df.copy()
        is_corr = False
        cols = list(df.columns)
    else:
        raise ValueError("Must provide either `df` or `corr_matrix`.")

    # 1. Ensure correlation matrix is positive semi-definite if passing precomputed matrix
    if is_corr:
        evals_raw = np.linalg.eigvalsh(input_data.values)
        if np.any(evals_raw <= 1e-7):
            input_data = pd.DataFrame(
                corr_nearest(input_data.values, threshold=1e-7),
                index=cols,
                columns=cols,
            )

    # 2. Preliminary fit to determine Kaiser factor count if not specified
    if n_factors is None:
        fa_initial = FactorAnalyzer(
            n_factors=len(cols),
            method="principal",
            rotation=None,
            is_corr_matrix=is_corr,
        )
        fa_initial.fit(input_data)
        evals_all, _ = fa_initial.get_eigenvalues()
        n_factors = int(np.sum(evals_all >= kaiser_cutoff))
        n_factors = max(1, min(n_factors, len(cols)))

    # 3. Fit target PCA with specified rotation
    fa = FactorAnalyzer(
        n_factors=n_factors,
        method="principal",
        rotation=rotation if n_factors > 1 else None,
        is_corr_matrix=is_corr,
    )
    fa.fit(input_data)

    # 4. Extract metrics
    evals, _ = fa.get_eigenvalues()
    pct_var, cum_var, _ = fa.get_factor_variance()
    comp_names = [f"Component {i + 1}" for i in range(n_factors)]

    variance_summary = pd.DataFrame(
        {
            "Eigenvalue": evals,
            "% of Variance": (evals / len(cols)) * 100.0,
            "Cumulative %": np.cumsum((evals / len(cols)) * 100.0),
        },
        index=[f"Component {i + 1}" for i in range(len(cols))],
    )

    # 5. Pattern Loadings, Factor Correlations (Phi), and Structure Matrix
    pattern_loadings = pd.DataFrame(
        fa.loadings_, index=cols, columns=comp_names
    )

    if fa.phi_ is not None and n_factors > 1:
        factor_correlations = pd.DataFrame(
            fa.phi_, index=comp_names, columns=comp_names
        )
        structure_matrix = pd.DataFrame(
            fa.structure_, index=cols, columns=comp_names
        )
    else:
        factor_correlations = pd.DataFrame(
            np.eye(n_factors), index=comp_names, columns=comp_names
        )
        structure_matrix = pattern_loadings.copy()

    # 6. Component assignments via max absolute loading
    primary_assignments = pattern_loadings.abs().idxmax(axis=1)
    clustering = {comp: [] for comp in comp_names}
    for var_name, comp_assigned in primary_assignments.items():
        clustering[comp_assigned].append(var_name)

    return {
        "n_factors_retained": n_factors,
        "eigenvalues": evals,
        "variance_summary": variance_summary,
        "pattern_loadings": pattern_loadings,
        "structure_matrix": structure_matrix,
        "factor_correlations": factor_correlations,
        "component_assignments": clustering,
    }


def _calc_alpha_from_corr(r: np.ndarray) -> float:
    """Computes standardized/ordinal alpha directly from a correlation matrix."""
    p = r.shape[0]
    if p <= 1:
        return np.nan
    sum_total = np.sum(r)
    trace_r = np.trace(r)
    return float((p / (p - 1.0)) * (1.0 - (trace_r / sum_total)))


def _calc_omega_from_data_or_corr(
    data: np.ndarray, is_corr: bool = False
) -> tuple[float, np.ndarray, np.ndarray]:
    """Extracts McDonald's Omega and 1-factor parameters using factor_analyzer."""
    p = data.shape[1] if not is_corr else data.shape[0]

    fa = FactorAnalyzer(
        n_factors=1,
        method="principal",
        rotation=None,
        is_corr_matrix=is_corr,
    )

    if is_corr:
        evals = np.linalg.eigvalsh(data)
        if np.any(evals <= 1e-6):
            data = corr_nearest(data, threshold=1e-6)

    fa.fit(data)
    loadings = np.clip(fa.loadings_.flatten(), -0.999, 0.999)

    # Align positive orientation
    if np.sum(loadings) < 0:
        loadings = -loadings

    uniquenesses = np.maximum(fa.get_uniquenesses(), 1e-4)

    sum_l = np.sum(loadings)
    denom = (sum_l**2) + np.sum(uniquenesses)
    omega = float((sum_l**2) / denom) if denom > 0 else 0.0

    return omega, loadings, uniquenesses


def compute_ordinal_reliability(
    df: pd.DataFrame,
    poly_corr_matrix: pd.DataFrame = None,
    n_bootstraps: int = 2000,
    alpha_ci: float = 0.05,
    random_state: int = 42,
) -> dict:
    """Computes Ordinal/Standard Cronbach's Alpha and McDonald's Omega

    with robust bootstrap CIs using pingouin, factor_analyzer, and scipy.
    """
    clean_df = (
        df.apply(pd.to_numeric, errors="coerce").dropna().astype(float).copy()
    )
    cols = list(clean_df.columns)
    data_arr = clean_df.to_numpy()
    conf_level = 1.0 - alpha_ci
    rng = np.random.default_rng(random_state)

    # 1. Standard Alpha & Parametric/Bootstrap CI (via Pingouin)
    std_alpha, std_alpha_ci = pg.cronbach_alpha(
        data=clean_df, ci=conf_level
    )

    # 2. Point Estimates for Ordinal Alpha and McDonald's Omega
    if poly_corr_matrix is not None:
        target_corr = poly_corr_matrix.to_numpy(dtype=float)
        ordinal_alpha = _calc_alpha_from_corr(target_corr)
        omega_val, lambda_loadings, unique_vars = _calc_omega_from_data_or_corr(
            target_corr, is_corr=True
        )
    else:
        target_corr = clean_df.corr().to_numpy()
        ordinal_alpha = _calc_alpha_from_corr(target_corr)
        omega_val, lambda_loadings, unique_vars = _calc_omega_from_data_or_corr(
            data_arr, is_corr=False
        )

    # 3. Standardized Bootstrap CIs for Ordinal Alpha & Omega (via SciPy)
    def _boot_ord_alpha_stat(sample_rows, axis=0):
        # sample_rows is received along axis from scipy bootstrap
        sample_corr = np.corrcoef(sample_rows, rowvar=False)
        return _calc_alpha_from_corr(sample_corr)

    def _boot_omega_stat(sample_rows, axis=0):
        try:
            om, _, _ = _calc_omega_from_data_or_corr(sample_rows, is_corr=False)
            return om
        except Exception:
            return np.nan

    # Vectorized bootstrap execution
    ord_alpha_boot = bootstrap(
        (data_arr,),
        statistic=_boot_ord_alpha_stat,
        confidence_level=conf_level,
        n_resamples=n_bootstraps,
        method="percentile",
        random_state=rng,
    )

    omega_boot = bootstrap(
        (data_arr,),
        statistic=_boot_omega_stat,
        confidence_level=conf_level,
        n_resamples=n_bootstraps,
        method="percentile",
        random_state=rng,
    )

    loadings_summary = pd.DataFrame(
        {
            "Standardized_Loading (λ)": lambda_loadings,
            "Error_Variance (θ)": unique_vars,
        },
        index=cols,
    )

    return {
        "ordinal_cronbach_alpha": ordinal_alpha,
        "ordinal_alpha_ci_95": (
            float(ord_alpha_boot.confidence_interval.low),
            float(ord_alpha_boot.confidence_interval.high),
        ),
        "composite_reliability_omega": omega_val,
        "omega_ci_95": (
            float(omega_boot.confidence_interval.low),
            float(omega_boot.confidence_interval.high),
        ),
        "standard_cronbach_alpha": float(std_alpha),
        "standard_alpha_ci_95": (float(std_alpha_ci[0]), float(std_alpha_ci[1])),
        "single_factor_loadings": loadings_summary,
    }


def compare_metric_groups(
    df: pd.DataFrame,
    grouping_series: pd.Series,
    group_name: str = "Group",
    padjust_method: str = "bonf",
) -> dict:
    """Performs non-parametric group comparisons and post-hoc pairwise tests

    using pingouin and pandas.
    """
    metrics = list(df.columns)
    combined_df = df.copy()
    combined_df[group_name] = grouping_series.values
    combined_df = combined_df.dropna(subset=[group_name])

    unique_groups = combined_df[group_name].unique()
    k_groups = len(unique_groups)

    if k_groups < 2:
        raise ValueError(
            f"Need at least 2 distinct groups to compare, found: {unique_groups}"
        )

    # 1. Descriptive Statistics via Pandas Aggregation
    def _iqr(x):
        return x.quantile(0.75) - x.quantile(0.25)

    desc_df = (
        combined_df.groupby(group_name)[metrics]
        .agg(
            N="count",
            Median="median",
            IQR=_iqr,
            Mean="mean",
            SD=lambda x: x.std(ddof=1),
        )
        .stack(level=1)
        .reset_index()
        .rename(columns={"level_1": "Metric"})
        .sort_values(by=["Metric", group_name])
    )

    # 2. Statistical Tests (Mann-Whitney U or Kruskal-Wallis via Pingouin)
    test_rows = []
    pairwise_frames = []

    for metric in metrics:
        if k_groups == 2:
            g1_name, g2_name = unique_groups[0], unique_groups[1]
            g1_vals = combined_df[combined_df[group_name] == g1_name][
                metric
            ].dropna()
            g2_vals = combined_df[combined_df[group_name] == g2_name][
                metric
            ].dropna()

            mwu_res = pg.mwu(g1_vals, g2_vals, alternative="two-sided")
            pval = mwu_res["p-val"].iloc[0]
            # Pingouin gives Common Language Effect Size (cles) and rank-biserial (RBC)
            rbc = mwu_res["RBC"].iloc[0]

            test_rows.append(
                {
                    "Metric": metric,
                    "Test": "Mann-Whitney U",
                    "Statistic (U)": mwu_res["U-val"].iloc[0],
                    "p_value": pval,
                    "Effect_Size (RBC)": rbc,
                    "Significant (p < .05)": pval < 0.05,
                }
            )

        else:
            # Kruskal-Wallis test
            kw_res = pg.kruskal(
                data=combined_df, dv=metric, between=group_name
            )
            pval = kw_res["p-unc"].iloc[0]
            h_stat = kw_res["H"].iloc[0]
            n_total = len(combined_df[metric].dropna())

            # Epsilon-squared effect size
            eps_sq = (
                max(0.0, (h_stat - k_groups + 1) / (n_total - k_groups))
                if n_total > k_groups
                else np.nan
            )

            test_rows.append(
                {
                    "Metric": metric,
                    "Test": "Kruskal-Wallis",
                    "Statistic (H)": h_stat,
                    "p_value": pval,
                    "Effect_Size (eta_sq)": eps_sq,
                    "Significant (p < .05)": pval < 0.05,
                }
            )

            # Post-hoc Pairwise Tests with Correction
            pw_res = pg.pairwise_tests(
                data=combined_df,
                dv=metric,
                between=group_name,
                parametric=False,
                padjust=padjust_method,
            )
            pw_res["Metric"] = metric
            pairwise_frames.append(pw_res)

    test_df = pd.DataFrame(test_rows)
    pairwise_df = pd.concat(pairwise_frames, ignore_index=True) if pairwise_frames else None

    return {
        "descriptive_stats": desc_df,
        "overall_tests": test_df,
        "pairwise_comparisons": pairwise_df,
    }


def compute_metric_distribution_summary(
    df: pd.DataFrame,
    scale_values: list = [0, 1, 2],
    include_parametric: bool = True,
) -> pd.DataFrame:
    """Computes frequency distributions, proportions, medians, IQRs, and

    parametric descriptives using vectorized pandas operations.
    """
    # 1. Ensure discrete categorical scale across all columns
    cat_df = df.apply(lambda col: pd.Categorical(col, categories=scale_values))

    # 2. Vectorized Counts and Percentages
    counts_df = (
        cat_df.apply(lambda col: col.value_counts(dropna=True))
        .T.add_prefix("Count_")
    )
    n_valid = df.notna().sum()
    pcts_df = (
        counts_df.div(n_valid, axis=0) * 100.0
    ).round(1).rename(columns=lambda c: c.replace("Count_", "Pct_"))

    # 3. Interleaved Frequency Columns: Count_0, Pct_0, Count_1, Pct_1, ...
    freq_cols = []
    for val in scale_values:
        freq_cols.extend([f"Count_{val}", f"Pct_{val}"])
    freq_summary = pd.concat([counts_df, pcts_df], axis=1)[freq_cols]

    # 4. Ordinal / Non-parametric Aggregations
    q25 = df.quantile(0.25)
    q75 = df.quantile(0.75)
    ord_stats = pd.DataFrame(
        {
            "N": n_valid,
            "Median": df.median(),
            "Q1 (25%)": q25,
            "Q3 (75%)": q75,
            "IQR": q75 - q25,
        },
        index=df.columns,
    )

    # 5. Assemble Summary Table
    summary_df = pd.concat([ord_stats[["N"]], freq_summary, ord_stats.drop(columns=["N"])], axis=1)

    # 6. Optional Parametric Descriptives
    if include_parametric:
        summary_df["Mean"] = df.mean().round(2)
        summary_df["SD"] = df.std(ddof=1).round(2)

    summary_df.index.name = "Metric"
    return summary_df


def _compute_gap_statistic(
    data: np.ndarray,
    k_range: list[int],
    n_bootstraps: int = 500,
    random_state: int = 42,
) -> tuple[dict, dict, int]:
    """Computes Tibshirani's Gap Statistic with 1-SE rule using scikit-learn."""
    rng = np.random.default_rng(random_state)
    n_samples, n_features = data.shape

    # Bounding box reference limits
    feat_min, feat_max = data.min(axis=0), data.max(axis=0)

    gap_stats = {}
    gap_std_errs = {}

    def _get_inertia(x, k, seed):
        km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(x)
        # Sum of pairwise distances in clusters equals 2 * N * inertia
        return km.inertia_

    for k in k_range:
        orig_inertia = _get_inertia(data, k, random_state)
        log_orig = np.log(max(orig_inertia, 1e-12))

        # Bootstrap uniform null references
        null_log_inertias = []
        for b in range(n_bootstraps):
            null_data = rng.uniform(feat_min, feat_max, size=(n_samples, n_features))
            null_inertia = _get_inertia(null_data, k, seed=None)
            null_log_inertias.append(np.log(max(null_inertia, 1e-12)))

        mean_null_log = np.mean(null_log_inertias)
        sd_k = np.std(null_log_inertias, ddof=1)
        gap_stats[k] = float(mean_null_log - log_orig)
        gap_std_errs[k] = float(sd_k * np.sqrt(1.0 + 1.0 / n_bootstraps))

    # 1-SE optimal rule: smallest k where Gap(k) >= Gap(k+1) - SE(k+1)
    k_list = sorted(gap_stats.keys())
    optimal_k = k_list[0]
    for i in range(len(k_list) - 1):
        curr_k, next_k = k_list[i], k_list[i + 1]
        if gap_stats[curr_k] >= (gap_stats[next_k] - gap_std_errs[next_k]):
            optimal_k = curr_k
            break
    else:
        optimal_k = max(gap_stats, key=gap_stats.get)

    return gap_stats, gap_std_errs, optimal_k


def compute_pre_clustering_triangulation(
    df_metrics: pd.DataFrame = None,
    poly_corr_matrix: pd.DataFrame = None,
    k_range: range = range(2, 6),
    n_bootstraps_gap: int = 500,
    random_state: int = 42,
    plot_results: bool = True,
    save_path: str = None,
) -> dict:
    """Performs pre-clustering triangulation across variables/metrics using

    Silhouette analysis, Gap Statistic, and Ward's dendrogram.
    """
    if poly_corr_matrix is not None:
        corr_df = poly_corr_matrix.copy()
    elif df_metrics is not None:
        corr_df = df_metrics.corr()
    else:
        raise ValueError("Must provide either `df_metrics` or `poly_corr_matrix`.")

    metric_names = list(corr_df.columns)
    p_vars = len(metric_names)
    valid_k_list = [k for k in k_range if k < p_vars]

    # 1. Angular correlation distance: D = sqrt(2 * (1 - R))
    r_clipped = np.clip(corr_df.values, -1.0, 1.0)
    dist_matrix = np.sqrt(np.maximum(2.0 * (1.0 - r_clipped), 0.0))
    np.fill_diagonal(dist_matrix, 0.0)
    condensed_dist = squareform(dist_matrix, checks=False)

    # 2. Silhouette Analysis (scikit-learn)
    silhouette_scores = {}
    for k in valid_k_list:
        labels = AgglomerativeClustering(
            n_clusters=k, metric="precomputed", linkage="average"
        ).fit_predict(dist_matrix)
        silhouette_scores[k] = float(
            silhouette_score(dist_matrix, labels, metric="precomputed")
        )

    optimal_k_silhouette = max(silhouette_scores, key=silhouette_scores.get)

    # 3. Gap Statistic (scikit-learn KMeans engine)
    gap_stats, gap_std_errs, optimal_k_gap = _compute_gap_statistic(
        data=dist_matrix,
        k_range=valid_k_list,
        n_bootstraps=n_bootstraps_gap,
        random_state=random_state,
    )

    # 4. Ward's Hierarchical Linkage (scipy)
    ward_linkage_matrix = linkage(condensed_dist, method="ward")

    # 5. Triangulation Summary Table
    triangulation_table = pd.DataFrame(
        {
            "Silhouette Score": [
                silhouette_scores.get(k, np.nan) for k in valid_k_list
            ],
            "Gap Statistic": [gap_stats.get(k, np.nan) for k in valid_k_list],
            "Gap Std Error": [
                gap_std_errs.get(k, np.nan) for k in valid_k_list
            ],
        },
        index=[f"k = {k}" for k in valid_k_list],
    )

    # 6. Visual Dashboard
    if plot_results:
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # Silhouette Curve
        axes[0].plot(
            valid_k_list,
            [silhouette_scores[k] for k in valid_k_list],
            "o-",
            color="#1f77b4",
            linewidth=2,
            markersize=7,
        )
        axes[0].axvline(
            optimal_k_silhouette,
            color="firebrick",
            linestyle="--",
            label=f"Max Silhouette (k={optimal_k_silhouette})",
        )
        axes[0].set_title("Silhouette Analysis", fontsize=12, pad=10)
        axes[0].set_xlabel("Number of Clusters (k)")
        axes[0].set_ylabel("Average Silhouette Width")
        axes[0].set_xticks(valid_k_list)
        axes[0].grid(axis="both", linestyle=":", alpha=0.6)
        axes[0].legend()

        # Gap Statistic Curve
        axes[1].errorbar(
            valid_k_list,
            [gap_stats[k] for k in valid_k_list],
            yerr=[gap_std_errs[k] for k in valid_k_list],
            fmt="s-",
            color="#2ca02c",
            ecolor="#2ca02c",
            capsize=4,
            linewidth=2,
            markersize=6,
        )
        axes[1].axvline(
            optimal_k_gap,
            color="firebrick",
            linestyle="--",
            label=f"Optimal Gap (k={optimal_k_gap})",
        )
        axes[1].set_title(f"Gap Statistic ($B={n_bootstraps_gap}$)", fontsize=12, pad=10)
        axes[1].set_xlabel("Number of Clusters (k)")
        axes[1].set_ylabel("Gap Statistic $Gap(k)$")
        axes[1].set_xticks(valid_k_list)
        axes[1].grid(axis="both", linestyle=":", alpha=0.6)
        axes[1].legend()

        # Ward's Dendrogram
        dendrogram(
            ward_linkage_matrix,
            labels=metric_names,
            orientation="right",
            ax=axes[2],
            leaf_font_size=9,
        )
        axes[2].set_title("Ward Hierarchical Dendrogram", fontsize=12, pad=10)
        axes[2].set_xlabel("Ward Dissimilarity Distance")

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300)
        plt.show()

    return {
        "triangulation_table": triangulation_table,
        "optimal_k_silhouette": optimal_k_silhouette,
        "optimal_k_gap": optimal_k_gap,
        "ward_linkage": ward_linkage_matrix,
        "distance_matrix": pd.DataFrame(
            dist_matrix, index=metric_names, columns=metric_names
        ),
    }


def project_and_cluster_publications(
    df_metrics: pd.DataFrame,
    pca_results: dict,
    n_clusters: int,
    clustering_method: str = "kmeans",
    publication_labels: list[str] = None,
    random_state: int = 42,
    plot_scatter: bool = True,
    save_path: str = None,
) -> dict:
    """Standardizes metrics, projects publications onto PCA axes via least-squares,

    and clusters them using scikit-learn.
    """
    n_obs, _ = df_metrics.shape
    pub_ids = (
        list(publication_labels)
        if publication_labels is not None
        else [f"Pub_{i+1:02d}" for i in range(n_obs)]
    )

    # 1. Preprocessing Pipeline: Impute -> Standardize
    preprocessor = make_pipeline(
        SimpleImputer(strategy="most_frequent"), StandardScaler()
    )
    z_scores = preprocessor.fit_transform(
        df_metrics.apply(pd.to_numeric, errors="coerce")
    )

    # 2. Project Component Scores via Least-Squares (Z = Scores @ Loadings.T)
    pattern_loadings = np.nan_to_num(
        pca_results["pattern_loadings"].to_numpy(dtype=float), nan=0.0
    )
    n_components = pattern_loadings.shape[1]
    comp_names = [f"PCA_Comp_{i+1}" for i in range(n_components)]

    # Solves argmin || z_scores - scores @ pattern_loadings.T ||
    projected_scores, _, _, _ = np.linalg.lstsq(
        pattern_loadings, z_scores.T, rcond=None
    )
    projected_scores = projected_scores.T  # Shape: (N, n_components)

    df_scores = pd.DataFrame(
        projected_scores, index=pub_ids, columns=comp_names
    )

    # 3. Clustering Execution
    if clustering_method.lower() == "kmeans":
        cluster_model = KMeans(
            n_clusters=n_clusters, random_state=random_state, n_init=20
        )
    elif clustering_method.lower() == "hierarchical":
        cluster_model = AgglomerativeClustering(
            n_clusters=n_clusters, linkage="ward"
        )
    else:
        raise ValueError(
            "clustering_method must be either 'kmeans' or 'hierarchical'"
        )

    cluster_labels = cluster_model.fit_predict(projected_scores)
    cluster_tags = [f"Cluster_{c+1}" for c in cluster_labels]

    # 4. Aggregated Summaries
    df_clustered = df_scores.copy()
    df_clustered.insert(0, "Assigned_Cluster", cluster_tags)
    df_clustered.insert(0, "Publication", pub_ids)

    cluster_centroids = (
        df_scores.groupby(cluster_tags).mean().sort_index()
    )

    df_original = df_metrics.copy()
    df_original.index = pub_ids
    df_original["Assigned_Cluster"] = cluster_tags
    metric_profiles = df_original.groupby("Assigned_Cluster").mean(
        numeric_only=True
    ).sort_index()

    # 5. Scatter Plot Dashboard
    if plot_scatter and n_components >= 2:
        plt.figure(figsize=(9, 6))
        palette = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
        ]

        for c_idx in range(n_clusters):
            c_tag = f"Cluster_{c_idx+1}"
            mask = np.array(cluster_tags) == c_tag
            plt.scatter(
                projected_scores[mask, 0],
                projected_scores[mask, 1],
                s=90,
                color=palette[c_idx % len(palette)],
                label=c_tag,
                alpha=0.85,
                edgecolors="black",
                linewidth=0.7,
            )

        for i, label in enumerate(pub_ids):
            plt.annotate(
                label,
                (projected_scores[i, 0], projected_scores[i, 1]),
                fontsize=8,
                xytext=(5, 4),
                textcoords="offset points",
                alpha=0.8,
            )

        plt.title(
            f"Publications Projected onto PCA Dimensions (k={n_clusters} Clusters)",
            fontsize=13,
            pad=12,
        )
        plt.xlabel(f"{comp_names[0]} (Factor Score)", fontsize=11)
        plt.ylabel(f"{comp_names[1]} (Factor Score)", fontsize=11)
        plt.axhline(0, color="gray", linestyle=":", linewidth=1, alpha=0.7)
        plt.axvline(0, color="gray", linestyle=":", linewidth=1, alpha=0.7)
        plt.grid(axis="both", linestyle=":", alpha=0.5)
        plt.legend(frameon=True, loc="best")
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300)
        plt.show()

    return {
        "projected_scores": df_scores,
        "clustered_publications": df_clustered,
        "cluster_centroids": cluster_centroids,
        "cluster_metric_profiles": metric_profiles,
    }


if __name__ == "__main__":

    # --- Execution Parameters ---
    excel_path = "/home/david/Git/SDLT_Number_Crunching/Data/adiya_isla_scores.xlsx"
    consensus_sheet_name = "Whole"
    n_bootstraps = 1000
    random_state = 42

    pipeline_results = {}

    # =========================================================================
    # STEP 1: Load and Clean Dataset
    # =========================================================================
    xls = pd.ExcelFile(excel_path)
    df_consensus = pd.read_excel(xls, sheet_name=consensus_sheet_name)

    # Extract publication identifiers for downstream clustering labels
    pub_labels = df_consensus["Publication"].astype(str).tolist() if "Publication" in df_consensus.columns else None

    # Filter and clean numerical rating columns (drop metadata columns like 'Publication')
    if "Publication" in df_consensus.columns:
        df_metrics = df_consensus.drop(columns=["Publication"])
    else:
        df_metrics = df_consensus.copy()

    # Ensure purely numeric float type
    df_metrics = df_metrics.apply(pd.to_numeric, errors="coerce").astype(float)
    pipeline_results["consensus_data"] = df_metrics


    # =========================================================================
    # STEP 2: Ordinal Correlation Matrices (Polychoric, Kendall, Spearman)
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 2: COMPUTING ORDINAL & POLYCHORIC CORRELATIONS")
    print("=" * 70)

    # Polychoric Matrix (projected to nearest positive definite)
    poly_corr = compute_polychoric_matrix(
        df=df_metrics,
        ensure_positive_definite=True,
    )
    pipeline_results["polychoric_matrix"] = poly_corr

    # Non-parametric Kendall's tau-b and Spearman's rho with permutation p-values and bootstrap CIs
    rank_corrs = compute_ordinal_correlations(
        df=df_metrics,
        n_bootstraps=n_bootstraps,
        random_state=random_state,
    )
    pipeline_results["ordinal_correlations"] = rank_corrs

    print("Polychoric Correlation Matrix (Nearest Positive Definite):")
    display(poly_corr.round(3))

    # Abels et al. screening: check if any metric fails to correlate (all |r| < 0.30)
    off_diag = poly_corr.values.copy()
    np.fill_diagonal(off_diag, 0.0)
    max_corr_per_item = pd.Series(np.max(np.abs(off_diag), axis=1), index=df_metrics.columns)
    low_corr_items = max_corr_per_item[max_corr_per_item < 0.30].index.tolist()

    if low_corr_items:
        print(f"\n[Warning] Items with max polychoric correlation < 0.30: {low_corr_items}")
    else:
        print("\nAll items meet the minimum correlation threshold (max |r| >= 0.30).")


    # =========================================================================
    # STEP 3: Matrix Factorability (KMO & Bartlett's Test with CIs)
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 3: MATRIX FACTORABILITY TESTS (KMO & BARTLETT'S SPHERICITY)")
    print("=" * 70)

    factorability = evaluate_kmo_bartlett_factorability(
        df=df_metrics,
        n_bootstraps=n_bootstraps,
        ci_level=0.95,
        random_state=random_state,
    )
    pipeline_results["factorability"] = factorability

    print(f"Overall KMO: {factorability['kmo_total']:.3f} | 95% CI: [{factorability['kmo_total_ci'][0]:.3f}, {factorability['kmo_total_ci'][1]:.3f}]")
    print(f"Bartlett's Chi-Square: {factorability['bartlett_chi_sq']:.2f} (p = {factorability['bartlett_p_val']:.4e})")
    print("\nItem-Level KMO with 95% Bootstrapped CIs:")
    display(factorability["kmo_per_item_cis"].round(3))


    # =========================================================================
    # STEP 4: Retention Diagnostics (Parallel Analysis & Scree Plot)
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 4: FACTOR RETENTION (HORN'S PARALLEL ANALYSIS & KAISER)")
    print("=" * 70)

    parallel_results = perform_horn_parallel_analysis(
        df=df_metrics,
        corr_matrix=poly_corr,
        n_simulations=1000,
        percentile=95.0,
        random_state=random_state,
        plot_comparison=True,
    )
    pipeline_results["parallel_analysis"] = parallel_results

    print("Horn's Parallel Analysis Eigenvalue Comparison:")
    display(parallel_results["comparison_table"].round(3))
    print(f"Optimal Factors Suggested by Parallel Analysis: {parallel_results['n_factors_parallel']}")


    # =========================================================================
    # STEP 5: PCA with Direct Oblimin Rotation (Abels et al. Model)
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 5: PCA WITH DIRECT OBLIMIN ROTATION")
    print("=" * 70)

    # Fits target PCA using factor_analyzer (uses parallel analysis factor count if desired, or defaults to Kaiser)
    retained_factors = parallel_results["n_factors_parallel"] if parallel_results["n_factors_parallel"] > 0 else None

    pca_results = perform_abels_pca_oblimin(
        corr_matrix=poly_corr,
        n_factors=retained_factors,
        kaiser_cutoff=1.0,
        rotation="oblimin",
    )
    pipeline_results["pca"] = pca_results

    print("Variance Explained Breakdown:")
    display(pca_results["variance_summary"].round(2))

    print("\nRotated Pattern Matrix (Table 3 Style):")
    display(pca_results["pattern_loadings"].round(3))

    print("\nStructure Matrix (Table 4 Style):")
    display(pca_results["structure_matrix"].round(3))

    print("\nComponent Correlations (Phi):")
    display(pca_results["factor_correlations"].round(3))


    # =========================================================================
    # STEP 6: Scale Internal Consistency & Reliability
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 6: ORDINAL SCALE RELIABILITY & COMPOSITE RELIABILITY")
    print("=" * 70)

    reliability = compute_ordinal_reliability(
        df=df_metrics,
        poly_corr_matrix=poly_corr,
        n_bootstraps=n_bootstraps,
        alpha_ci=0.05,
        random_state=random_state,
    )
    pipeline_results["reliability"] = reliability

    print(f"Ordinal Cronbach's Alpha: {reliability['ordinal_cronbach_alpha']:.3f} | 95% CI: [{reliability['ordinal_alpha_ci_95'][0]:.3f}, {reliability['ordinal_alpha_ci_95'][1]:.3f}]")
    print(f"McDonald's Omega (Composite): {reliability['composite_reliability_omega']:.3f} | 95% CI: [{reliability['omega_ci_95'][0]:.3f}, {reliability['omega_ci_95'][1]:.3f}]")
    print(f"Standard Cronbach's Alpha: {reliability['standard_cronbach_alpha']:.3f} | 95% CI: [{reliability['standard_alpha_ci_95'][0]:.3f}, {reliability['standard_alpha_ci_95'][1]:.3f}]")
    print("\nSingle-Factor CFA Standardized Loadings (Unidimensional Model):")
    display(reliability["single_factor_loadings"].round(3))


    # =========================================================================
    # STEP 7: Compute Pre-Clustering Triangulation
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 7: PRE-CLUSTERING TRIANGULATION (SILHOUETTE, GAP STATISTIC, WARD)")
    print("=" * 70)

    clustering_triangulation = compute_pre_clustering_triangulation(
        poly_corr_matrix=poly_corr,
        k_range=range(2, 6),
        n_bootstraps_gap=500,
        random_state=random_state,
        plot_results=True,
    )
    pipeline_results["clustering_triangulation"] = clustering_triangulation

    print("Triangulation Diagnostics Summary:")
    display(clustering_triangulation["triangulation_table"].round(3))


    # =========================================================================
    # STEP 8: Project & Cluster Publications
    # =========================================================================
    print("\n" + "=" * 70)
    print("STEP 8: PUBLICATION PROJECTION & CLUSTERING")
    print("=" * 70)

    optimal_k = clustering_triangulation["optimal_k_silhouette"]

    pca_clusters = project_and_cluster_publications(
        df_metrics=df_metrics,
        pca_results=pipeline_results["pca"],
        n_clusters=optimal_k,
        clustering_method="kmeans",
        publication_labels=pub_labels,
        random_state=random_state,
        plot_scatter=True,
    )
    pipeline_results["publication_clusters"] = pca_clusters

    print(f"Cluster Assignments (k={optimal_k}):")
    display(pca_clusters["clustered_publications"].head(10))

    print(f"\nMean Metric Profiles (0-2 Scale) Across Clusters:")
    display(pca_clusters["cluster_metric_profiles"].round(2))