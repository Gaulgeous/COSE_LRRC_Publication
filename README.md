# Empirical Legislative Analysis of Digital Forensics Literature

An open-source statistical analysis pipeline for evaluating academic publications scored using an ordinal Likert rubric. This repository accompanies the (pending) publication:

> **"Reconciling Performance and Legal Admissibility: Evaluating Forensic Evidence Acquisition Through an Empirical Legislative Lens"**

The codebase analyses how digital forensics literature adheres to legislative and regulatory requirements for forensic evidence admissibility in a court of law.

---

## Overview

The workflow is split into two sequential components:

1. **Score Consolidation (Python — `src/`):**
   - Aggregates multi-annotator, dual-round inter-rater reliability (IRR) scores from raw spreadsheets (`Data/final_scores.xlsx`) into a consolidated workbook (`Data/consolidated_results.xlsx`).
   - Prioritises Round 2 consensus ratings with fallback to Round 1, resolves ties deterministically, and filters incomplete entries.

2. **Statistical Analysis & Visualisation (R — `r_scripts/`):**
   - **Category Profile Means & Centroid Distances:** Computes mean scorecard profiles across publication categories and evaluates pairwise Euclidean distances to identify sub-domains with similar legislative coverage.
   - **Inter-Rater Reliability (IRR):** Computes ordinal Krippendorff's alpha across all criteria sheets, bootstrapped over 1,000 resamples for 95% confidence intervals and reproducibility.
   - **Principal Component Analysis (PCA):** Evaluates factorability (Bartlett's test, KMO), dimension retention (Kaiser-Guttman & Horn's parallel analysis), oblique direct oblimin rotation, and scale internal consistency (Cronbach's alpha, ordinal alpha, McDonald's omega).
   - **Publication Mapping & Clustering:** Projects all publications onto the extracted principal component space, validates cluster counts (Elbow, Silhouette, Gap statistic, and `NbClust` consensus), profiles cluster centroids, and generates publication figures.

---

## Repository Structure

```text
.
├── Data/
│   ├── final_scores.xlsx              # Raw multi-round annotator Likert scores
│   ├── consolidated_results.xlsx      # Consolidated workbook produced by Python pipeline
│   └── old_consolidated_results.xlsx  # Archived prior consolidation output
├── Images/                            # Exported PDF figures (PCA plots, donuts, eigenvalues)
├── r_scripts/
│   ├── main.R                         # Primary R analysis pipeline entrypoint
│   ├── helpers.R                      # Label normalisation, PSD guards, Excel ingestion, preliminary stats
│   ├── statistical_means.R            # Category profile means and Euclidean distance heatmaps
│   ├── krippendorff_alpha.R           # Bootstrapped ordinal Krippendorff's alpha IRR pipeline
│   ├── simple_pca.R                   # Streamlined polychoric/Kendall PCA, parallel analysis, reliability
│   ├── pca_helpers.R                  # Exact/multi-start direct oblimin PCA utilities
│   ├── clustering_helpers.R           # PCA score projection, cluster triangulation, NbClust, profiling
│   ├── visualisations.R               # 2D/3D PCA plots, correlation heatmaps, donut & eigenvalue charts
│   ├── r_scripts.Rproj                # RStudio project file (anchors working directory and here())
│   ├── renv.lock                      # Locked R package versions (R 4.3.3, 76 CRAN packages)
│   ├── .Rprofile                      # Auto-activates renv when R starts in r_scripts/
│   └── renv/                          # renv bootstrap scripts and configuration
├── src/
│   └── consolidate_results.ipynb      # Notebook to consolidate multi-round annotator scores
├── pyproject.toml                     # Python project metadata and dependencies
├── uv.lock                            # Locked Python dependency versions
├── .python-version                    # Target Python version (3.13)
├── run_r_scripts.sh                   # Convenience script to restore renv and run main.R
├── LICENSE                            # MIT License
└── README.md                          # Project documentation
```

---

## Prerequisites & Environment Setup

Both language environments use lockfiles so dependencies are reproducible without manual package hunting.

### 1. R Environment (`renv`)
* **Required:** [R](https://cran.r-project.org/) (`4.3.3` recommended; system build tools and libraries such as `gcc`/`gfortran`, `libcurl`, `libssl`, and `libxml2` may be required on Linux to compile CRAN packages).
* **Restore dependencies:**
  ```bash
  cd r_scripts
  Rscript -e "renv::restore()"
  ```
  *(If opening `r_scripts/r_scripts.Rproj` in RStudio, `renv` activates automatically via `.Rprofile` and prompts to run `renv::restore()`.)*

### 2. Python Environment (`uv`)
* **Required:** Python `>= 3.13` and [`uv`](https://docs.astral.sh/uv/) (recommended), or `pip`.
* **Sync dependencies:**
  ```bash
  uv sync
  ```
  This creates `.venv/` with exact versions from `uv.lock` (`pandas`, `openpyxl`, `numpy`, `ipykernel`, etc.).

---

## Running the Code

### Step 1: Consolidate Raw IRR Scores (Python)
1. Ensure raw annotator data is present at `Data/final_scores.xlsx`.
2. Open `src/consolidate_results.ipynb` in VS Code or Jupyter (selecting the `.venv` Python kernel).
3. Run the notebook cells (specifically the multi-round consolidation cell) to generate `Data/consolidated_results.xlsx`.
   * *Note:* `consolidate_results.ipynb` expects to be executed with the repository root as the working directory (check `excel_path` in cell 2 if running from `src/`).

### Step 2: Run Statistical Analysis & Generate Plots (R)
Using the convenience script from the repository root:
```bash
bash run_r_scripts.sh
```
Or manually from `r_scripts/`:
```bash
cd r_scripts
Rscript main.R
```

### Viewing Output Plots
* **Batch execution (`Rscript main.R`):** Because `Rscript` runs non-interactively, all generated plots are written sequentially to **`r_scripts/Rplots.pdf`**.
* **Interactive execution (VS Code R extension / RStudio):** Opening `r_scripts/main.R` and running lines interactively renders plots directly in the VS Code **Plot Viewer** (with `REditorSupport.r` / `httpgd`) or the RStudio **Plots** pane.
* **Reproducibility:** `main.R` fixes `set.seed(42)` and runs `n_boot <- 1000` bootstrap/Monte Carlo iterations for Krippendorff's alpha and Horn's parallel analysis.

---

## License

This project is licensed under the **MIT License** — the standard permissive open-source license widely used for peer-reviewed academic and scientific software. See [LICENSE](LICENSE) for full terms.
