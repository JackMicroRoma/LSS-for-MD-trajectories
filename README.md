# LSS Conformational Sampling Pipeline for Ion Channels

LSS-based conformational sampling pipeline for ion channel proteins from molecular dynamics trajectories. Processes all-atom MD trajectories through SNRV, MDN, and WGANGP models to generate enhanced structural ensembles. Includes monomer extraction, LSS training notebooks, and PCA/RMSD/RMSF analysis scripts comparing MD and predicted trajectories.

---

## Overview

This project applies the **Latent Space Sampling (LSS)** framework to a diverse set of ion channels and membrane proteins to generate enhanced structural ensembles that capture slow conformational dynamics beyond standard MD timescales.

The pipeline consists of three sequential neural network models:

- **SNRV** (Slow Relaxation Variable network) — learns collective variables (CVs) from backbone pairwise distances
- **MDN** (Mixture Density Network) — propagator that models CV dynamics
- **WGANGP** (Wasserstein GAN with gradient penalty) — decoder that reconstructs synthetic 3D trajectories conditioned on CVs

For multimeric channels, individual monomers are extracted as independent units retaining only mainchain atoms (N, Cα, C, O), enabling tractable training on consumer hardware.

---

## Proteins

| Protein | Type | Monomers |
|---------|------|----------|
| aBK | Homotetrameric K⁺ channel | 4 |
| cBK | Homotetrameric K⁺ channel | 4 |
| KcsA | Homotetrameric K⁺ channel | 4 |
| Kv4.2 | Homotetrameric K⁺ channel | 4 |
| MthK | Homotetrameric K⁺ channel | 4 |
| TREK-2 | Homodimeric K⁺ channel | 2 |
| NOMPC | Mechanosensitive channel | 4 |
| 1YGG, 3HRF, 3II0, 3VMM, 4RK9 | Monomeric proteins | 1 |
| AGO_mut, AGO_wt | Monomeric proteins | 1 |
| ca_zhujinwei_mut, ca_zhujinwei_wt | Monomeric proteins | 1 |

---

## Repository Structure

```
LSS-main/
├── data_2T/                        # Input data (not tracked by git)
│   └── {protein}/
│       ├── monomers/               # Extracted monomer GROs and XTCs
│       │   ├── MonomerA_mc.gro
│       │   └── MonomerA/
│       │       └── *.xtc
│       └── lss_outputs/            # LSS-predicted structures
│           └── {protein}_MonomerA/
│               └── {protein}_MonomerA_run{0-9}.pdb
│
├── pca_2T/                         # Analysis outputs (not tracked by git)
│   ├── scripts/                    # Master copy of all analysis scripts
│   └── {protein}/                  # Per-protein working directory
│       ├── pairs.tsv               # Pairing table: GRO / XTC / pred PDB
│       ├── run_all.sh              # Full pipeline launcher
│       ├── pca_data/               # GROMACS covariance + eigenvector data
│       ├── rmsd_data/              # RMSD time series
│       ├── rmsf_analysis/          # Per-residue RMSF data
│       ├── eigenvec_results/       # Eigenvector comparison matrices
│       └── results/                # Final plots and summary tables
│
├── ADP_backbone_LSS.ipynb          # LSS training notebook
├── WLALL_backbone_LSS.ipynb        # LSS sampling notebook
└── README.md
```

---

## Dependencies

- Python ≥ 3.9
- [GROMACS](https://www.gromacs.org/) (`gmx_mpi`)
- [MDAnalysis](https://www.mdanalysis.org/)
- NumPy, SciPy, Matplotlib
- PyTorch (for LSS training)

Install Python dependencies:

```bash
pip install mdanalysis numpy scipy matplotlib torch
```

---

## Workflow

### 1. Monomer extraction

Extract individual monomers from each multimeric channel trajectory, retaining mainchain atoms only:

```bash
python find_monomers.py
```

### 2. LSS training and sampling

Run the Jupyter notebooks for each monomer:

- `ADP_backbone_LSS.ipynb` — trains SNRV + MDN + WGANGP
- `WLALL_backbone_LSS.ipynb` — samples synthetic trajectories

LSS outputs (predicted structures) are saved to `data_2T/{protein}/lss_outputs/`.

### 3. Analysis pipeline setup

All analysis scripts live in `pca_2T/scripts/`. Setup creates a self-contained working directory for each protein with a `pairs.tsv` pairing table that maps each MD trajectory to its corresponding LSS prediction:

```bash
cd pca_2T/scripts/

# Multimeric channels (4 monomers × 10 runs = 40 pairs)
./setup_protein.sh aBK

# Monomeric proteins (N runs = N pairs)
./setup_protein_mono.sh 1YGG

# Or set up all proteins at once
./run_all_proteins.sh
```

### 4. Run analysis

```bash
cd pca_2T/aBK/
./run_all.sh
```

This executes the full pipeline automatically:

| Step | Script | Output |
|------|--------|--------|
| Covariance overlap | `run_covar_overlap.sh` | Overlap score per pair |
| Eigenvector comparison | `run_compare_eigenvec.sh` | 10×10 cosine similarity matrices |
| PCA plots | `run_pca_analysis.sh` | Variance + distribution plots |
| KL divergence (PCA) | `py_eigenvec_stats_kl.py` | KL divergence of PC projections |
| RMSD | `run_rmsd.sh` | RMSD time series per pair |
| RMSD analysis | `py_rmsd_analysis.py` | Free energy + KL divergence |
| PCA free energy | `py_pca_fe_analysis.py` | Free energy landscapes |
| Eigenvector ratio | `py_eigenvec_ratio.py` | Diagonal alignment ratio |
| RMSF | `run_rmsf.sh` | Per-residue flexibility |
| RMSF analysis | `py_rmsf.py` | RMSF comparison + KL divergence |

---

## Analysis Description

Each pair compares one MD trajectory against one LSS-predicted trajectory for the same monomer and run index. The following metrics are computed:

**Covariance overlap** — measures how well the predicted ensemble reproduces the dominant collective motions of the MD trajectory (range 0–1, higher is better).

**Eigenvector cosine similarity** — 10×10 matrix of dot products between the top 10 MD and predicted eigenvectors. Diagonal dominance indicates the prediction captures the same modes in the same order.

**Eigenvector diagonal-alignment ratio** — for each column of the similarity matrix, the ratio min(h,v)/max(h,v) where h is the column index and v is the row of the maximum value. A ratio of 1.0 means perfect diagonal alignment.

**KL divergence** — computed for PC1/PC2 projections and RMSD distributions, quantifying how different the MD and predicted conformational distributions are (lower is better).

**RMSF** — per-residue root mean square fluctuation compared between MD and prediction, reporting both the profile correlation and KL divergence of the distributions.

---

## Notes

- `pairs.tsv` is auto-generated by the setup scripts and maps each `pair_idx` to its GRO topology, XTC trajectory, and predicted PDB file. Shell scripts and Python scripts read from this table — no hardcoded filenames anywhere.
- For proteins where the LSS prediction PDB includes a terminal OXT atom not present in the GROMACS mainchain selection, a 1–2 atom mismatch is handled automatically by truncating both eigenvectors to the shorter length before computing dot products.
- The `#covar.log.N#` files that accumulate in working directories are Emacs backup files and can be safely deleted: `rm '#covar.log.'*'#'`

