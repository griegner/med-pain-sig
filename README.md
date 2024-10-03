## med-pain-sig
> meditation effects on fMRI-based pain signatures

**publication**
```
@article{riegner2024mindfulness,
  title={Mindfulness meditation and placebo modulate distinct multivariate neural signatures to reduce pain},
  author={Riegner, Gabriel and Dean, Jon and Wager, Tor D and Zeidan, Fadel},
  journal={Biological Psychiatry},
  year={2024},
  publisher={Elsevier}
}
```

**project organization**
```
.
├── README.md                   <- this README file
├── config
│   ├── bold-dem.csv            <- study 1 demographics
│   ├── bold-filtered.csv       <- study 1 bold runs to process
│   ├── cbf-dem.csv             <- study 2 demographics
│   ├── cbf-filtered.csv        <- study 2 cbf runs to process
│   └── smk-config.yml          <- configuration file for workflow/snakemake
├── data
│   ├── bold-ps-agg.csv         <- study 1 pain signature response data
│   ├── cbf-ps-agg.csv          <- study 2 pain signature response data
│   ├── bold/                   <- study 1 bold glm estimates [-]
│   ├── cbf/                    <- study 2 cerebral blood flow estimates [-]
│   └── pain-sigs/              <- CANlab pain-predictive signatures [-]
├── figures/                    <- figures derived from notebooks/
├── notebooks/
│   ├── fig01.ipynb             <- plot neuroimaging data
│   ├── fig02to05.ipynb         <- plot and analyze pain signature responses
│   ├── fig02to05.py            <- helper functions for fig02to05.ipynb
│   └── sup01.ipynb             <- plot and analyze supplementary data
├── results/                    <- result tables derived from notebooks/
└── workflow
    ├── envs                    <- required python packages to reproduce full analysis pipeline
    ├── scripts
    │   ├── 00_bold-to-pe.py    <- calculate glm parameter estimates for study 1
    │   ├── 00_cbf-to-ps.py     <- apply pain signatures to study 2
    │   ├── 01_bold-to-ps.py    <- apply pain signatures to study 1
    │   └── apply_ps.py
    └── snakefile               <- documentation of snakemake workflow

[-] not included in this repository
```

The `data/*csv` files contain the long-form pain signature and behavioral data, from which the `figures/` and `results/` can be reproduced with these [notebooks/](notebooks/). Running the '*ipynb' files requires both Python and R packages, so a pre-installed distribution of the [conda package manager](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#) (Anaconda, Miniconda, or Miniforge) is needed.  

The pain signatures can be accessed on the [CANLAB website](https://sites.google.com/dartmouth.edu/canlab-brainpatterns/multivariate-brain-signatures).

The steps to reproduce the full analysis pipeline (from [ASLPrep](https://aslprep.readthedocs.io/en/latest/) and [fMRIPrep](https://fmriprep.org/en/stable/) outputs) are outlined in `workflow/snakefile`, but access to the ~500GB of neuroimaging data is not included here.

**installation**

clone this repository
```
git clone https://github.com/griegner/med-pain-sig.git  
cd med-pain-sig
```

install required dependencies
```
conda env create --file requirements.yml  
conda activate med-pain-sig
```
