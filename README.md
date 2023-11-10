## med-pain-sig
> meditation effects on fMRI-based pain signatures

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
│   ├── bold/                   <- study 1 bold glm estimates *
│   ├── cbf/                    <- study 2 cerebral blood flow estimates *
│   └── pain-sigs/              <- CANlab pain-predictive signatures *
├── figures/                    <- figures derived from notebooks/
├── notebooks/
│   ├── fig01.ipynb         <- plot neuroimaging data
│   ├── fig02to05.ipynb     <- plot and analyze pain signature responses
│   ├── fig02to05.py        <- helper functions for fig02to05.ipynb
│   └── sup01.ipynb         <- plot and analyze supplimentary data
├── results/
│   ├── bold-ps-agg.csv         <- study 1 pain signature response data
│   └── cbf-ps-agg.csv          <- study 2 pain signature response data
└── workflow
    ├── envs                    <- required python packages to reproduce analyses
    │   └── main.yml            <- conda env create --name med_pain_sig --file workflow/envs/
    ├── scripts
    │   ├── 00_bold-to-pe.py    <- calculate glm parameter estimates for study 1
    │   ├── 00_cbf-to-ps.py     <- apply pain signatures to study 2
    │   ├── 01_bold-to-ps.py    <- apply pain signatures to study 1
    │   └── apply_ps.py
    └── snakefile               <- documentation of snakemake workflow

* not version controlled
```