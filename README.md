# Circadian benchmarking for CircaN R package

Code for the benchmarking of the CircaN package previous to submission.
The complete file structure is:

``` bash
├── data
│   ├── geo_submission
│   └── raw
├── docs
│   └── manuscript
├── references
│   ├── Hughes
│   ├── RNAseq
│   ├── Salva
│   └── Subsampling
├── reports
├── results
│   ├── Hughes
│   ├── NLS_benchmarking
│   ├── pval_qval
│   ├── qPCR
│   ├── RNAseq
│   ├── Salva
│   └── Subsampling
└── src
    ├── CircaN
    └── dependencies
``` 

\* Most of these directories contain only data files, and so are not available in the git repo. A level 2 tree is provided to help understand the project structure.       
\*\* Please note that I reestructured all the directories to make them fit better with the Cookiecutter Data Science standard. The paths to files etc may have changed. Following is a list of changes I made:

- SRC/ is now src/
- SW/ is now src/dependencies/
- FIGURES/ is now reports/figures/
- GEO/ is now data/geo_submission/
- DOC/ is now references/
- RESULTS is now results/
- RAW is now data/raw/
- PAPER is now docs/manuscript/

I maintained the results directory as is instead of moving the corresponding files into reports and data for convenience.