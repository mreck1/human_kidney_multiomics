# Human kidney multiomics
This repository includes R/python3 scripts to reproduce all main and supplementary figures shown in the manuscript 'Multiomic analysis of human kidney disease identifies a tractable inflammatory, pro-fibrotic tubular cell phenotype'.

The main datasets (including raw and processed data in Seurat format) are available on GEO: GSE254185, GSE253439, GSE252584, GSE254187.

Intermediate data needed for some parts of the analysis is available on Zenodo:

https://zenodo.org/records/10849854?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjQxMjEyNmNhLTlkZjItNDNjMy04YzRhLWViZGFiYzA1NWFiYiIsImRhdGEiOnt9LCJyYW5kb20iOiI0YzA3OGNkM2NhMTdmZjExZGJkYzYyOGIwZDcxMzdhYyJ9.yRj9yQNZCMsoNp_PhC3IwPVVsKCMWbuMyYYq10fRJ06_UHzgLCBNUAb6tGmrzG8LPKJN2WioCwdC1arqKN_PxQ

utils.R includes helper function required to execute the code and needs to be loaded before. R package versions are available in the session info document. To run the RNA velocity analysis a standard installation of scVelo is required.
