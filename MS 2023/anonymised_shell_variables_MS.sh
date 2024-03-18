#!/bin/bash

### Global settings
## CPU threads - some software, e.g. samtools, count the number of CPUs used as n+1. Therefore, if n=6, they will actually use 7.
CPU=7

## set flag for max heat size for java executables
JRAM=-Xmx8000M #8GB

## Base directories
bd=~/Documents/MS/            # MS basedir
bda=${bd}analysis/  # analysis
bdr=${bd}raw/       # raw MS files  
bdpdb=${bd}pdb/     # protein databases
bdmzid=${bda}mzid/

### Software paths
# msgf+
msgfplus=~/software/MSGFPlus_v20230112/MSGFPlus.jar
# ThermoRawFileParser.exe
thermoRFP=~/software/ThermoRawFileParser1.4.2/ThermoRawFileParser.exe

### Input and output directories
## raw ThermoFisher Orbitrap Exploris 480 file directories
raw_plk89=${bdr}PLK089__D/
raw_plk89F=${bdr}PLK089__D_FAIMS/
raw_plk99=${bdr}PLK099__D/

## raw ThermoFisher Orbitrap Exploris 480 file directories
bda_mM=${bda}mzML/ # define basedir for mzMLs
mzml_plk89=${bda_mM}plk89/
mzml_plk89F=${bda_mM}plk89_FAIMS/
mzml_plk99=${bda_mM}plk99/

## protein databases
# dirs
pdb_top5_dir=${bdpdb}top5/
pdb_hp_dir=${bdpdb}hp/
pdb_hpsp_dir=${bdpdb}hp_sp/
# files
pdb_top5=${pdb_top5_dir}top5-crap-mycoplasma.fasta
pdb_hp=${pdb_hp_dir}uniprot_combined-crap-mycoplasma.fasta
pdb_hpsp=${pdb_hpsp_dir}uniprot_swissprot-crap-mycoplasma.fasta

## mzid - just proteins coded by ChIP-seq-identified genes
bda_m5=${bda}mzid/top5/ # define basedir for top5 mzids
mzid_plk89_top5=${bda_m5}plk89/
mzid_plk89F_top5=${bda_m5}plk89_FAIMS/
mzid_plk99_top5=${bda_m5}plk99/

## mzid - whole proteome (Uniprot Swissprot + TrEMBL)
bda_hp=${bda}mzid/hp/ # define basedir for hp mzids
mzid_plk89_hp=${bda_hp}plk89/
mzid_plk89F_hp=${bda_hp}plk89_FAIMS/
mzid_plk99_hp=${bda_hp}plk99/

## mzid - whole proteome (Uniprot Swissprot only)
bda_hpsp=${bda}mzid/hpsp/ # define basedir for hp mzids
mzid_plk89_hpsp=${bda_hpsp}plk89/
mzid_plk89F_hpsp=${bda_hpsp}plk89_FAIMS/
mzid_plk99_hpsp=${bda_hpsp}plk99/