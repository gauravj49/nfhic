#!/bin/bash

# USAGE: bash scripts/convert_hicPro_to_bedpe_format.sh /media/rad/HDD1/hic/anja/tcelldev/results/hic_results/matrix/raw/TcellDev-HSC-HSPCtoTcell-Rep3_mm_hic_pe_70000.matrix /media/rad/HDD1/hic/anja/tcelldev/results/hic_results/matrix/raw/TcellDev-HSC-HSPCtoTcell-Rep3_mm_hic_pe_70000_abs.bed /media/rad/HDD1/hic/anja/tcelldev/results/hic_results/matrix/raw/TcellDev-HSC-HSPCtoTcell-Rep3_mm_hic_pe_70000.bedpe

# Input argument and relevant directories
inMFile=${1:-"/media/rad/HDD1/hic/anja/tcelldev/results/hic_results/matrix/raw/TcellDev-HSC-HSPCtoTcell-Rep3_mm_hic_pe_70000.matrix"}
inBFile=${2:-"/media/rad/HDD1/hic/anja/tcelldev/results/hic_results/matrix/raw/TcellDev-HSC-HSPCtoTcell-Rep3_mm_hic_pe_70000_abs.bed"}
outFile=${3:-"/media/rad/HDD1/hic/anja/tcelldev/results/hic_results/matrix/raw/TcellDev-HSC-HSPCtoTcell-Rep3_mm_hic_pe_70000.bedpe"}

# Get the Hic data for each cellline
# 4) Convert raw matrices (hicpro format) to bedpe format
# 4.1) from another script
python scripts/conv_hicpro_mat.py ${inMFile} ${inBFile} > ${outFile}

# 4.2) Put a dot (.) in the 7th column
python -c "import sys; import pandas as pd; import datatable as dt; input_file=sys.argv[1]; outtxt_file=sys.argv[2]; bedpeDT = dt.fread(input_file, sep='\t', header=False, nthreads=16); bedpeDF = bedpeDT.to_pandas(); bedpeDF.insert(6, 'name', '.'); bedpeDF.to_csv(outtxt_file, header=False, index=False, sep='\t');" ${outFile} ${outFile}.tmp && mv ${outFile}.tmp ${outFile}


