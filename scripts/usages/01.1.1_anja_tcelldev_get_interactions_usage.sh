# pwd
cd /media/rad/HDD1/hic/tcelldev

# Parameters for the script
species="mouse"
user="anja"
projName="tcelldev"
outdir="/media/rad/HDD1/hic"
jobdir="/home/rad/users/gaurav/projects/workflows/nfhic"

#########################################################################################
# 1) RUN THE PIPELINE
#########################################################################################
# Parameters to run the pipeline
projDir="${outdir}/${user}/${projName}"
fastqDir="${projDir}/fastq"; mkdir -p ${fastqDir}

nextflow run /home/rad/users/gaurav/projects/workflows/nfhic --reads "${fastqDir}/*_R{1,2}.fastq.gz" --genome GRCm38 --restriction_site C^TAG --bin_size 10000,20000,40000,70000,100000,200000,500000,1000000 --saveReference -name "${projName}HiC" --max_memory '2.GB' --skipIce

#########################################################################################
# 2) PREPROCESS OUTPUT OF THE PIPELINE FOR DOWNSTREAM ANALYSIS
#########################################################################################
# 2.1) Run ICE on raw matrices
rawMatrixDir=${projDir}/results/hic_results/matrix/raw  ; mkdir -p ${rawMatrixDir}
iceMatrixDir=${projDir}/results/hic_results/matrix/iced ; mkdir -p ${iceMatrixDir}
for inputRawMatrix in ${rawMatrixDir}/*.matrix;
do 
  bname=$(basename ${inputRawMatrix} .matrix)
  echo "Running ICE for ${bname} ..."
  ice --filter_low_counts_perc 0.02  --filter_high_counts_perc 0 --max_iter 100 --eps 0.1 --remove-all-zeros-loci --output-bias 1 --verbose 1 --results_filename ${iceMatrixDir}/${bname}_iced.matrix ${inputRawMatrix}
done

# 2.2) Convert ICED matrices to h5 format
echo "Converting to h5..."
totalSamples=$(ls ${rawMatrixDir}/*.matrix | wc -l)
i=0
for inputRawMatrix in ${rawMatrixDir}/*.matrix;
do
  bname=$(basename ${inputRawMatrix} .matrix)
  h5outname=${iceMatrixDir}/${bname}_iced.h5
  iceMatFile=${iceMatrixDir}/${bname}_iced.matrix
  hicproBed=${rawMatrixDir}/${bname}_abs.bed
  i=$((i+1)); echo "${i} of ${totalSamples}: ${bname}"
  hicConvertFormat  -m ${iceMatFile} --bedFileHicpro ${hicproBed} --inputFormat hicpro --outputFormat h5 -o ${h5outname}
  echo ""
done

# 2.3) Get the TADs
# --minDepth
# Minimum window length (in bp) to be considered to the left and to the right of each Hi-C bin. This number should be at least 3 times as large as the bin size of the Hi-C matrix.
# --maxDepth
# Maximum window length to be considered to the left and to the right of the cut point in bp. This number should around 6-10 times as large as the bin size of the Hi-C matrix.
# --step
# Step size when moving from –minDepth to –maxDepth. Note, the step size grows exponentially as minDeph + (step * int(x)**1.5) for x in [0, 1, …] until it reaches maxDepth. For example, selecting step=10,000, minDepth=20,000 and maxDepth=150,000 will compute TAD-scores for window sizes: 20,000, 30,000, 40,000, 70,000 and 100,000

echo "Get the TADs..."
totalSamples=$(ls ${iceMatrixDir}/*.h5 | wc -l)
i=0
tadsDir=${projDir}/results/hic_results/tads ; mkdir -p ${tadsDir}
tadsLogDir=${tadsDir}/logs ; mkdir -p ${tadsLogDir}
tadLogsFile=${tadsLogDir}/${projName}_TAD_calling_hicFindTADs.log; rm ${tadLogsFile}
echo "-------------------------------------------------" >> ${tadLogsFile} 2>&1
for h5matrix in $(ls ${iceMatrixDir}/*.h5);
do
  # Get the variables and parameters
  bname=$(basename ${h5matrix} _iced.h5)             # TcellDev-HSC-HSPCtoTcell-Rep3_mm_hic_pe_70000
  binsize=$(echo ${bname}|rev|cut -d '_' -f 1 | rev) # 70000 (binsize from the filename)
  mindepth=$(bc <<< 3*${binsize})
  maxdepth=$(bc <<< 8*${binsize})
  outfilename=${tadsDir}/${bname}_iced_TADs
  # Find TADs
  i=$((i+1)); echo "${i} of ${totalSamples}: ${bname}"
  hicFindTADs -m ${h5matrix} --outPrefix ${outfilename} --minDepth  ${mindepth} --maxDepth ${maxdepth} --step ${binsize} --thresholdComparisons 0.05 --delta 0.01 --correctForMultipleTesting fdr -p 8
  # Add to log file
  echo "- ${i} of ${totalSamples}: ${bname}" >> ${tadLogsFile} 2>&1
  echo "hicFindTADs -m ${h5matrix} --outPrefix ${outfilename} --minDepth  ${mindepth} --maxDepth ${maxdepth} --step ${binsize} --thresholdComparisons 0.05 --delta 0.01 --correctForMultipleTesting fdr -p 8" >> ${tadLogsFile} 2>&1
  echo "-------------------------------------------------" >> ${tadLogsFile} 2>&1
  echo "" >> ${tadLogsFile} 2>&1
done

# 2.4) Get AB compartments
echo "Get AB compartments..."
totalSamples=$(ls ${iceMatrixDir}/*.h5 | wc -l)
i=0
pcaDir=${projDir}/results/hic_results/compartments/pca ; mkdir -p ${pcaDir}
compartmentsDir=${projDir}/results/hic_results/compartments ; mkdir -p ${compartmentsDir}
compartmentsLogDir=${compartmentsDir}/logs ; mkdir -p ${compartmentsLogDir}
compartmentLogsFile=${compartmentsLogDir}/${projName}_Compartments.log; rm ${compartmentLogsFile}
echo "#-------------------------------------------------" >> ${compartmentLogsFile} 2>&1
for h5matrix in $(ls ${iceMatrixDir}/*.h5);
do
  # Get the variables and parameters
  bname=$(basename ${h5matrix} _iced.h5)             # TcellDev-HSC-HSPCtoTcell-Rep3_mm_hic_pe_70000
  pca1outfile=${pcaDir}/${bname}_iced_pca1.bedgraph
  pca2outfile=${pcaDir}/${bname}_iced_pca2.bedgraph
  outfilename=${compartmentsDir}/${bname}_iced_compartments_global_signal.png
  outmatrix=${compartmentsDir}/${bname}_iced_compartments.npz
   # Computes PCA eigenvectors for a Hi-C matrix
  i=$((i+1)); echo "${i} of ${totalSamples}: ${bname}"
  hicPCA --matrix ${h5matrix} --outputFileName ${pca1outfile} ${pca2outfile}
  # Compute the global compartmentalization signal
  hicCompartmentalization --obsexp_matrices ${h5matrix} --pca ${pca1outfile} -o ${outfilename} --outputMatrix ${outmatrix}
  # Add to log file
  echo "#- ${i} of ${totalSamples}: ${bname}" >> ${compartmentLogsFile} 2>&1
  echo "#-------------------------------------------------" >> ${compartmentLogsFile} 2>&1
  echo "hicPCA –-matrix ${h5matrix} --outputFileName  ${pca1outfile} ${pca2outfile}" >> ${compartmentLogsFile} 2>&1
  echo "" >> ${compartmentLogsFile} 2>&1
  echo "hicCompartmentalization --obsexp_matrices ${h5matrix} --pca ${pca1outfile} -o ${outfilename} --outputMatrix ${outmatrix}" >> ${compartmentLogsFile} 2>&1
  echo "" >> ${compartmentLogsFile} 2>&1
  echo "#-------------------------------------------------" >> ${compartmentLogsFile} 2>&1
done
