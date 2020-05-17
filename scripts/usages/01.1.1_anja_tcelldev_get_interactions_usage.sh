# pwd
cd /media/rad/HDD1/hic/tcelldev
nextflow run /home/rad/users/gaurav/projects/workflows/nfhic --reads 'fastq/*_R{1,2}.fastq.gz' --genome GRCm38 --restriction_site C^TAG --bin_size 10000,20000,40000,70000,100000,200000,500000,1000000 --saveReference -name tcelldevHiC --max_memory '2.Gb'

