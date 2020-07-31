# pwd
cd /media/rad/HDD1/hic/tcelldev

# Parameters for the script
species="mouse"
user="anja"
projName="tcelldev"
outdir="/media/rad/HDD1/hic"
jobdir="/home/rad/users/gaurav/projects/workflows/nfhic"

# Intersect CIS with TADs (70kb)
sort -u /media/rad/HDD1/hic/anja/tcelldev/results/analysis/cis/cisfinal.bed | sort -k1,1V -k2,2g -k3,3g > tmp && mv tmp /media/rad/HDD1/hic/anja/tcelldev/results/analysis/cis/cisfinal.bed


# Add chr to tads file
for f in $(ls /media/rad/HDD1/hic/anja/tcelldev/results/analysis/input/tads/70kb/*.bed)
do 
  echo ${f}
  # sed -i 's/^/chr/' ${f}
  sort -k1,1V -k2,2g -k3,3g ${f} -o ${f}
done

# head /media/rad/HDD1/hic/anja/tcelldev/results/analysis/input/tads/70kb/TcellDev-HSC-HSPCtoTcell-Rep3_mm_hic_pe_70000_iced_TADs_boundaries.bed
# chr1    4235000 4305000 B15318  -0.501479907564 .
# chr1    5145000 5215000 B15331  -0.653814654936 .
# chr1    7245000 7315000 B15361  -0.501315085020 .

# head /media/rad/HDD1/hic/anja/tcelldev/results/analysis/cis/cisfinal.bed
# chr1    13306488        13311900
# chr1    15803386        15809782
# chr1    21958793        21968141

intersectBed -a /media/rad/HDD1/hic/anja/tcelldev/results/analysis/input/tads/70kb/TcellDev-HSC-HSPCtoTcell-Rep3_mm_hic_pe_70000_iced_TADs_boundaries.bed -b /media/rad/HDD1/hic/anja/tcelldev/results/analysis/cis/cisfinal.bed -wo > HSCtads70kb_cisfinal.txt

# head HSCtads70kb_cisfinal.txt
# chr1    60935000        61005000        B16128  -0.706124819622 .       chr1    60989614        60998962        9348
# chr1    63035000        63105000        B16158  -0.781011604133 .       chr1    63049617        63052569        2952
# chr1    69755000        69825000        B16254  -0.900527237053 .       chr1    69762952        69775744        12792
# chr1    80745000        80815000        B16411  -0.522267024797 .       chr1    80754223        80761111        6888

# Remove columns for input file. This is temporary only. I`ll introduce header which will avoid this problem
cut -f1-3,7-9 HSCtads70kb_cisfinal.txt > HSCtads70kb_cisfinal_filtered.txt
# head HSCtads70kb_cisfinal_filtered.txt
# chr1    60935000        61005000        chr1    60989614        60998962
# chr1    63035000        63105000        chr1    63049617        63052569
# chr1    69755000        69825000        chr1    69762952        69775744

# Get the TAD-CIS distance plot
name       <- args[1];
out        <- args[2];
mainTitle  <- args[3];
span       <- args[4]; # encode or genomewide
retBinInfo <- args[5]; # 0 = draw plots, 1 = return bins and counts 

Rscript /home/rad/users/gaurav/projects/workflows/nfhic/scripts/R_drawDistancePlots.R /media/rad/HDD1/hic/anja/tcelldev/results/analysis/HSCtads70kb_cisfinal_filtered.txt /media/rad/HDD1/hic/anja/tcelldev/results/analysis/HSCtads70kb_cisfinal_filtered_plots.png HSCtads70kb_cisfinal genomewide 0
