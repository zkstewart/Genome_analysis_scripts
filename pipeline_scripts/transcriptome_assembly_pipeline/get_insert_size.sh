#!/bin/bash -l
#PBS -N insert_size
#PBS -l walltime=01:00:00
#PBS -l mem=10G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

# SCRIPT SETUP START
# >> Setup: Module loads
module load atg/picard/2.2.2

# >> Setup: Input file locations
SAMDIR=/home/stewarz2/plant_group/plant_rnaseq/analysis/avo/star_map/mac80
SAMFILE=Aligned.out.sam

# >> Setup: Prefixes
SPECIES=avo

# >> Setup: Amount of reads to sample
SAMPLESIZE=100000

# >> Setup: Generate helper script automatically
PICARDPYTHONHELPER="imetrics_rnaseq_densepeak.py"
echo "
#! python3
from pathlib import Path
imetrics_file = Path(r\"${SPECIES}.subset${SAMPLESIZE}.imetrics\")
nums = []
cols = []
skip = True
with open(imetrics_file, 'r') as file_in:
        for line in file_in:
                # Skip to relevant section of file
                if line.startswith('insert_size'):
                        skip = False
                        continue
                if skip == True:
                        continue
                if '\t' not in line:
                        continue
                # Hold onto data
                col = line.rstrip('\r\n ').split('\t')
                cols.append([int(col[0]), int(col[1])])
                nums.append(int(col[1]))
maximumest_index = [0, 0]
for i in range(len(cols)):
        if i < 5 or i + 6 > len(nums): # Skip edges
                continue
        summed_density = 0
        for x in range(i - 5, i + 6):
                summed_density += cols[x][1]
        if summed_density > maximumest_index[1]:
                maximumest_index = [cols[i][0], summed_density]
print(maximumest_index[0])
" > ${PICARDPYTHONHELPER}
sed -i '1d' ${PICARDPYTHONHELPER}
# SCRIPT SETUP END

# RUN START
# >> STEP 1: Subset SAM file
head -n ${SAMPLESIZE} ${SAMDIR}/${SAMFILE} > ${SPECIES}.subset${SAMPLESIZE}.sam

# >> STEP 2: Run picard to derive statistics
picard CollectInsertSizeMetrics H=${SPECIES}.subset${SAMPLESIZE}.histo I=${SPECIES}.subset${SAMPLESIZE}.sam O=${SPECIES}.subset${SAMPLESIZE}.imetrics

# >> STEP 3: Extract insert size details from imetrics output file
INSERTSIZE=$(python ${PICARDPYTHONHELPER})

# >> STEP 4: Report insert size
echo "INSERT_SIZE: ${INSERTSIZE}" > ${SPECIES}.insert_size.txt
# RUN END
