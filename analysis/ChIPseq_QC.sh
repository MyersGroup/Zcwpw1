# Perform QC and filtering on BAMs,

export PATH=/homes/wells/data-on-saxony/single-cell/sequencing/software/fastqc/FastQC/:$PATH
export PATH=/homes/wells/data-on-saxony/single-cell/sequencing/software/bedtools/bedtools2/bin/:$PATH
export PATH=/homes/wells/data-on-saxony/single-cell/sequencing/software/samtools/samtools-1.7/:$PATH
export PATH=/homes/wells/data-on-saxony/single-cell/sequencing/software/:$PATH

# install mosdepth
# wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
# tar jxf htslib-1.9.tar.bz2
#
# ./configure --prefix=/homes/wells/data-on-saxony/single-cell/sequencing/software/htslib_1.9/
# make
# make install
# # finish install mosdepth

export LD_LIBRARY_PATH=/homes/wells/data-on-saxony/single-cell/sequencing/software/htslib_1.9/lib/:$LD_LIBRARY_PATH


###### FastQC
for sample in ${samples[@]}
do
(
  fastqc $path/WTCHG_${expid}_${sample}_1.fastq.gz -o qc
  fastqc $path/WTCHG_${expid}_${sample}_2.fastq.gz -o qc
) &
done



# zcat $path/WTCHG_${expid}_${sample}_1.fastq.gz | echo $((`wc -l`/4))
# same as fastqc number

###### Statistics on Lib Complexity from Unmapped reads

for sample in ${samples[@]}
do
(
  java -jar ../software/picard.jar EstimateLibraryComplexity INPUT=$path/WTCHG_${expid}_$sample.bam OUTPUT=qc/WTCHG_${expid}_${sample}_LibComplex.txt
) &
done

###### Statistics on mapping

for sample in ${samples[@]}
do
(
  samtools flagstat $path/WTCHG_${expid}_$sample.bam > qc/WTCHG_${expid}_${sample}_FlagStat.txt
) &
done


###### remove duplicates

for sample in ${samples[@]}
do
(
  java -jar ../software/picard.jar MarkDuplicates \
  INPUT=$path/WTCHG_${expid}_$sample.bam \
  OUTPUT=undup/WTCHG_${expid}_$sample.bam \
  METRICS_FILE=undup/WTCHG_${expid}_${sample}_metrics.txt \
  REMOVE_DUPLICATES=true
) &
done


###### filter low quality

for sample in ${samples[@]}
do
(
  samtools view undup/WTCHG_${expid}_${sample}.bam -q 20 -b -o HQmap/WTCHG_${expid}_${sample}.bam
  samtools index HQmap/WTCHG_${expid}_${sample}.bam
)&
done


###### use bed index to find genome wide coverage

for sample in ${samples[@]}
do
goleft indexcov --sex chrX,chrY -d plots/${sample}_indexcov HQmap/WTCHG_${expid}_${sample}.bam
done


###### Create blacklists of high coverage reigons

for sample in ${samples[@]}
do
Rscript gen_blacklist.R plots/${sample}_indexcov/ ${sample} plots/
done


###### Combine blacklists

for sample in ${samples[@]}
do
cat ${sample}_blacklist.bed >> ${expid}_combined_blacklist.bed
done

sort -k1,1V -k2,2n -u ${expid}_combined_blacklist.bed -o ${expid}_combined_blacklist.bed


###### Filter high cov reigons

for sample in ${samples[@]}
do
(bedtools intersect -v -a HQmap/WTCHG_${expid}_${sample}.bam -b ${expid}_combined_blacklist.bed > higcovfiltered/${sample}.bam
  samtools index higcovfiltered/${sample}.bam) &
  done

###### how many reads are left at the end

for sample in ${samples[@]}
do
(
  echo $sample
  samtools view -c -F 260 higcovfiltered/${sample}.bam
  #samtools idxstats higcovfiltered/${sample}.bam | awk '{print SUM += $3} END { print SUM}'
  #samtools flagstat higcovfiltered/${sample}.bam
)
done

for sample in ${samples[@]}
do
(
samtools view higcovfiltered/${sample}.bam -L beds/chr19_sample.bed -b -o sample_regions/sample_${expid}_${sample}_chr19.bam
samtools index sample_regions/sample_${expid}_${sample}_chr19.bam
) &
done



###### Compare depth/coverage distribution for raw & deduped

for sample in ${samples[@]}
do
(
MOSDEPTH_PRECISION=5 mosdepth -n -t 4 mosdepth/${expid}_${sample} higcovfiltered/${sample}.bam
MOSDEPTH_PRECISION=5 mosdepth -n -t 4 mosdepth/${expid}_raw_${sample} $path/WTCHG_${expid}_$sample.bam
) &
done

# need to add to list of input sample id's here manually :(
Rscript plot_mosdepth.R ${expid} "278163|276139|217108|220144"
