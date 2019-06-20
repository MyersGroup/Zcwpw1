computeMatrix scale-regions \
    --metagene \
    -S deeptools/bigwigs/WTCHG_538916_221156_vs_WTCHG_538916_217108_CPM_log2.bw \
      deeptools/bigwigs/WTCHG_538916_223180_vs_WTCHG_538916_220144_CPM_log2.bw \
      deeptools/bigwigs/SRA_Altemose2015_SRR5627138_AND_SRA_Altemose2015_SRR5627139_vs_SRA_Altemose2015_SRR5627140_CPM_log2.bw \
      deeptools/bigwigs/SRA_Altemose2015_SRR5627152_vs_SRA_Altemose2015_SRR5627143_CPM_log2.bw \
      deeptools/bigwigs/SRA_Altemose2015_SRR5627149_vs_SRA_Altemose2015_SRR5627143_CPM_log2.bw \
      -R Homo_sapiens.GRCh38.94.gtf \
      -m 8000 \
      -a 2000 \
      -b 2000 \
      -bs 100 \
      -p 15 \
      --skipZeros \
      --missingDataAsZero \
      -o gtf_exon_test.gz


plotProfile \
    -m gtf_exon_test.gz \
    --numPlotsPerRow 2 \
    --samplesLabel ChpHA_ZHA-vs-In_ZHA_CPM ChpHA_ZHA_hP9V5-vs-Input hPrdm9_Nterm-vs-Input H3K4wPrdm9_v_in H3K36wPrdm9_v_in \
    --yMax 0.15 0.2 0.1 0.3 0.05 \
    -out gtf_exon_test.pdf
