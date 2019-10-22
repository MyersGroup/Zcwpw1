Rscript pipelines/MultiProfilePlot.R bwplots/AlleleSpecificity.pdf 'Zcwpw1 binds Prdm9 sites in an allele specific manner' 8 8 \
bwprofilesNorm/WTCHG_538916_223180_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
bwprofilesNorm/WTCHG_538916_223180_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627145_AND_NA15-SRR5627144_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
bwprofilesNorm/WTCHG_538916_224192_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
bwprofilesNorm/WTCHG_538916_224192_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627145_AND_NA15-SRR5627144_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
'1) Zcwpw1 cotransfected with human Prdm9, at human Prdm9 sites' \
'2) Zcwpw1 cotransfected with human Prdm9, at chimp Prdm9 sites' \
'3) Zcwpw1 cotransfected with chimp Prdm9, at human Prdm9 sites' \
'4) Zcwpw1 cotransfected with chimp Prdm9, at chimp Prdm9 sites'

# Use Motif Centered and Stranded for Human allele
Rscript pipelines/MultiProfilePlot.R bwplots/AlleleSpecificityMCT.pdf 'Zcwpw1 binds Prdm9 sites in an allele specific manner (Motif Centered and Stranded MCS)' 10 7 \
bwprofilesNorm/WTCHG_538916_223180_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL_MotifCenteredStranded.profile \
bwprofilesNorm/WTCHG_538916_223180_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627145_AND_NA15-SRR5627144_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
bwprofilesNorm/WTCHG_538916_224192_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL_MotifCenteredStranded.profile \
bwprofilesNorm/WTCHG_538916_224192_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627145_AND_NA15-SRR5627144_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
'1) Zcwpw1 cotransfected with human Prdm9, at human Prdm9 MCS sites' \
'2) Zcwpw1 cotransfected with human Prdm9, at chimp Prdm9 sites' \
'3) Zcwpw1 cotransfected with chimp Prdm9, at human Prdm9 MCS sites' \
'4) Zcwpw1 cotransfected with chimp Prdm9, at chimp Prdm9 sites'


# vs input
Rscript pipelines/MultiProfilePlot.R bwplots/AlleleSpecificity_vInput.pdf 'Zcwpw1 binds Prdm9 sites in an allele specific manner' 8 8 \
bwprofilesNorm/WTCHG_538916_223180_VS_WTCHG_538916_220144_AT_SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
bwprofilesNorm/WTCHG_538916_223180_VS_WTCHG_538916_220144_AT_SingleBasePeaks.NA15-SRR5627145_AND_NA15-SRR5627144_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
bwprofilesNorm/WTCHG_538916_224192_VS_WTCHG_538916_220144_AT_SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
bwprofilesNorm/WTCHG_538916_224192_VS_WTCHG_538916_220144_AT_SingleBasePeaks.NA15-SRR5627145_AND_NA15-SRR5627144_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
'1) Zcwpw1 cotransfected with human Prdm9, at human Prdm9 sites' \
'2) Zcwpw1 cotransfected with human Prdm9, at chimp Prdm9 sites' \
'3) Zcwpw1 cotransfected with chimp Prdm9, at human Prdm9 sites' \
'4) Zcwpw1 cotransfected with chimp Prdm9, at chimp Prdm9 sites'
