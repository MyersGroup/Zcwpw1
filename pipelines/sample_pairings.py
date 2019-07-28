
# sample, input, description
# when there are two biological replicates, seperate them using the special string '_AND_'
# when there are pseudoreplicates, just enter the common sample name, assumes the PR1 PR2 file ending convention was used


sample_pairings = [("WTCHG_538916_221156", "WTCHG_538916_217108","P1 ChipHA_ZHA vs In_ZHA"),
              ("WTCHG_538916_223180", "WTCHG_538916_220144","P1 ChipHA_ZHA_hP9V5 vs In_ZHA_hP9V5"),
              ("WTCHG_538916_224192", "WTCHG_538916_220144","P1 ChipHA_ZHA_cP9V5 vs In_ZHA_hP9V5"),
              
              # combined inputs
              ("WTCHG_538916_221156", "WTCHG_538916_217108_AND_WTCHG_538916_220144","P1 ChipHA_ZHA vs In_ZHA+In_ZHA_hP9V5"),
              ("WTCHG_538916_223180", "WTCHG_538916_217108_AND_WTCHG_538916_220144","P1 ChipHA_ZHA_hP9V5 vs In_ZHA+In_ZHA_hP9V5"),
              ("WTCHG_538916_224192", "WTCHG_538916_217108_AND_WTCHG_538916_220144","P1 ChipHA_ZHA_cP9V5 vs In_ZHA+In_ZHA_hP9V5"),
              
              # chip vs chip
              ("WTCHG_538916_223180", "WTCHG_538916_221156","P1 ChipHA_ZHA_hP9V5 vs ChipHA_ZHA"),
              ("WTCHG_538916_224192", "WTCHG_538916_221156","P1 ChipHA_ZHA_cP9V5 vs ChipHA_ZHA"),


              # Human
              ("NA15-SRR5627138_AND_NA15-SRR5627139", "NA15-SRR5627140","NA ChipGFP_YFPhP9_rep1+2 vs In_YFPhP9"),

              ("NA15-SRR5627141", "NA15-SRR5627140","NA ChipH3K4_YFPhP9 vs In_YFPhP9"),
              #("NA15-SRR5627141", "NA15-SRR5627136_AND_NA15-SRR5627135","NA ChipH3K4_YFPhP9 vs ChipH3K4_UT"),

              ("NA15-SRR5627146", "NA15-SRR5627143","NA ChipHA_hP9HA vs In_hP9HA"),
              ("NA15-SRR5627147", "NA15-SRR5627143","NA ChipHA_hP9V5 vs In_hP9HA"),
              ("NA15-SRR5627146_AND_NA15-SRR5627147", "NA15-SRR5627143","NA ChipHA_hP9HA+V5 vs In_hP9HA"),

              # Chimp
              ("NA15-SRR5627144", "NA15-SRR5627143","NA ChipHA_cP9HA vs In_hP9HA"),
              ("NA15-SRR5627145", "NA15-SRR5627143","NA ChipHA_cP9V5 vs In_hP9HA"),
              ("NA15-SRR5627145_AND_NA15-SRR5627144", "NA15-SRR5627143","NA ChipHA_cP9HA+V5 vs In_hP9HA"),

              # Histone Mods
              ## H3K4
              ("NA15-SRR5627152_AND_NA15-SRR5627153", "NA15-SRR5627143","NA ChipH3K4_hP9HA_rep1+2 vs In_hP9HA"),
              ("NA15-SRR5627150_AND_NA15-SRR5627151", "NA15-SRR5627142","NA ChipH3K4_UT_rep1+2 vs In_UT"),
              ("NA15-SRR5627152_AND_NA15-SRR5627153", "NA15-SRR5627150_AND_NA15-SRR5627151","NA ChipH3K4_hP9HA_rep1+2 vs ChipH3K4_UT_rep1+2"),
              
              ## H3K36
              ("NA15-SRR5627149", "NA15-SRR5627143","NA ChipH3K36_hP9HA vs In_hP9HA"),
              ("NA15-SRR5627148", "NA15-SRR5627142","NA ChipH3K36_UT vs In_UT"),
              ("NA15-SRR5627149", "NA15-SRR5627148","NA ChipH3K36_hP9HA vs ChipH3K36_UT"),


              # Project2 Chip H3K4 & ZF only
              ("WTCHG_659233_256196", "WTCHG_659233_255184","P2 ChipV5_hP9V5zfonly vs In_hP9V5zfonly"),

              ("WTCHG_659233_253160", "WTCHG_659233_252148","P2 ChipHA_ZHA_cP9V5 vs In_ZHA_cP9V5"),
              ("WTCHG_659233_254172", "WTCHG_659233_252148","P2 ChipH3K4_ZHA_cP9V5 vs In_ZHA_cP9V5"),

              ("WTCHG_659233_250124", "WTCHG_659233_249112","P2 ChipHA_ZHA_hP9V5 vs In_ZHA_hP9V5 (Repeat)"),
              ("WTCHG_659233_251136", "WTCHG_659233_249112","P2 ChipH3K4_ZHA_hP9V5 vs In_ZHA_hP9V5"),

              ("WTCHG_659233_234122", "WTCHG_659233_233110","P2 ChipHA_ZHA_hP9V5zfonly vs In_ZHA_hP9V5zfonly"), 
              ("WTCHG_659233_235134", "WTCHG_659233_233110","P2 ChipH3K4_ZHA_hP9V5zfonly vs In_ZHA_hP9V5zfonly")] 
