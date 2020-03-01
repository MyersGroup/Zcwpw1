
rule all:
    input:
        "results/Domain_Architechture.md",
        "results/sanger.md",
        "results/dmc1_counts.md"


rule domain_architechture:
  input:
    "data/alignment/Simons_Zcwpw1_alignment.rds",
    "data/alignment/simons_alignment.fa",
    "data/alignment/JSD_simons_alignment.txt",
    rmd = "analysis/Domain_Architechture.Rmd"
  output:
    md = "results/Domain_Architechture.md"
  shell:
    """
    R -e "knitr::knit('{input.rmd}', '{output.md}')"
    """


rule sanger:
  input:
    "data/cellular/WT_allele_ZCW_read_as_is.ab1",
    "data/cellular/KO_allele_ZCW_read_as_reverse_complement.ab1",
    rmd = "analysis/sanger.Rmd"
  output:
    md = "results/sanger.md"
  shell:
    """
    R -e "knitr::knit('{input.rmd}', '{output.md}')"
    """


rule dmc1_counts:
  input:
    "data/cellular/zcw_dmc1.xlsx",
    "data/cellular/Second prdm9 KO DMC1 foci.xlsx",
    "data/cellular/Rad51.xlsx",
    rmd = "analysis/dmc1_counts.Rmd"
  output:
    md = "results/dmc1_counts.md"
  shell:
    """
    R -e "knitr::knit('{input.rmd}', '{output.md}')"
    """
