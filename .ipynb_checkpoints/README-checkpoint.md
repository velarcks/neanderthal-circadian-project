Neanderthal Circadian Project
=============================
` Keila S. Velazquez-Arcelay `


Scripts
-------

bin/
* circadian_genes_candidate.py
  - Description: Imports 5 different sources of evidence to detect genes with circadian function. Generates a set of low, medium, and high confidence circadian genes. High confidence circadian genes are defined as genes that were expertly curated by Dr. Douglas McMahon's team or genes that contain evidence from 3 out of the other 4 sources. Medium confidence genes are defined as genes with evidence from 2 out of 4 sources. Low confidence genes contain evidence only in 1 source.
  - Output: circadian_genes_candidate.tsv

* circadian_gene_loci.py
  - Description: Create hg19 loci for each gene in the previously generated circadian set. Include only medium and high confidence genes. Also creates file with 1Mb gene flanking loci for each gene.
  - Output: circadian_genes.bed; data/circadian_genes_flanking.bed

* circadian_variants_gene.py
  - Description: Extract circadian variants inside gene regions from 1KGP variants dataset.
  - Output: circadian_variants_gene.bed

* circadian_tss_hg38.py
  - Description: Create promoter loci for each circadian gene based on the TSS. TSS were retrieved from Biomart. A promoter region is defined as 5kb up- and 1k downstream of the TSS. The output is in genome build hg38. Liftover to hg19 using: "liftOver circadian_promoter_hg38.bed /dors/capra_lab/data/ucsc/liftOver/hg38ToHg19.over.chain.gz circadian_promoter_hg38.liftover.hg19.bed circadian_promoter_hg38.liftover.hg19_unlifted.bed"
  - Output: circadian_promoter_hg38.bed

* circadian_variants_promoter.py
  - Description: Extract circadian variants in promoter regions from 1KGP dataset.
  - Output: circadian_variants_promoters.bed

* circadian_variants_ccres.py
  - Description: Extract variants flanking circadian genes from a set of cCREs (The ENCODE Project Consortium, 2020).
  - Output: circadian_variants_ccres.bed

* circadian_variants_introgressed.py
  - Description: Extract circadian variants from a dataset of archaic introgressed variants published by Browning et al, 2018. The variants are classified into: 
  - Output: circadian_variants_introgressed.bed

* parse_kuhlwilm19.py NOT_AVAILABLE
  - Description: Identifies human specific and archaic specific variants from a set of variants published by Kuhlwilm and Boeckx 2019. The estimation here doesn't exactly match the counts from the published data.
  - Output: kuhlwilm19_archaic_fixed.tsv; kuhlwilm19_human_fixed.tsv; kuhlwilm19_human_fixed_tony.tsv

* circadian_variants_fixed.py
  - Description: Extract circadian variants from human specific and archaic specific variants published in Kuhlwilm et al. 2019
  - Output: circadian_variants_hhmc.bed; circadian_variants_ahmc.bed

* circadian_variants_fixed_promoter.py
  - Description: Extract circadian variants in promoter regions using circadian promoter regions file.
  - Output: circadian_variants_hhmc_promoters.bed; circadian_variants_ahmc_promoters.bed

* circadian_variants_fixed_ccres.py
  - Description: Extract human specific and Neanderthal specific (Kuhlwilm et al. 2019) circadian variants in regulatory regions.
  - Output: circadian_variants_hhmc_ccres.bed; circadian_variants_ahmc_ccres.bed

* predixcan_define_dr.py
  - Description: Extracts list of divergently regulated genes by extracting genes that have a p-value = 0
  - Output: predixcan_dr_circadian_altai.tsv; predixcan_dr_circadian_vindija.tsv; predixcan_dr_circadian_denisova.tsv;  predixcan_dr_circadian.tsv; predixcan_dr_circadian_archaics.tsv; 
 
-------

* gtex_extract_eqtls.py
  - Description: Extract circadian and introgressed eQTLs from GTEx v8, hg19.
  - Output: gtex_v8_hg19_circadian.bed; gtex_v8_hg19_introgressed.bed
 
* gtex_introgressed_counts.py
  - Description: Export counts of unique eQTLs in circadian genes and GTEx tissues.
  - Output: gtex_introgressed_tissue_counts.tab; gtex_introgressed_gene_counts.tab


-------

statistics/
* enrichment_fixed_circadian_genes.py
     - Description: Testing for enrichment of archaic and human specific variants (Kuhlwilm et al. 2019) inside circadian genes using Fisher's exact test.
     - Output: enrichment_fixed_circadian_genes.txt

* enrichment_fixed_circadian_promoter.py
     - Description: Testing for enrichment of archaic and human specific variants (Kuhlwilm et al. 2019) in the promoter region of circadian genes using Fisher's exact test.
     - Output: enrichment_fixed_circadian_promoter.txt

* enrichment_fixed_circadian_ccres.py
     - Description: Testing for enrichment of archaic and human specific variants (Kuhlwilm et al. 2019) in the regulatory regions of circadian genes using Fisher's exact test.
     - Output: enrichment_fixed_circadian_ccres.txt

* enrichment_circadian_introgressed_eqtls.py
  - Description: Hypergeometric distribution of introgressed variants in circadian genes with evidence of being eQTL in GTEx. 
  - Output: enrichment_circadian_introgressed_eqtls.txt

* enrichment_predixcan_dr_analysis.py NOT_AVAILABLE
  - Description: 
  - Output: 

* enrichment_ai_genomatnn.py | enrichment_ai_maladapt.py | enrichment_ai_func.py
  - Description: Tests for overrepresentation of circadian variants in regions predicted to be under adaptive introgression by the genomatnn and MaLAdapt methods.
  - Output: enrichment_ai_genomatnn.txt; enrichment_ai_maladapt.txt



Notebooks
---------

* Fig3ab_upsetr_sav_dr.ipynb
  - Description: Figure 3, panels A and B. A) Archaic splice-altering variants in 4 archaic individuals. B) Archaic divergently regulated gene/tissue pairs in 3 archaic individuals.

* Fig3c_predixcan_dr_piechart.ipynb
  - Description: Figure 3, panel C. Proportion of circadian genes that contain evidence of regulatory differences between archaics and humans from splice-altering variants (SAV) and divergent regulation (DR) predictions.

* Fig4_predixcan_distribution.ipynb
  - Description: Figure 4. Distribution of PrediXcan core clock gene regulation predictions. Figure S2. Distribution of PrediXcan circadian gene regulation predictions in divergently regulated genes.

* Fig5b_fishers_introgressed_by_tissue.ipynb
  - Description: Figure 5, panel B. Testing for the probability of success of introgressed circadian eQTLs in each GTEx tissue using Fisher's exact test.

* Fig6_ukbiobank_morningness_1180.ipynb
  - Description: Figure 6. Plot the cumulative fraction of introgressed circadian loci that increase morningness propensity. The Morning/evening person trait data was extracted from the UKBiobank GWAS by the Neale Lab.

* Fig7b_opentargets_traits_per_snp.ipynb
  - Description: Figure 7, panel B. Plot the counts of genome-wide association traits per SNP in a set of introgressed circadian and introgressed non-circadian variants. The traits were scraped from Open Targets Genetics.

* FigS3_ai_hbars.ipynb
  - Description: Figure S3. Enrichment analysis to test the probability of success of adaptive ingrogression regions containing introgressed circadian variants.

