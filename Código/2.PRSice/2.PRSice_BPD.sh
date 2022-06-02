#########################################
#### PRSice: 0017_HYPERACTIVITY (HA) ####
#########################################
#a=/media/mireia/HardDiskHDD/Documentos/Universidad/Master_Bioinformatica_UOC/Curso_21_22/SEMESTRE_2/TFM/PEC3/GWAS_IMPUTATION/GWAS_GRANADA/1_GWAS_GR_HELIX_protocol/GWAS_PROCESSING/3_PRSs_PRSice/3.1_Inputs
#b=/media/mireia/HardDiskHDD/Documentos/Universidad/Master_Bioinformatica_UOC/Curso_21_22/SEMESTRE_2/TFM/PEC3/GWAS_IMPUTATION/GWAS_GRANADA/1_GWAS_GR_HELIX_protocol/GWAS_PROCESSING/3_PRSs_PRSice/3.2_Outputs


a=/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/GWAS_IMPUTATION/GWAS_GRANADA/1_GWAS_GR_HELIX_protocol/GWAS_PROCESSING/3_PRSs_PRSice/3.1_Inputs
b=/mnt/9e33c15a-acae-4445-96c2-5ebe80a1c36a/MIREIA/GWAS_IMPUTATION/GWAS_GRANADA/1_GWAS_GR_HELIX_protocol/GWAS_PROCESSING/3_PRSs_PRSice/3.2_Outputs

plink1.9 --bfile $a/PUBMEP_HRC.merged.sex.nohh.rs \
    --maf 0.05 \
    --mind 0.1 \
    --geno 0.1 \
    --hwe 1e-6 \
    --make-just-bim \
    --make-just-fam \
    --out $b/PUBMEP_HRC.merged.sex.nohh.rs.qc

## Addiong column of p-values (all 1) to enable PRSice to analyze our data using PGS Catalog PRS as bas data input.
# Add column of 1s
awk -F '\t' -v OFS='\t' '{ $(NF+1) = 0; print }' Obesity_PGS002033.txt  > Obesity_PGS002033_pvals.txt
# Change column names
sed -e '1s/effect_weight/BETA/' -e '1s/0/P-value/' Obesity_PGS002033_pvals.txt > Obesity_PGS002033_pvals_named.txt

Rscript $a/PRSice.R --dir . \
--prsice $a/PRSice_linux \
--base $a/Obesity_PGS002033_pvals_named.txt \
--target $a/PUBMEP_HRC.merged.sex.nohh.rs \
--thread 10 \
--snp rsID \
--chr chr_name \
--bp chr_position \
--A1 effect_allele \
--A2 other_allele \
--stat BETA \
--pvalue P-value \
--binary-target T \
--all-score \
--out $b/PUBMEP_BPD \


## NOTES:
# if there are spaces after " \ "next command will not be recogized
# --fastscore T \  : true if you want to get different thresholds
# --all-score : to get PRS for each of the thesholds indicated
# --out : indicate name of the output files
# by default: makes average, if --score sum indicated, sum of effect sizes
# We have to specify the names of the variables in the GWAS reference as:
    # --snp SNP \
    # --chr CHR \
    # --bp BP \
    # --A1 A1 \ #effect allele
    # --A2 A2 \ #non effect allele
    # --stat OR \
    # --pvalue P

# A temporal folder that PRSice needs to run is automatically and by default created
# Therefore, one can remove them by the following line of command
rm -r lib
