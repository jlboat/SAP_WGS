VCF="propinquum.all"

~/gatk-4.2.0.0/gatk --java-options "-Xmx50g" VariantFiltration \
    -V ${VCF}.vcf.gz \
    --filter-name "QD2" \
    --filter-expression "QD < 2.0" \
    --filter-name "F-NEG" \
    --filter-expression "InbreedingCoeff < 0.0" \
    --filter-name "QUAL30" \
    --filter-expression "QUAL < 30.0" \
    --filter-name "SOR3" \
    --filter-expression "SOR > 3.0" \
    --filter-name "FS60" \
    --filter-expression "FS > 60.0" \
    --filter-name "MQ40" \
    --filter-expression "MQ < 40.0" \
    --filter-name "MQRankSum-12.5" \
    --filter-expression "MQRankSum < -12.5" \
    --filter-name "ReadPosRankSum-8" \
    --filter-expression "ReadPosRankSum < -8.0" \
    -O ${VCF}.hard_filtered.vcf

bgzip -c --threads 16 ${VCF}.hard_filtered.vcf > ${VCF}.hard_filtered.vcf.gz
bash ~/scripts/indexVcf.bash -i ${VCF}.hard_filtered.vcf.gz -t 16

bcftools view --types snps --min-alleles 2 --max-alleles 2 --threads 16 -Oz -f .,PASS ${VCF}.hard_filtered.vcf.gz > ${VCF}.quality.vcf.gz
