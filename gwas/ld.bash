#!/bin/bash

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize variables:
INPUT=""
OUTPUT=""

while getopts "h?i:o:" opt; do
    case "${opt}" in
    h|\?)
        printf "\nEx. bash ld.bash -i snps.txt -o output_dir n\n\
                -i      input SNP file -- one per line (required)\n\
                -o      output directory for files (required)\n\
";
        exit 0;
        ;;
    i)  INPUT=${OPTARG}
        ;;
    o)  OUTPUT=${OPTARG}
        ;;
    esac
done
shift $((OPTIND -1))

if [ "${INPUT}" == "" ]
  then
    printf "\nbash ld.bash -i snps.txt -o output_dir    -i REQUIRED\n\n"
    exit 0;
fi

if [ "${OUTPUT}" == "" ]
then
    printf "\nbash ld.bash -i snps.txt -o output_dir    -o REQUIRED\n\n"
    exit 0;
fi

if [ ! -e ${OUTPUT} ]
then
    mkdir -p ${OUTPUT}
fi

for i in $(seq 1 $(wc -l ${INPUT} | cut -f 1 -d ' '))
do
    SNP_NAME=$(head -n ${i} ${INPUT} | tail -n 1)
    echo $SNP_NAME
    ~/Plink/plink --bfile named_data \
        --r2 --ld-snp ${SNP_NAME} --ld-window-kb 10000 --ld-window 99999 --ld-window-r2 0.15
    mv plink.ld ./${OUTPUT}/${SNP_NAME}.plink.ld
done
