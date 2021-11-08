import sys
import gzip
import argparse

def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description="This script was " + 
            "designed to generate a table of variant calls from a VCF.\n\n")

    parser.add_argument("--vcf", type=str, required=True, \
    help="Input VCF file.\n",
    action="store")

    parser.add_argument("--output", type=str, required=True, \
    help="Additional output column name -- recommend chromosome name or population name", 
    action="store")

    return parser.parse_args()


def check_filename(filename):
    if filename.endswith("gz"):
        return "gzvcf"
    elif filename.endswith("vcf"):
        return "vcf"
    else:
        sys.write.err("Unrecognized VCF format\n")
        sys.exit(1)


def vcf_to_cmplot(line):
    split_line = line.split()
    chrom = split_line[0].lower().replace("chr","")
    if chrom.startswith("0"):
        chrom = chrom.replace("0","")
    position = split_line[1]
    return f"Chr{chrom}_{position}\t{chrom}\t{position}\t1"


if __name__ == "__main__":
    args = parse_arguments()
    file_type = check_filename(args.vcf)
    output_extra_field = args.output
    print(f"SNP\tchr\tpos\t{output_extra_field}")
    if file_type == "vcf":
        with open(args.vcf) as f:
            for line in f:
                if (not line.startswith("#")) and ("super" not in line.lower()):
                    print(vcf_to_cmplot(line))
    elif file_type == "gzvcf":
        with gzip.open(args.vcf, mode='rt') as f:
            for line in f:
                if (not line.startswith("#")) and ("super" not in line.lower()):
                    print(vcf_to_cmplot(line))
