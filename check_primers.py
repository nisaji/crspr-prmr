import sys
from pyfaidx import Fasta
import re
import primer3
from Bio import SeqIO
from Bio.Seq import Seq


primers_gff = sys.argv[1]
genome_fasta = Fasta(sys.argv[2])
genome_gff = sys.argv[3]

primers_gff_lines = [line.rstrip() for line in open(primers_gff)]
genome_gff_lines = [line.rstrip() for line in open(genome_gff)]


def extract_gene_seq(line):
    global gene_seq
    # print(line)
    direction = line.split("\t")[6]
    attributes = line.split("\t")[8]
    gene_id = attributes.split(";")[0].split(".exon")[0].replace("ID=", "")
    gene_id = "ID=" + gene_id + ";"
    for line in genome_gff_lines:
        if gene_id in line:
            colums = line.split("\t")
            scaffold = colums[0]
            gene_start = colums[3]
            gene_end = colums[4]
            gene_seq = genome_fasta[scaffold][int(
                int(gene_start) - 1): int(gene_end)].seq
            if direction == "-":
                gene_seq = str(Seq(gene_seq).reverse_complement())
    return(gene_seq)


def extract_primers_and_gRNA(line):
    global left_primer_seq
    global right_primer_seq
    global gRNA_seq
    global product_loc
    global gRNA_id
    # print(line)
    colums = line.split("\t")
    scaffold = colums[0]
    gRNA_start = colums[3]
    gRNA_end = colums[4]
    attributes = colums[8]
    gRNA_id = attributes.split(";")[0].replace("ID=", "")
    if "DownStreamExonLeftPrimer" in attributes and "DownStreamExonRightPrimer" in attributes:
        product_loc = "down"
    elif "UpStreamExonLeftPrimer" in attributes and "DownStreamExonRightPrimer" in attributes:
        product_loc = "up"
    gRNA_seq = genome_fasta[scaffold][int(
        int(gRNA_start) - 1): int(gRNA_end)].seq
    left_primer_seq = re.split('[=;]', attributes)[5]
    right_primer_seq = re.split('[=;]', attributes)[7]
    return(gRNA_seq, left_primer_seq, right_primer_seq, product_loc, gRNA_id)


def check_product_size(gene_seq, gRNA_seq, left_primer_seq, right_primer_seq):
    global product_size
    pcr_out = primer3.designPrimers({
        'SEQUENCE_ID': "ID",
        'SEQUENCE_TEMPLATE': gene_seq,
    },
        {
        "PRIMER_OPT_SIZE": 24,
        "PRIMER_MIN_SIZE": 20,
        "PRIMER_MAX_SIZE": 27,
        "PRIMER_OPT_TM": 65,
        "PRIMER_MIN_TM": 60,
        "PRIMER_MAX_TM": 72,
        "PRIMER_FIRST_BASE_INDEX": 1,
        "PRIMER_MAX_DIFF_TM": 2,
        "PRIMER_MIN_GC": 45,
        "PRIMER_OPT_GC_PERCENT": 50,
        "PRIMER_MAX_GC": 65,
        "PRIMER_NUM_RETURN": 5,
        "SEQUENCE_PRIMER": left_primer_seq,
        "SEQUENCE_PRIMER_REVCOMP": right_primer_seq,
    })
    product_size = pcr_out['PRIMER_PAIR_0_PRODUCT_SIZE']
    return(product_size)


def main():
    for line in primers_gff_lines:
        # extract_primers_and_gRNA(line)
        extract_gene_seq(line)
        extract_primers_and_gRNA(line)
        check_product_size(gene_seq, gRNA_seq,
                           left_primer_seq, right_primer_seq)
        print(gRNA_id + ":" + product_loc + ":" + str(product_size))


if __name__ == "__main__":
    main()
