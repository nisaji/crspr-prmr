from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from joblib import Parallel, delayed
import subprocess
import re
import os
from pyfaidx import Fasta
import linecache
import sys
import shutil
import primer3
from concurrent.futures import ProcessPoolExecutor

"""
  usage : introncheck.py genome.fasta transdecoder.gff Label.reference.database.gff3

"""

genome = Fasta(sys.argv[1])
genome_gff = sys.argv[2]
gRNA_gff = sys.argv[3]

scaffold_list = []

intron_len = 1500
gRNA_length = 23
introncheck_tmp = "./introncheck_tmp/"

combined_exon_seq_list = []
exon_dict = {}
exon_seq_dict = {}
exon_with_gRNA_dict = {}
exon_loc_dict = {}
exon_list = []
mRNA_dict = {}

# global combined_exon_seq_list
# global exon_dict
# global exon_seq_dict
# global exon_with_gRNA_dict
# global exon_loc_dict
# global exon_list


def gff_split(genome_gff, gRNA_gff):
    # Make scaffold(first colum of gff) list
    genome_line = [line.rstrip() for line in open(genome_gff)]
    scaffold = ""
    for line in genome_line:
        seqname = line.split("\t")[0]
        if seqname not in scaffold_list:
            scaffold_list.append(seqname)

    # Split genome_gff file in introncheck_tmp
    for line in genome_line:
        seqname = line.split("\t")[0]
        if seqname != scaffold:
            scaffold = seqname
            f = open(introncheck_tmp
                     + sys.argv[2].split(".gff")[0] + seqname + ".gff3", mode='a')
            f.write(line + "\n")
        else:
            f.write(line + "\n")

    # Split gRNA_gff file in introncheck_tmp
    gRNA_line = [line.rstrip() for line in open(gRNA_gff)]
    gRNA_scaffold = ""
    for line in gRNA_line:
        seqname = line.split("\t")[0]
        if seqname != gRNA_scaffold:
            scaffold = seqname
            f = open(introncheck_tmp +
                     sys.argv[3].split(".gff")[0] + scaffold + ".gff3", mode='a')
            f.write(line + "\n")
        else:
            f.write(line + "\n")
    return()


# Check intron between 2 exons and make exon1:exon2 dictionaries
def intron_check(genome_line):
    for line in genome_line:
        colums = line.split("\t")
        seqname = colums[0]
        feature = colums[2]
        Start = colums[3]
        End = colums[4]
        attributes = colums[8]
        ParentID = attributes.split(';')[1].replace("Parent=", "")
        if feature == "intron":
            intron_seq = genome[seqname][int(Start):int(End)].seq
            # Search exon1
            if len(intron_seq) < int(intron_len):
                exon1 = ""
                exon2 = ""
                exonParent = ParentID
                exon1End = int(Start) - 1
                exon2Start = int(End) + 1
                for line in genome_line:
                    feature = line.split("\t")[2]
                    End = line.split("\t")[4]
                    Parent = line.split("\t")[8].split(
                        ';')[1].replace("Parent=", "")
                    if Parent == exonParent and End == str(exon1End) and feature == "exon":
                        exon1 = line.split("\t")[8].split(';')[
                            0].replace('ID=', '')
                # After searching exon1, search exon2
                for line in genome_line:
                    feature = line.split("\t")[2]
                    Start = line.split("\t")[3]
                    Parent = line.split("\t")[8].split(
                        ';')[1].replace("Parent=", "")
                    if Parent == exonParent and str(Start) == str(exon2Start) and feature == "exon":
                        exon2 = line.split("\t")[8].split(';')[
                            0].replace('ID=', '')
                exon_dict[exon1] = exon2
    return(exon_dict)


def select_exon_with_gRNA(exon_dict, scaffold, genome_line):
    # Create exon list wich include gRNA
    gRNA_line = [line.rstrip() for line in open(
        introncheck_tmp + sys.argv[3].split(".gff")[0] + scaffold + ".gff3")]
    for line in gRNA_line:
        exon = line.split("\t")[8].split("Parent=")[1]
        if exon not in exon_list:
            exon_list.append(exon)
    # Convert exon_list -> exon_set
    exon_set = set(exon_list)
    #  Search both exons in exon_dict are included in exon_list above(exon in exon_list inlude gRNA)
    for exon1, exon2 in exon_dict.items():
        if exon1 in exon_set and exon2 in exon_set:
            # check mRNA direction
            mRNA = exon1.split(".exon")[0]
            for line in genome_line:
                colums = line.split("\t")
                feature = colums[2]
                direction = colums[6]
                mRNAID = colums[8].split(";")[0].replace("ID=", "")
                if feature == "mRNA" and mRNAID == mRNA:
                    mRNA_direction = direction
            if mRNA not in exon_with_gRNA_dict:
                exon_with_gRNA_dict[mRNA] = []
            exon_with_gRNA_dict[mRNA].append(
                exon1 + ':' + exon2 + ":" + mRNA_direction)
    return(exon_with_gRNA_dict)


def prepare_combined_exon_seq_list(scaffold, exon_with_gRNA_dict, genome_line):
    exon_loc_dict = {}
    exon_seq_dict = {}
    seqname = scaffold
    # Create exon_loc_dict for all exons in scaffold
    for line in genome_line:
        colums = line.split("\t")
        feature = colums[2]
        if feature == "exon":
            exonID = colums[8].split(";")[0].replace("ID=", "")
            Start = colums[3]
            End = colums[4]
            exon_seq = genome[seqname][int(int(Start) - 1): int(End)].seq
            exon_seq_dict[exonID] = exon_seq
            exon_loc_dict[exonID] = Start + "," + End
    # Items in exon_with_gRNA_dict.items is as below.
    # 'STRG.843.1.p1': ['STRG.843.1.p1.exon2:STRG.843.1.p1.exon3:+', ... ,'STRG.842.1.p3': ['STRG.842.1.p3.exon2:STRG.842.1.p3.exon3:+']
    # Item_key = mRNA , item_value = exon_sets,
    for mRNA, exon_sets in exon_with_gRNA_dict.items():
        # Make mRNA_list and use it later for check which exon is primer made.
        if mRNA not in mRNA_dict:
            mRNA_dict[mRNA] = 0
        for exon_set in exon_sets:
            exon1 = exon_set.split(":")[0]
            exon2 = exon_set.split(":")[1]
            direction = exon_set.split(":")[2]
            exon1seq = exon_seq_dict[exon1]
            exon1start = exon_loc_dict[exon1].split(",")[0]
            exon2seq = exon_seq_dict[exon2]
            exon2start = exon_loc_dict[exon2].split(",")[0]
            # Check direction of exons and which exon is upstream one in exon_with_gRNA_dict.
            if direction == "+":
                if exon1start < exon2start:
                    combined_exon_seq = exon1seq + exon2seq
                    upstream_exon = exon1
                    upstream_exon_start = exon_loc_dict[exon1].split(",")[
                        0]
                    upstream_exon_end = exon_loc_dict[exon1].split(",")[1]
                    upstream_exon_seq = exon1seq
                    downstream_exon = exon2
                    downstream_exon_seq = exon2seq
                    downstream_exon_start = exon_loc_dict[exon2].split(",")[
                        0]
                    downstream_exon_end = exon_loc_dict[exon2].split(",")[
                        1]
                else:
                    combined_exon_seq = exon2seq + exon1seq
                    upstream_exon = exon2
                    upstream_exon_start = exon_loc_dict[exon2].split(",")[
                        0]
                    upstream_exon_end = exon_loc_dict[exon2].split(",")[1]
                    upstream_exon_seq = exon2seq
                    downstream_exon = exon1
                    downstream_exon_seq = exon1seq
                    downstream_exon_start = exon_loc_dict[exon1].split(",")[
                        0]
                    downstream_exon_end = exon_loc_dict[exon1].split(",")[
                        1]
            elif direction == "-":
                # If direction is "-", each exon_seq will be reverse complementally.
                if exon1start > exon2start:
                    combined_exon_seq = exon2seq + exon1seq
                    combined_exon_seq = str(Seq(
                        combined_exon_seq).reverse_complement())
                    upstream_exon = exon1
                    upstream_exon_start = exon_loc_dict[exon1].split(",")[
                        0]
                    upstream_exon_end = exon_loc_dict[exon1].split(",")[1]
                    upstream_exon_seq = exon1seq
                    upstream_exon_seq = str(Seq(
                        upstream_exon_seq).reverse_complement())
                    downstream_exon = exon2
                    downstream_exon_seq = exon2seq
                    downstream_exon_seq = str(Seq(
                        downstream_exon_seq).reverse_complement())
                    downstream_exon_start = exon_loc_dict[exon2].split(",")[
                        0]
                    downstream_exon_end = exon_loc_dict[exon2].split(",")[
                        1]
                else:
                    combined_exon_seq = exon1seq + exon2seq
                    combined_exon_seq = str(Seq(
                        combined_exon_seq).reverse_complement())
                    upstream_exon = exon2
                    upstream_exon_start = exon_loc_dict[exon2].split(",")[
                        0]
                    upstream_exon_end = exon_loc_dict[exon2].split(",")[1]
                    upstream_exon_seq = exon2seq
                    upstream_exon_seq = str(Seq(
                        upstream_exon_seq).reverse_complement())
                    downstream_exon = exon1
                    downstream_exon_seq = exon1seq
                    downstream_exon_seq = str(Seq(
                        downstream_exon_seq).reverse_complement())
                    downstream_exon_start = exon_loc_dict[exon1].split(",")[
                        0]
                    downstream_exon_end = exon_loc_dict[exon1].split(",")[
                        1]
            combined_exon_seq_list.append(direction + ":" + upstream_exon + ":" + upstream_exon_start + ": " + upstream_exon_end + ": " + upstream_exon_seq +
                                          ":" + downstream_exon + ":" + downstream_exon_start + ":" + downstream_exon_end + ":" + downstream_exon_seq + ":" + combined_exon_seq)
    return(combined_exon_seq_list)


def make_primers(combined_exon_seq_list, scaffold):
    gRNA_line = [line.rstrip() for line in open(introncheck_tmp +
                                                sys.argv[3].split(".gff")[0] + scaffold + ".gff3")]
    with open(introncheck_tmp +
              sys.argv[3].replace(".gff3", "") + scaffold + ".primer.gff3", mode='w') as f:
        for combined_exon_seq in combined_exon_seq_list:
            elems = combined_exon_seq.split(":")
            direction = elems[0]
            upstream_exon = elems[1]
            upstream_exon_start = elems[2]
            upstream_exon_end = elems[3]
            upstream_exon_seq = elems[4]
            downstream_exon = elems[5]
            downstream_exon_start = elems[6]
            downstream_exon_end = elems[7]
            downstream_exon_seq = elems[8]
            combined_exon_seq = elems[9]
            # Fecth mRNA name from exon name because it will be used for mRNA_dict to check wich mRNA is primer made
            mRNA = upstream_exon.split(".exon")[0]
            # Check downstream_gRNA first.
            # Count gRNA num in downstream_exon.
            downstream_gRNA_num = 0
            for line in gRNA_line:
                ParentID = line.split("\t")[8].split(";")[
                    1].replace("Parent=", "")
                if ParentID == downstream_exon:
                    downstream_gRNA_num += 1
            # For each downstream exon, select most middle gRNA in exon.
            # Internal exon is most small value in abs(left_len - rigth_len)
            # Left_len : left area from gRNA in exon , right_len : right_area from gRNA in exon.
            downstream_gRNA_check_num = 0
            downstream_gRNA = ""
            downstream_gRNA_start = 0
            downstream_gRNA_end = 0
            for line in gRNA_line:
                colums = line.split("\t")
                gRNAID = colums[8].split(";")[0].replace("ID=", "")
                exonID = colums[8].split(";")[1].replace("Parent=", "")
                if exonID == downstream_exon and downstream_gRNA_check_num == 0:
                    downstream_gRNA_check_num += 1
                    downstream_gRNA_start_tmp = colums[3]
                    downstream_gRNA_end_tmp = colums[4]
                    downstream_gRNA_start = downstream_gRNA_start_tmp
                    downstream_gRNA_end = downstream_gRNA_end_tmp
                    downstream_left_len = int(
                        downstream_gRNA_start) - int(downstream_exon_start)
                    downstream_right_len = int(
                        downstream_exon_end) - int(downstream_gRNA_end)
                    downstream_len_diff = abs(
                        downstream_left_len - downstream_right_len)
                    downstream_gRNA = gRNAID
                elif exonID == downstream_exon and downstream_gRNA_check_num < downstream_gRNA_num:
                    downstream_gRNA_check_num += 1
                    downstream_gRNA_start_tmp = colums[3]
                    downstream_gRNA_end_tmp = colums[4]
                    downstream_left_len_tmp = int(downstream_gRNA_start_tmp) - \
                        int(downstream_exon_start)
                    downstream_right_len_tmp = int(
                        downstream_exon_end) - int(downstream_gRNA_end_tmp)
                    downstream_len_diff_tmp = abs(
                        downstream_left_len_tmp - downstream_right_len_tmp)
                    if downstream_len_diff > downstream_len_diff_tmp:
                        downstream_len_diff = downstream_len_diff_tmp
                        downstream_gRNA = gRNAID
                        downstream_gRNA_start = downstream_gRNA_start_tmp
                        downstream_gRNA_end = downstream_gRNA_end_tmp
                elif downstream_gRNA_check_num == downstream_gRNA_num:
                    # print("break")
                    break

            # Select primers in downstream_exon with primer3
            downstream_left_start = 0
            downstream_left_length = int(
                downstream_gRNA_start) - int(downstream_exon_start)
            downstream_right_start = int(downstream_gRNA_end) + 1
            downstream_right_length = int(
                downstream_exon_end) - int(downstream_right_start)
            if (downstream_left_length + gRNA_length + downstream_right_length) > 82 and downstream_left_length > 30 and downstream_right_length > 30:
                try:
                    downstream_primer_out = primer3.designPrimers({
                        'SEQUENCE_ID': downstream_exon,
                        'SEQUENCE_TEMPLATE': downstream_exon_seq,
                    },
                        {
                        "PRIMER_OPT_SIZE": 24,
                        "PRIMER_MIN_SIZE": 20,
                        "PRIMER_MAX_SIZE": 27,
                        "PRIMER_PRODUCT_SIZE_RANGE": [28, 200],
                        "PRIMER_OPT_TM": 65,
                        "PRIMER_MIN_TM": 60,
                        "PRIMER_MAX_TM": 72,
                        "PRIMER_FIRST_BASE_INDEX": 1,
                        "PRIMER_MAX_DIFF_TM": 2,
                        "PRIMER_MIN_GC": 45,
                        "PRIMER_OPT_GC_PERCENT": 50,
                        "PRIMER_MAX_GC": 65,
                        "PRIMER_NUM_RETURN": 5,
                        "SEQUENCE_PRIMER_REVCOMP": '',
                        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [downstream_left_start, downstream_left_length, downstream_right_start, downstream_right_len],
                    })
                    # downstream_exon_left_primer = downstream_primer_out["PRIMER_LEFT_0"]
                    downstream_exon_left_primer_seq = downstream_primer_out["PRIMER_LEFT_0_SEQUENCE"]

                    # downstream_exon_right_primer = downstream_primer_out["PRIMER_RIGHT_0"]
                    downstream_exon_right_primer_seq = downstream_primer_out[
                        "PRIMER_RIGHT_0_SEQUENCE"]
                    downstream_exon_product_size = downstream_primer_out['PRIMER_PAIR_0_PRODUCT_SIZE']
                    downstream_gRNA_gff = str(scaffold + "\t" + "CRISPR-Local" + "\t" + "gRNA" + "\t" + downstream_gRNA_start + "\t" + downstream_gRNA_end + "\t" + "." + "\t" + direction + "\t" + "." + "\t" + "ID=" + downstream_gRNA +
                                              ";" + "Parent=" + downstream_exon + ";" + "DownStreamExonLeftPrimer=" + downstream_exon_left_primer_seq + ";" "DownStreamExonRightPrimer=" + downstream_exon_right_primer_seq)
                    f.write(downstream_gRNA_gff + "\n")
                except KeyError:
                    continue
                # After making primer on downstream_exon, make primer on upstream_exon
                # Count gRNA on upstream_exon
                upstream_gRNA_num = 0
                for line in gRNA_line:
                    ParentID = line.split("\t")[8].split(";")[
                        1].replace("Parent=", "")
                    if ParentID == upstream_exon:
                        upstream_gRNA_num += 1
                # Search suitable gRNA on upstream_exon
                upstream_gRNA_check_num = 0
                upstream_gRNA = ""
                for line in gRNA_line:
                    colums = line.split("\t")
                    gRNAID = colums[8].split(";")[0].replace("ID=", "")
                    exonID = colums[8].split(";")[1].replace("Parent=", "")
                    if exonID == upstream_exon and upstream_gRNA_check_num == 0:
                        upstream_gRNA_check_num += 1
                        upstream_gRNA_start_tmp = colums[3]
                        upstream_gRNA_end_tmp = colums[4]
                        upstream_gRNA_start = upstream_gRNA_start_tmp
                        upstream_gRNA_end = upstream_gRNA_end_tmp
                        upstream_left_length = int(
                            upstream_gRNA_start) - int(upstream_exon_start)
                        upstream_gRNA = gRNAID
                    elif exonID == upstream_exon and upstream_gRNA_check_num < upstream_gRNA_num:
                        upstream_gRNA_check_num += 1
                        upstream_gRNA_start_tmp = colums[3]
                        upstream_gRNA_end_tmp = colums[4]
                        upstream_left_length_tmp = int(upstream_gRNA_start) - \
                            int(upstream_exon_start)
                        if upstream_left_length < upstream_left_length_tmp:
                            upstream_left_length = upstream_left_length_tmp
                            upstream_gRNA = gRNAID
                            upstream_gRNA_start = upstream_gRNA_start_tmp
                            upstream_gRNA_end = upstream_gRNA_end_tmp
                    elif upstream_gRNA_check_num == upstream_gRNA_num:
                        # print("break")
                        break
                try:
                    combined_min_product_size = int(
                        downstream_exon_product_size * 1.2)
                    combined_max_product_size = int(
                        downstream_exon_product_size * 1.5)
                    combined_exon_length = len(combined_exon_seq)
                    downstream_right_primer_region_start = (int(upstream_exon_end) - int(
                        upstream_exon_start)) + (int(downstream_gRNA_end) - int(downstream_exon_start))
                    downstream_right_primer_region_length = combined_exon_length - \
                        downstream_right_primer_region_start
                    upstream_left_primer_start = 0
                    upstream_primer_out = primer3.designPrimers({
                        'SEQUENCE_ID': upstream_exon,
                        'SEQUENCE_TEMPLATE': combined_exon_seq,
                    },
                        {
                        "PRIMER_OPT_SIZE": 24,
                        "PRIMER_MIN_SIZE": 20,
                        "PRIMER_MAX_SIZE": 27,
                        "PRIMER_PRODUCT_SIZE_RANGE": [combined_min_product_size, combined_max_product_size],
                        "PRIMER_OPT_TM": 65,
                        "PRIMER_MIN_TM": 60,
                        "PRIMER_MAX_TM": 72,
                        "PRIMER_FIRST_BASE_INDEX": 1,
                        "PRIMER_MAX_DIFF_TM": 2,
                        "PRIMER_MIN_GC": 45,
                        "PRIMER_OPT_GC_PERCENT": 50,
                        "PRIMER_MAX_GC": 65,
                        "PRIMER_NUM_RETURN": 5,
                        "SEQUENCE_PRIMER_REVCOMP": downstream_exon_right_primer_seq,
                        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [int(upstream_left_primer_start), int(upstream_left_length), int(downstream_right_primer_region_start), int(downstream_right_primer_region_length)],
                    })
                    upstream_exon_left_primer_seq = upstream_primer_out["PRIMER_LEFT_0_SEQUENCE"]
                    upstream_gRNA_gff = str(scaffold + "\t" + "CRISPR-Local" + "\t" + "gRNA" + "\t" + upstream_gRNA_start + "\t" + upstream_gRNA_end + "\t" + "." + "\t" + direction + "\t" + "." + "\t" + "ID=" +
                                            upstream_gRNA + ";" + "Parent=" + upstream_exon + ";" + "UpStreamExonLeftPrimer=" + upstream_exon_left_primer_seq + ";" "DownStreamExonRightPrimer=" + downstream_exon_right_primer_seq)
                    f.write(upstream_gRNA_gff + "\n")
                    mRNA_dict[mRNA] = 1
                except KeyError:
                    continue


def clean_make_primer(scaffold):
    with open(introncheck_tmp +
              sys.argv[3].replace(".gff3", "") + scaffold + ".primer.clean.gff3", mode='w') as f:
        gRNA_line = [line.rstrip() for line in open(
            introncheck_tmp + sys.argv[3].replace(".gff3", "") + scaffold + ".primer.gff3")]
        for line in gRNA_line:
            line_num = gRNA_line.index(line)
            colums = line.split("\t")
            attributes = colums[8]
            if "Down" in attributes:
                downstream_exon = attributes.split(
                    ";")[1].replace("Parent=", "")
                downstream_exon_num = downstream_exon.split("exon")[1]
                upstream_exon = downstream_exon.split(
                    'exon')[0] + "exon" + str(int(downstream_exon_num) - 1)
                next_line = gRNA_line[line_num + 1]
                next_line_colums = next_line.split("\t")
                attributes = next_line_colums[8]
                if "Up" in attributes and upstream_exon in attributes:
                    f.write(line + "\n" + next_line + "\n")


def process_scaffold(scaffold):
    # This function is series of processes to each scaffolds
    global combined_exon_seq_list
    global exon_dict
    global exon_seq_dict
    global exon_with_gRNA_dict
    global exon_loc_dict
    global exon_list
    global mRNA_dict
    print(scaffold)
    genome_line = [line.rstrip() for line in open(introncheck_tmp +
                                                  sys.argv[2].split(".gff")[0] + str(scaffold) + ".gff3")]
    intron_check(genome_line)
    select_exon_with_gRNA(exon_dict, scaffold, genome_line)
    prepare_combined_exon_seq_list(
        scaffold, exon_with_gRNA_dict, genome_line)
    make_primers(combined_exon_seq_list, scaffold)

    primer_set_rate = str(
        sum(value == 0 for value in mRNA_dict.values()) / len(mRNA_dict))
    print("%s primers set rate:%s" % (scaffold, primer_set_rate))
    return()


def main():
    if os.path.exists(introncheck_tmp) != True:
        os.mkdir(introncheck_tmp)
    else:
        shutil.rmtree(introncheck_tmp)
        os.mkdir(introncheck_tmp)
    print("spliting gff files...")
    gff_split(genome_gff, gRNA_gff)
    # print(scaffold_list)
    # Multi thread process
    with ProcessPoolExecutor(max_workers=7) as process:
        processes = process.map(process_scaffold, scaffold_list)
        for process in processes:
            process
    # Merge all gRNA/primers  set in one *.primer.gff3 and remove all intron_tmp files
    # merged_output = sys.argv[3].split(".gff")[0] + ".primers.gff"
    # print("merging each gRNA output")
    # subprocess.call("cat *.primers > %s" % merged_output)
    # subprocess.call("mv %s ../introncheck_out/%s" %
    #                 (merged_output, merged_output))
    # shutil.rmtree(introncheck_tmp)


if __name__ == "__main__":
    main()
