import re
import sys
import subprocess

sgRNA_database = sys.argv[1]
exon_set = set()
out_gff = sgRNA_database.split(".txt")[0] + ".gff"

lines = [line.rstrip() for line in open(sgRNA_database)]


# sgRNAs with no offtarget
def crispr_local_to_gff(lines):
    with open(out_gff, mode="w") as f:
        num = 1
        for line in lines:
            direction = "-"
            chromStart = 0
            chromEnd = 0
            colums = line.split("\t")
            try:
                loc = re.search("STRG\.\d*\.\d*\.p\d*.exon\d*",
                                colums[8]).group()
            except AttributeError:
                # print(colums[8])
                break

            if loc not in exon_set:
                exon_set.add(loc)
                num = 1
            else:
                num += 1
            loci = colums[1]
            seqname = loci.split(":")[0]
            chromStrand = loci.split(":")[1]

            if "+" in chromStrand:
                chromStart = int(chromStrand.replace("+", "").replace("-", ""))
                chromEnd = int(chromStart) + 22
                direction = "+"
            elif "-" in chromStrand:
                chromStart = int(chromStrand.replace(
                    "+", "").replace("-", ""))-1
                chromEnd = int(chromStart) + 22
                direction = "-"
            out_line = seqname + "\t" + "CRISPR-Local" + "\t" + "sgRNA" + "\t" + str(chromStart) + "\t" + str(
                chromEnd) + "\t" + "." + "\t" + direction + "\t" + "." + "\t" + "ID=" + loc + ".sgRNA" + str(num) + ";" + "Parent=" + loc
            f.write(out_line + "\n")


def main():
    crispr_local_to_gff(lines)


if __name__ == "__main__":
    main()
