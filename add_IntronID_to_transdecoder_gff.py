import re

# lines = [line.rstrip() for line in open(
#     './Label.reference.database.unique.sort.sort.gff3')]
lines = [line.rstrip() for line in open(
    './transcripts.fasta.transdecoder.genome.introns.gff3')]
# lines = [line.rstrip() for line in open(
#     './Label.reference.database.unique.sort.gff3')]

# sgRNA = "none"

Parent_ID_dict = {}
num = 1
for line in lines:
    colums = line.split("\t")
    seqname = colums[0]
    Start = colums[3]
    End = colums[4]
    direction = colums[6]
    attributes = colums[8]
    intron_ID = attributes.split('=')[2]
    Parent_ID = attributes.split(';')[1]

    if Parent_ID not in Parent_ID_dict:
        Parent_ID_dict[Parent_ID] = 1
    else:
        Parent_ID_dict[Parent_ID] += 1
    print(seqname + "\t" + "transdecoder" + "\t" + "intron" + "\t" + str(Start) + "\t" +
          str(End) + "\t" + "." + "\t" + direction + "\t" + "." + "\t" + "ID=" + intron_ID + ".intron" + str(Parent_ID_dict[Parent_ID]) + ";" + Parent_ID)
