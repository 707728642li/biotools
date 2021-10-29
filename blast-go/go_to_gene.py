# 2020-03-28 Li
# For python3 
# Usage: python this.py anno_file out_file
# Delimiter="\t"

import sys

anno_file = sys.argv[1]
out_file = sys.argv[2]

with open(out_file,"w") as out:
    for line in open(anno_file):
        try:
            gene_name = line.strip().split("\t")[0]
            for each in line.strip().split("\t")[4].split(";"):
                if each!="":
                    out.write(each.replace(" ","")+"\t"+gene_name+"\n")
        except:
            pass

# print("Finished! File: {}".format(out_file))
