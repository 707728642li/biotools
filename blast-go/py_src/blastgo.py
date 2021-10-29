# 2020-03-20 Li
# For python3 
# Usage: python this.py blast_result uniprot_sprot.dat.gz output_name
# uniprot_sprot.dat.gz is originally gz file!!!

import sys
import gzip
import collections

def get_blast_item(blast_result):
    blast_dict = {}  
    with open(blast_result, 'r') as f:
        for line in f:
            each = line.strip().split("\t")
            try:
                blast_dict[each[0]] = each[1].split("|")[1:]
            except:
                pass
        return blast_dict

def get_annotation_item(interpro_database):
    n=0
    info_dict = collections.defaultdict(dict)
    with gzip.open(interpro_database,'rt') as f:
        for line in f:
            index1 = line.strip().split()[0]
            if index1 == 'ID':
                try:
                    name = line.strip().split()[1]
                except:
                    name = "unexception"+str(n)                
            if index1 == 'DE':
                info_dict[name]["DE"] =info_dict[name].get("DE",[]) +[i for i in line.strip().split()[2:] if "ECO" not in i]
            try:
                if line.strip().split()[1]=="GO;":
                    info_dict[name]["GO;"] = info_dict[name].get("GO;",[])+[line.strip().split()[2]]
            except:
                pass
            try:
                if line.strip().split()[1]=="InterPro;":
                    info_dict[name]["InterPro;"] = info_dict[name].get("InterPro;",[])+[line.strip().split()[-1]]
            except:
                pass
            try:
                if line.strip().split()[1]=="KEGG;":
                    info_dict[name]["atID"] =info_dict[name].get("atID",[])+[line.strip().split()[2].split(":")[-1].replace(";","")]
            except:
                pass
            n+=1
            if n%10000==0:
                print("\r %i lines have been processed!"%n,end = "")
        return info_dict


in_1 = sys.argv[1]
in_2 = sys.argv[2]
out_put = sys.argv[3]

blast_res = get_blast_item(in_1)
info_res = get_annotation_item(in_2)

with open(out_put,"w") as result:
    for key in blast_res.keys():
        try:
            if info_res[blast_res[key][-1]]!={}:
                gene_id = key
                protein_id = "|".join(blast_res[key]+info_res[blast_res[key][-1]]["atID"])
                domain = "|".join(info_res[blast_res[key][-1]]["InterPro;"])
                description = "Full_name="+" ".join(info_res[blast_res[key][-1]]["DE"]).replace("Full=","")
                GO =  " ".join(info_res[blast_res[key][-1]]["GO;"])
                result.write(gene_id+"\t"+protein_id+"\t"+domain+"\t"+description+"\t"+GO+"\n")
        except:
            pass

print("\nFinished! file: {}".format(out_put))

#FINISH

