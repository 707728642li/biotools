# 2020-03-27 Li
# For python3 
# Usage: python this.py obo_file output_name

import sys
obo = sys.argv[1]
out = sys.argv[2]

go_dict = {}    
for line in open(obo):
    try:
        each = line.strip().split(": ",1)  
        if each[0]=="id":
            name = each[1]
        elif each[0]=="name":
            go_dict[name] = each[1]
        elif each[0]=="namespace":
            go_dict[name] = [each[1][0].upper()]+[go_dict.get(name,"")]
    except:
        pass
with open(out,"w") as result:
    for key in go_dict.keys():
        result.write(key+'\t'+":".join(go_dict[key])[:50][::-1].split(" ",1)[-1][::-1]+"\n")#[:50] to avoid long sentences

# print("finished! File: {}".format(out))
