import sys
import pandas as pd
import itertools

table = pd.read_csv(sys.argv[1])

sp_list = list(table.species)

# for i in itertools.combinations(sp_list,2):
#     print(".".join(list(i)))

for i in range(0,len(sp_list)):
    for j in range(0,len(sp_list)):
        print("{}.{}".format(sp_list[i].replace(' ', '_'),sp_list[j].replace(' ', '_')))

