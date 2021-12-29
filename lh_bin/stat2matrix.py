import csv
from collections import defaultdict
from pandas import DataFrame
import sys

#https://stackoverflow.com/questions/3269769/convert-three-column-text-file-to-matrix

rdr = csv.reader(open(sys.argv[1]), delimiter='\t', skipinitialspace=True)
datacols = defaultdict(list)

# skip header
# rdr.next()
for qry, ref, result in rdr:
    datacols['qry'].append(qry)
    datacols['ref'].append(ref)
    datacols['result'].append(float(result))

df = DataFrame(datacols)
df2 = df.pivot(index='qry', columns='ref', values='result')
df2.to_csv(sys.argv[1] + '.xls', sep='\t')

