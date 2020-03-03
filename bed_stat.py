#!/usr/bin/env python
import sys
import pandas as pd
import logging

#http://www.datasciencemadesimple.com/descriptive-summary-statistics-python-pandas/

size_list = []
with open(sys.argv[1]) as fh:
    for line in fh:
        mylist = line.rstrip().split()
        size_list.append(int(mylist[2]) - int(mylist[1]))

df = pd.DataFrame(size_list)

logging.warning('Finished Reading df')

pd.set_option('display.float_format', lambda x: '%.0f' % x)

print(df.describe())

