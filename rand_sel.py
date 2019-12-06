#!/usr/bin/env python
from random import randint
import sys

select_number = int(sys.argv[1])
filename = sys.argv[2]


def get_rand_list(total_size):
    num_dict = dict()
    while len(num_dict) < total_size:
        num_dict[randint(0, line_count-1)] = 1
    return num_dict.keys()


with open(filename) as fh:
    fh_list = fh.readlines()
    line_count = len(fh_list)
    line_list = get_rand_list(select_number)
    for line_id in line_list:
        print(fh_list[line_id], end = "")
