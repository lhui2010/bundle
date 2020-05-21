import sys

flanking = 10000

with open(sys.argv[1]) as fh:
    for line in fh:
        (name, chrid, start, end) = line.rstrip().split()
        start = str(int(start) - flanking)
        end = str(int(end) + flanking)
        print("\t".join([name, chrid, start, end]))
