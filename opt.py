import argparse

parser = argparse.ArgumentParser(description='Liftover from aligned.bed and old gff')
parser.add_argument('GFF', type=str, nargs = 1,
                    help='old.gff')
parser.add_argument('BED', type=str, nargs = 1,
                    help='aligned.bed')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')


args = parser.parse_args()

#if args.pos_arg > 10:
#        parser.error("pos_arg cannot be larger than 10")

print("Argument values:")
print(args.GFF)
print(args.BED)
