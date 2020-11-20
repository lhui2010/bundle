

#all

#parallel run

#train

#liftover
#Require RaGOO
def liftover():
    pass

#function

def main():
    """
    """
    actions = (
        ('compile', 'extract telomere length and ccn'),
        ('traits', 'make HTML page that reports eye and skin color'),
        # Age paper plots
        ('qc', 'plot distributions of basic statistics of a sample'),
        ('correlation', 'plot correlation of age vs. postgenomic features'),
        ('heritability', 'plot composite on heritability estimates'),
        ('regression', 'plot chronological vs. predicted age'),
        ('ccn', 'plot several ccn plots including chr1,chrX,chrY,chrM'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

    print(__file__)
    print(__doc__)
    exit()
    prog_name = "busco_wrapper"
    usage = "run busco on selected GENOME"

    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("GENOME", help="Genome to be evalutated in fasta format")
    parser.add_argument("-t", "--threads", default=64, type=int, help="flanking distance default (1000)")
    args = parser.parse_args()  

    busco(args.GENOME)
#    flanking_distance = args.flanking

if __name__ == "__main__":
    main()
