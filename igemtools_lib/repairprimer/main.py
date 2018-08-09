from .tools import *
from sys import argv

def main(*argv):
    if len(argv) <= 1:
        return
    filename = argv[1]
    outfile = argv[2]

    final_bases = []
    seqs = read_fa(filename)
    compatible = False

    while not compatible:
        working_bases = [
            generate_bases()
            for i in range(len(seqs) - 1)
        ]

        if is_compatible(seqs, working_bases):
            compatible = True
            final_bases = working_bases

    write_fa(outfile, final_bases)

    print("Done!")

if __name__ == "__main__":
    main(*argv)
