"""Build a solvated Trp-cage AMBER pair for parm2pq."""

import sys

from PQAnalysis.io.protein import prepare_protein

if __name__ == "__main__":
    print(
        prepare_protein(sys.argv[1] if len(sys.argv) > 1 else "1L2Y.pdb", ".")
    )
