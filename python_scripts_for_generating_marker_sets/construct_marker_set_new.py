"""
Script used to construct the marker gene set for each taxon.

input: a taxon name, that could map to the names in the tree.
output: a dict that is ultimately saved in pickle format, to be written into the marker set tsv file.
"""

import pickle
import argparse
import os
import json
import sys
from calculate_marker_set_utils_new import *
import ete3
from ete3 import Tree
import pandas as pd
import multiprocessing as mp

# from colocalize_marker_sets_new import *


def main():
    # Still have time to modify. To use the same structure as VirRep.
    parser = argparse.ArgumentParser(
        description="Construct a marker set for the taxon given (may also request taxonomic level given)"
    )
    parser.add_argument(
        "--taxon",
        type=str,
        required=False,
        default="/data/A2K_v2.0/config_files/taxon.txt",
        help="taxon to calculate marker set, only incorporated when node-wise argument is False.",
    )
    parser.add_argument(
        "--Stramenopiles",
        help="Caculate the marker set for Stramenopiles, as there are many singleton groups in the tree.",
        required=False,
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=False,
        help="output path of pickle file storing taxon marker set dict, \
            only incorporated when node-wise argument is False.",
    )
    parser.add_argument(
        "--input-tree",
        type=str,
        required=True,
        help="The input tree file that stores the topology of the phylogenetic tree, must be imputed.",
    )
    parser.add_argument(
        "--node-wise",
        "-n",
        action="store_true",
        default=False,
        required=False,
        help="Determine whether do taxon-key guided marker set calculations or iterate over all nodes \
        and store marker sets. Once this argument is set, taxon key, level and output will no longer be needed. \
            Must be imputed. [True|False]",
    )
    parser.add_argument(
        "--overlap-map",
        type=str,
        required=False,
        default="/data/2025Proj_Wint_Spr/merged_pfam_protein/overlap_map.pkl",
        help="A file that contains the overlap map of marker genes, pickle format.",
    )
    parser.add_argument(
        "--protein2pfam",
        type=str,
        required=False,
        default="/data/pfam2protein_0.6.tsv",
        help="A file that contains the mapping from protein to pfam IDs, pickle format.",
    )
    parser.add_argument(
        "--output-tree",
        type=str,
        required=False,
        default=None,
        help="The output directory to store the node-wise marker sets. \
            Only incorporated when node-wise argument is True.",
    )
    parser.add_argument(
        "--single-copy-threshold",
        type=float,
        required=False,
        default=0.97,
        help="The threshold for single-copy marker genes.",
    )

    parser.add_argument(
        "--complex-relationship",
        type=str,
        help="A Json file that contains endosymbiotic \
                        relationships among groups' plastids.\
                        e.g.: {'green_algae': ['Euglenozoa', 'Chlorarachniophyta', 'Dinophyceae'], \
                        'Bacillariophyta': ['Dinophyceae'], 'Dictyochophyceae': ['Centrohelida']}.",
        default="/home/student.unimelb.edu.au/yuhtong/andy/data/A2K_v2.0/config_files/endosymbiosis_dict.json",
    )

    parser.add_argument(
        "--hmmsearch-out-dir",
        type=str,
        default="/home/student.unimelb.edu.au/yuhtong/andy/data/A2K_v2.0/marker_gene_db_from_NCBI/hmmsearch_out_dir",
        help="Directory to store the hmmsearch output files. Contain pfam and protein searchouts.",
    )
    parser.add_argument(
        "--merge-ms-methods",
        type=str,
        choices=["universal", "node-wise"],
        default="universal",
        help="Method to merge marker sets. 'universal' for universally conserved pfam-protein pairs, 'node-wise' for node conserved pfam-protein pairs.",
    )
    arguments = parser.parse_args()

    # assign arguments to variables and preprocess them.
    HMM_PFAM_TABLE = os.path.join(
        arguments.hmmsearch_out_dir, "merged_hmmsearch_output_pfam.txt"
    )
    HMM_PROTEIN_TABLE = os.path.join(
        arguments.hmmsearch_out_dir, "merged_hmmsearch_output_protein.fixed.changed.txt"
    )

    HMM_FEAT = (HMM_PFAM_TABLE, HMM_PROTEIN_TABLE)

    tree = Tree(arguments.input_tree, format=1)  # Load the tree in Newick format.

    pruned_tree = prune_cyan(input_tree=tree)

    tree = preprocess_branches(input_tree=pruned_tree)

    with open(arguments.overlap_map, "rb") as f:
        overlap_map = pickle.load(f)

    if arguments.node_wise:
        if not arguments.output_tree:
            print(
                "Error: --output-tree argument must be provided when node-wise is True."
            )
            sys.exit(1)
        if not arguments.hmmsearch_out_dir:
            print(
                "Error: --hmmsearch-out-dir argument must be provided when node-wise is True."
            )
            sys.exit(1)

        # Now tree is without bacterial leaves, we can calculate marker sets for plastids only.
        endosymbiosis_dict = (
            json.load(open(arguments.complex_relationship, "r"))
            if arguments.complex_relationship
            else None
        )
        for key, value in endosymbiosis_dict.items():
            # Convert the list of strings to a set for faster lookup.
            endosymbiosis_dict[key] = [group.strip() for group in value.split(",")]

        Calculate_nodes_marker_set(
            HMM_FEAT,
            tree=tree,
            merge_method=arguments.merge_ms_methods,
            outfile=arguments.output_tree,
            single_copy_threshold=arguments.single_copy_threshold,
            overlap_map=overlap_map,
            protein2pfam=arguments.protein2pfam,
        )
        # Final output: a phylogenetic tree (unvisualizable), with each node annotated with
        # single-copy marker genes (uncolocalized).
        if arguments.Stramenopiles:
            print("Calculating marker sets for Stramenopiles...")
            Calculate_stramenopiles_marker_set(
                tree=arguments.output_tree,
                output=arguments.output_dir,
                hmm_pfam_proteins=HMM_FEAT,
            )

    else:
        # If not node-wise, we need to read the taxon file and level.
        if not arguments.taxon or not arguments.output_dir:
            print(
                "Error: --taxon, --level and --output arguments must be provided when node-wise is False."
            )
            sys.exit(1)

        if not os.path.exists(arguments.taxon):
            print(f"Error: Taxon file {arguments.taxon} does not exist.")
            sys.exit(1)

        # Read the taxon file and level.
        taxon = arguments.taxon  # the taxon must be one valid taxon from the tree.
        output = arguments.output_dir
        # Check if the taxon is valid.
        if not os.path.exists(arguments.input_tree):
            print(f"Error: Input tree file {arguments.input_tree} does not exist.")
            sys.exit(1)

        # This doesn't need to be parallelized, as it is only for one taxon.
        Calculate_lineage_marker_set(
            taxon=taxon,
            output=output,
            tree=tree,
            hmm_pfam_proteins=HMM_FEAT,
            endosymbiosis_dict=(
                json.load(open(arguments.complex_relationship, "r"))
                if arguments.complex_relationship
                else None
            ),
            single_copy_threshold=arguments.single_copy_threshold,
            overlap_map=overlap_map,
            protein2pfam=arguments.protein2pfam,
            merge_method=arguments.merge_ms_methods,
        )

    print("Marker set calculation completed successfully.")

    return


if __name__ == "__main__":
    main()
