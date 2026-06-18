"""
Script used to construct the marker gene set for each taxon.

input: a taxon name, that could map to the names in the tree.
output: a dict that is ultimately saved in pickle format, to be written into the marker set tsv file.
"""

import pickle as pkl
import argparse
from calculate_marker_set_utils import *
import os
import sys
from colocalize_marker_sets import *


def main():
    # Still have time to modify. To use the same structure as VirRep.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Construct a marker set for the taxon given (may also request taxonomic level given)",
    )
    parser.add_argument(
        "--taxon",
        metavar="TAXON",
        type=Path,
        required=False,
        help="taxon to calculate marker set, only incorporated when node-wise argument is False.",
    )
    parser.add_argument(
        "--level",
        metavar="LEVEL",
        type=int,
        required=False,
        help="taxon hierarchy level [1-phylum|2-class|3-order|4-family|5-genus|6-species], \
            only incorporated when node-wise argument is False.",
    )
    parser.add_argument(
        "--output",
        metavar="OUTPUT",
        type=str,
        required=False,
        help="output path of pickle file storing taxon marker set dict, \
            only incorporated when node-wise argument is False.",
    )
    parser.add_argument(
        "--input-tree",
        metavar="INPUT_TREE",
        type=str,
        required=True,
        help="The input tree file that stores the topology of the phylogenetic tree, must be imputed.",
    )
    parser.add_argument(
        "--node-wise",
        "-n",
        metavar="NODEWISE",
        type=str,
        default="False",
        required=True,
        help="Determine whether do taxon-key guided marker set calculations or iterate over all nodes \
        and store marker sets. Once this argument is set, taxon key, level and output will no longer be needed. \
            Must be imputed. [True|False]",
    )

    parser.add_argument(
        "--output-tree",
        metavar="OUTPUT_TREE",
        type=str,
        required=False,
        help="The output directory to store the node-wise marker sets. \
            Only incorporated when node-wise argument is True.",
    )

    parser.add_argument(
        "--species-genome-dict",
        metavar="SPECIES_GENOME_DICT",
        type=str,
        required=False,
        help="The dictionary that stores the species' genome content. \
            Only incorporated when node-wise argument is True.",
    )

    arguments = parser.parse_args()
    # If the args say we need all organisms' conserved genes, calculate the gene set using all genomes.
    if arguments.node_wise == "False":
        if str(arguments.taxon) == "all":
            marker_dict = get_universal_markers(arguments.input_tree)  # a dict.
            universal_marker_set = list(marker_dict["marker set"])
            # Have an alternative way to store and colocalize universal marker sets, update it to binny_ChloroScan.
            # save it in pickle format.
            with open(arguments.output, "wb") as f:
                pkl.dump(universal_marker_set, f)
            print("Universal marker sets saved.")
        else:
            marker_dict = dict()
            genomes_incl = get_all_descendants(
                str(arguments.taxon), str(arguments.input_tree)
            )
            marker_set = find_conserved_markers(genomes_incl)  # a set.
            print(
                f"marker sets for {arguments.taxon}: {marker_set}\n length of marker sets: {len(marker_set)}"
            )
            marker_dict["taxon level"] = str(arguments.level)
            marker_dict["taxon"] = str(arguments.taxon)
            marker_dict["marker_set"] = marker_set

        if os.path.exists(arguments.output):
            print("The previous file exists, now remove it.")
            os.remove(arguments.output)
        # Now we work this out.
        with open(arguments.output, "wb") as f:
            pkl.dump(marker_dict, f)
        print("New marker sets saved.")
    elif arguments.node_wise == "True":
        with open(arguments.species_genome_dict, "rb") as f:
            species_genome_dict = pkl.load(f)
        print(
            "Calculate lineage-specific marker sets for each internal node in input tree."
        )
        Calculate_node_marker_set(
            species_genome_dict=species_genome_dict,
            tree_obj=str(arguments.input_tree),
            outfile=str(arguments.output_tree),
        )
    # Here we want to use the tree only and calculate the marker sets for each of the node.
    # To let the outcome represent the node correctly, we want to use the highest rank as the name to summarize.
    else:
        print("The argument for node-wise is not expected, exiting.")
        sys.exit(1)
    return


if __name__ == "__main__":
    main()
