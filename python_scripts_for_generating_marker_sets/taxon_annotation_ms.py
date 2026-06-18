import argparse
from ete3 import *
import pickle as pkl
from os import listdir
import os
import sys
from pathlib import Path
import ete3
from ete3.coretype.tree import TreeError
import pandas as pd
from calculate_marker_set_utils import (
    get_all_descendants,
    find_conserved_markers,
    unify_species_keys,
    iterative_match_species_to_key,
)
from colocalize_marker_sets import *


def colocalize_universal_marker_sets(
    input_tree: ete3.coretype.tree,
    universal_marker_set: list,
    locus_to_species: dict,
    genbank_files: list,
) -> dict:
    # This function will colocalize the universal marker sets for each taxon.
    return


def store_universal_marker_sets(marker_set: dict, output: str) -> None:
    return


def Pick_Node_marker_sets(
    input_tree: ete3.coretype.tree,
    taxon_key: str,
    endosymbiosis_map: dict,
    full_taxon_hierarchy: str,
    taxon_rank: int,
    discarded_species: list,
    locus_to_species: dict,
    genbanks_dir: Path,
) -> dict:
    # For 457 nodes, iterate to find the node that exactly shows the feature desired for each taxon key.
    """
    Here is the thing: we want to find which node shows extremely high-level
    of monophyly of the taxon chosen (looked up tree before), with a taxon level written in taxon_key line.

    Simply, for each node, check the all_leaves in relation to the node itself and to the whole tree, to
    identify the monophyly of node descendants. (Skip the root for this).

    Then, if the nodes' descendants show non-monophyletic group in relation to one tree, and did not show the
    monophyly for one node, the taxon group is removed in selection. If the taxon key shows the highest monophyly
    degree in one node, then we gonna choose the node.

    If we've chosen the node and its corresponding marker set, we now gonna construct a dictionary to
    save all essential elements for the tsv. Lets begin.

    """
    # Core is the iteration over nodes.
    MARKER = "marker_set.pkl"
    NODE = "node.pkl"

    tsv_file_contents = dict()

    Target_Taxon_members = [
        leaf.name for leaf in input_tree.get_leaves() if taxon_key in leaf.name
    ]
    if (
        taxon_key == "Rhodophyta"
        or taxon_key == "Streptophyta"
        or taxon_key == "Florideophyceae"
    ):  # They all got misclassified species.
        Num_Target_Taxon = len(Target_Taxon_members) - 1
    else:
        Num_Target_Taxon = len(Target_Taxon_members)
    print("taxon is {}".format(taxon_key))

    # Important: postorder checks smaller clades first and go with larger clades later.
    for node in input_tree.traverse("postorder"):
        if node.is_leaf():
            continue
        elif any("marker_set" in elem for elem in node.__dict__.items()):
            marker_set = node.marker_set
            if node.is_root():
                continue
            else:
                # We already know that there is one genome mapping to Grateloupa that is misclassified.
                picked_taxon = [
                    leaf.name for leaf in node.get_leaves() if taxon_key in leaf.name
                ]
                Num_Node_Genomes = len(picked_taxon)
                if (
                    Num_Node_Genomes == Num_Target_Taxon
                ):  # Locate branches with the target taxon included.
                    try:
                        monophyly_tuple = node.check_monophyly(
                            picked_taxon, target_attr="name"
                        )
                        print(monophyly_tuple, "\n", len(node.get_leaves()))
                        left_genomes = monophyly_tuple[-1]
                        if (
                            len(left_genomes) == 0
                        ):  # Most taxa in taxon list incl. Rhodophyta.
                            # print("The node is monophyletic.")
                            print(marker_set.split("_"))
                            tsv_file_contents[taxon_key] = dict()
                            tsv_file_contents[taxon_key]["rank"] = taxon_rank
                            tsv_file_contents[taxon_key]["name"] = taxon_key
                            tsv_file_contents[taxon_key][
                                "lineage"
                            ] = full_taxon_hierarchy
                            tsv_file_contents[taxon_key]["descendant_number"] = len(
                                picked_taxon
                            )
                            tsv_file_contents[taxon_key]["marker_set_number"] = len(
                                marker_set.split("@")
                            )
                            marker_set = marker_set.split("@")
                            final_colocalized_marker_sets = Colocalize_marker_sets(
                                picked_leaves=picked_taxon,
                                marker_set=marker_set,
                                locus_to_species=locus_to_species,
                                genbank_files=genbanks_dir,
                            )
                            tsv_file_contents[taxon_key]["set_number"] = (
                                final_colocalized_marker_sets.count("set")
                            )
                            tsv_file_contents[taxon_key][
                                "marker_set"
                            ] = final_colocalized_marker_sets
                            break
                        else:
                            waiver_count, is_complete_waiver = endosymbiotic_wavier(
                                taxon_key, left_genomes, endosymbiosis_map
                            )
                            count_misclassified, contains_discard = (
                                count_misclassified_species(
                                    left_genomes, discarded_species
                                )
                            )
                            sum_irregular_descendants = (
                                waiver_count + count_misclassified
                            )
                            if is_complete_waiver and waiver_count == len(left_genomes):
                                tsv_file_contents[taxon_key] = dict()
                                tsv_file_contents[taxon_key]["rank"] = taxon_rank
                                tsv_file_contents[taxon_key]["name"] = taxon_key
                                tsv_file_contents[taxon_key][
                                    "lineage"
                                ] = full_taxon_hierarchy
                                tsv_file_contents[taxon_key]["descendant_number"] = len(
                                    picked_taxon
                                )
                                tsv_file_contents[taxon_key]["marker_set_number"] = len(
                                    marker_set.split("@")
                                )
                                marker_set = marker_set.split("@")
                                final_colocalized_marker_sets = Colocalize_marker_sets(
                                    picked_leaves=picked_taxon,
                                    marker_set=marker_set,
                                    locus_to_species=locus_to_species,
                                    genbank_files=genbanks_dir,
                                )
                                tsv_file_contents[taxon_key]["set_number"] = (
                                    final_colocalized_marker_sets.count("set")
                                )
                                tsv_file_contents[taxon_key][
                                    "marker_set"
                                ] = final_colocalized_marker_sets
                                break

                            elif (
                                not is_complete_waiver and contains_discard
                            ) and sum_irregular_descendants == len(
                                left_genomes
                            ):  # Chlorophyta.
                                # Construct the dictionary for tsv file.
                                tsv_file_contents[taxon_key] = dict()
                                tsv_file_contents[taxon_key]["rank"] = taxon_rank
                                tsv_file_contents[taxon_key]["name"] = taxon_key
                                tsv_file_contents[taxon_key][
                                    "lineage"
                                ] = full_taxon_hierarchy
                                tsv_file_contents[taxon_key]["descendant_number"] = len(
                                    picked_taxon
                                )
                                tsv_file_contents[taxon_key]["marker_set_number"] = len(
                                    marker_set.split("@")
                                )
                                marker_set = marker_set.split("@")
                                final_colocalized_marker_sets = Colocalize_marker_sets(
                                    picked_leaves=picked_taxon,
                                    marker_set=marker_set,
                                    locus_to_species=locus_to_species,
                                    genbank_files=genbanks_dir,
                                )
                                tsv_file_contents[taxon_key]["set_number"] = (
                                    final_colocalized_marker_sets.count("set")
                                )
                                tsv_file_contents[taxon_key][
                                    "marker_set"
                                ] = final_colocalized_marker_sets
                                break
                            elif contains_discard and count_misclassified == len(
                                left_genomes
                            ):  # Sar group.
                                # Construct the dictionary for tsv file.
                                tsv_file_contents[taxon_key] = dict()
                                tsv_file_contents[taxon_key]["rank"] = taxon_rank
                                tsv_file_contents[taxon_key]["name"] = taxon_key
                                tsv_file_contents[taxon_key][
                                    "lineage"
                                ] = full_taxon_hierarchy
                                tsv_file_contents[taxon_key]["descendant_number"] = len(
                                    picked_taxon
                                )
                                tsv_file_contents[taxon_key]["marker_set_number"] = len(
                                    marker_set.split("@")
                                )
                                marker_set = marker_set.split("@")
                                final_colocalized_marker_sets = Colocalize_marker_sets(
                                    picked_leaves=picked_taxon,
                                    marker_set=marker_set,
                                    locus_to_species=locus_to_species,
                                    genbank_files=genbanks_dir,
                                )
                                tsv_file_contents[taxon_key]["set_number"] = (
                                    final_colocalized_marker_sets.count("set")
                                )
                                tsv_file_contents[taxon_key][
                                    "marker_set"
                                ] = final_colocalized_marker_sets
                                break

                            else:
                                continue

                    except TreeError:
                        continue
                else:
                    continue

        else:
            continue

    # Colocalize marker sets for each taxon.
    # tsv_file_contents is dictionary that saves each line to write. One key is marker set, based on the selected genes,
    # we can isolate the taxa inside the tree, mapping it to the list_of_genbanks_plastome, scan for gene positions.

    print(tsv_file_contents)
    return tsv_file_contents


def endosymbiotic_wavier(
    taxon_key: str, left_descendants: list, endosymbiosis_map: dict
):
    # Compare the left_descendants, whether these leaves are in waiver map.
    if taxon_key not in endosymbiosis_map.keys():
        return False
    endosymbiotic_group = endosymbiosis_map[taxon_key]
    count = 0
    for descendant in left_descendants:
        if any(waiver_taxon in descendant.name for waiver_taxon in endosymbiotic_group):
            count += 1
    is_waiver = True if count == len(left_descendants) else False
    return count, is_waiver


def count_misclassified_species(left_descendants: list, discarded_species: list) -> int:
    # breakpoint()
    count_discarded = len(
        [
            leaf.name
            for leaf in left_descendants
            if any(misclassified in leaf.name for misclassified in discarded_species)
        ]
    )
    return count_discarded, bool(count_discarded)


# Comment: using root tree cannot work out monophyly, as for groups such as Chlorophyta,
# the mosaic grouping makes this quite hard. The only solution is to look at the polyphyly's residual groups.
# If they are in consensus with what evolutionary studies and the tree conveys, then still make Chlorophyta computable.
# Remember to make sure that misannotated or misclassified taxa should be firstly discarded.


# Convenience is that ete3 check monophyly also returns extra leaves that are not in the same group. Alright, we might need a dict to map to these info.
def write_marker_sets_tsv(marker_sets: dict, output: str) -> None:
    # This function will write all marker_sets to a file with fixed name: taxon_marker_sets_lineage_sorted.tsv.
    FILE = "taxon_marker_sets_lineage_sorted.tsv"
    marker_set_df = pd.DataFrame()
    list_info = []
    for marker_set in marker_sets:
        for _, marker_set_dict in marker_set.items():
            list_info.append(marker_set_dict)

    marker_set_df = pd.DataFrame(list_info)
    marker_set_df.to_csv(output, sep="\t", index=False, header=False)
    return


def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Give a list of taxa and iterate over nodes to pick the right one.",
    )

    parser.add_argument("--taxon-list", "-t", type=Path, required=True, default=None)

    parser.add_argument(
        "--universal-marker-set",
        "-u",
        type=bool,
        required=False,
        default=False,
        help="Whether to use the universal marker set or not.",
    )

    parser.add_argument(
        "--universal-marker-set-output",
        type=str,
        required=False,
        help="Output path of the tsv file storing the universal marker sets.",
    )

    parser.add_argument("--input-tree", "-i", type=Path, required=True, default=None)

    parser.add_argument(
        "--endosymbiosis-map",
        "-e",
        type=Path,
        required=False,
        help="A pickle dict object storing endosymbiotic groups within bigger algal groups.",
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=False,
        help="Output path of the tsv file storing all marker sets.",
    )

    parser.add_argument(
        "--discarded",
        "-d",
        type=Path,
        required=False,
        help="Pre-identified misclassified algal species.",
        default=Path("/data/discarded.txt"),
    )

    parser.add_argument(
        "--species-to-locus",
        "-s",
        type=Path,
        required=False,
        default=Path("/data/locus_to_species.pkl"),
        help="Directory storing locus tag mapping to species information.",
    )

    parser.add_argument(
        "--genbank-dir",
        "-g",
        type=Path,
        required=False,
        default=Path("/data/list_of_genbanks_of_algal_plastomes"),
        help="Directory storing genbank files",
    )

    arguments = parser.parse_args()
    Whole_tree = Tree(str(arguments.input_tree), format=1)
    taxon_df = pd.read_csv(str(arguments.taxon_list))

    if arguments.universal_marker_set:
        universal_marker_sets = colocalize_universal_marker_sets(
            input_tree=Whole_tree,
            universal_marker_set=arguments.universal_marker_set,
            locus_to_species=arguments.species_to_locus,
            genbanks_dir=arguments.genbank_dir,
        )
        store_universal_marker_sets(
            marker_set=universal_marker_sets,
            output=arguments.universal_marker_set_output,
        )
        return
        # just calculate this using all genomes.
    with open(str(arguments.endosymbiosis_map), "rb") as em:
        endosymbiosis_map_dict = pkl.load(em)
    with open(str(arguments.discarded), "r") as d:
        discards = d.readlines()

    # Load the species_to_locus dict.
    with open(str(arguments.species_to_locus), "rb") as stl:
        species_to_locus = pkl.load(stl)

    # List the genbank files.
    genbank_dir = [
        f"{str(arguments.genbank_dir)}/{genbank_entry}"
        for genbank_entry in listdir(str(arguments.genbank_dir))
    ]
    discarded_species = [line.strip() for line in discards]
    Collection_of_marker_sets_writable = []
    for taxon in taxon_df.iterrows():
        taxon_key = taxon[1][0].split(";")[-1]
        full_taxon_hierarchy = str(taxon[1][0])
        rank = taxon[1][2]
        Collection_of_marker_sets_writable.append(
            Pick_Node_marker_sets(
                input_tree=Whole_tree,
                taxon_key=taxon_key,
                endosymbiosis_map=endosymbiosis_map_dict,
                full_taxon_hierarchy=full_taxon_hierarchy,
                taxon_rank=rank,
                discarded_species=discarded_species,
                locus_to_species=species_to_locus,
                genbanks_dir=genbank_dir,
            )
        )
    write_marker_sets_tsv(Collection_of_marker_sets_writable, arguments.output)
    return


if __name__ == "__main__":
    main()
