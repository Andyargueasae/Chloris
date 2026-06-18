from heapq import merge
from ete3 import Tree
import ete3.coretype
from numpy import argmin
import pandas as pd
import ete3.coretype.tree
from collections import Counter
import re
import multiprocessing as mp
import os
import tqdm
import networkx as nx
import sys
import pickle
from colocalize_marker_set_util_new import *

PFAM_ACCESSION = "query_accession"
PROTEIN_ACCESSION = "query_name"


# We now moved our gene content from species_genome_dict to tabular hmmsearch output.
# The genes are predicted by prodigal and then annotated by hmmsearch.
class TreeNode(ete3.coretype.tree.TreeNode):
    def __init__(self, node):
        super().__init__()
        self.node = node

    def node_marker_sets(
        self,
        hmm_gene_table: tuple,
        single_copy_threshold=0.97,
        overlap_map=None,
        protein2pfam=None,
        merge_method="universal",
    ):
        descendants = self.node.get_leaves()
        genbank_ids = [leaf.name.split("@")[-1] for leaf in descendants]
        # using the hmm_gene_tables, find marker genes single-copy in these descendants.
        single_copy_genes = find_single_copy_genes(
            genbank_ids, hmm_gene_table, single_copy_threshold
        )
        # After finding single-copy genes, output the one from full-length proteins and process them to
        # merge them, then output the merged version.
        assert (
            type(single_copy_genes) == tuple
        ), "single_copy_genes should be a tuple of (pfam_markers, protein_markers)"
        pfam_markers, protein_markers = single_copy_genes
        merged_pfam_protein_set = merge_pfam_proteins(
            pfam_markers,
            protein_markers,
            genbank_ids,
            hmm_gene_table,
            overlap_map,
            protein2pfam,
        )

        return protein_markers, pfam_markers, merged_pfam_protein_set

    def taxon_marker_sets(
        self,
        taxon_name,
        picked_taxon,
        endosymbiosis_dict=None,
        merge_method="universal",
        hmm_gene_table: tuple = None,
        single_copy_threshold=0.97,
        overlap_map=None,
        protein2pfam=None,
    ):
        """
        Calculate marker sets for a taxon, which is defined as a monophyletic group in the tree.
        The taxon is defined by the chosen_node, which is a TreeNode object.
        """
        # first, check in your node: whether all leaves are of same taxon (monophyly).
        if taxon_name == "all":
            # this is the root node, skip the monophyly check.
            print("Calculating universal marker sets for all plastids in the tree.")
        else:
            picked_taxon = [leaf.name for leaf in picked_taxon]
            # print(ta)
            monophyly_tuple = self.node.check_monophyly(
                picked_taxon, target_attr="name"
            )
            left_genomes = monophyly_tuple[-1]
            # if the node is monophyletic, left genomes is 0.
            # else there are some endosymbiotic lineages, prune them and do calculations.
            if len(left_genomes) != 0:
                keep_genomes = [
                    leaf.name
                    for leaf in self.node.get_leaves()
                    if leaf not in left_genomes
                ]
                self.node.prune(keep_genomes, preserve_branch_length=True)
                assert (
                    len(self.node.check_monophyly(picked_taxon, target_attr="name")[-1])
                    == 0
                ), "The node is not monophyletic after pruning, please check your input tree."
        descendants = self.node.get_leaves()
        genbank_ids = [leaf.name.split("@")[-1] for leaf in descendants]
        # using the hmm_gene_tables, find marker genes single-copy in these descendants.
        single_copy_genes = find_single_copy_genes(
            genbank_ids, hmm_gene_table, single_copy_threshold
        )
        # After finding single-copy genes, output the one from full-length proteins and process them to
        # merge them, then output the merged version.
        assert (
            type(single_copy_genes) == tuple
        ), "single_copy_genes should be a tuple of (pfam_markers, protein_markers)"
        pfam_markers, protein_markers = single_copy_genes
        if merge_method == "universal":
            # merge the pfam and protein markers by universal conserved pairs.
            merged_pfam_protein_set = merge_pfam_proteins(
                pfam_markers,
                protein_markers,
                genbank_ids,
                hmm_gene_table,
                overlap_map,
                protein2pfam=protein2pfam,
            )
        else:
            merged_pfam_protein_set = merge_pfam_proteins(
                pfam_markers,
                protein_markers,
                genbank_ids,
                hmm_gene_table,
                overlap_map,
                protein2pfam=None,
            )
        return protein_markers, pfam_markers, merged_pfam_protein_set


def merge_pfam_proteins(
    pfam_markers: set,
    protein_markers: set,
    genbank_ids: list,
    hmm_gene_table: tuple,
    overlap_map: dict,
    protein2pfam=None,
):
    # iterate over each protein marker (as they generally covers more of genome),
    # find the corresponding pfam markers that overlaps with the protein markers.
    pfam_table, protein_table = hmm_gene_table
    pfam_df = pd.read_csv(pfam_table, sep="\t", header=0)
    protein_df = pd.read_csv(protein_table, sep="\t", header=0)

    # This is where the problem arises.

    pfam_df["genome"] = pfam_df["genome"].apply(
        lambda x: x.split(".")[0]
    )  # Extract genbank id.
    protein_df["genome"] = protein_df["genome"].apply(
        lambda x: x.split(".")[0]
    )  # Extract genbank id.

    pfam_filtered = pfam_df[
        ["genome", "target_name", "query_accession", "start_pos", "end_pos"]
    ]
    protein_filtered = protein_df[
        ["genome", "target_name", "query_name", "start_pos", "end_pos"]
    ]

    # Filter the pfam_df and protein_df by genbank_ids.
    pfam_filtered = pfam_filtered[pfam_filtered["genome"].isin(genbank_ids)]
    protein_filtered = protein_filtered[protein_filtered["genome"].isin(genbank_ids)]
    if overlap_map is None:
        print(
            "Overlap map is not provided, recommended to run merge marker sets scripts in 2025Proj_Wint_Spr."
        )
        sys.exit(1)
    overlap_map_filtered = {
        k.split(".")[0]: v
        for k, v in overlap_map.items()
        if k.split(".")[0] in genbank_ids
    }
    threshold = 0.60  # try a lower number: 0.8/0.7?

    essential_pair, _ = identify_conserved_pairs(overlap_map_filtered, threshold)
    # breakpoint()
    # think of adding this old approach to the database.
    if protein2pfam is not None:
        # merge the pfam and protein markers by protein2pfam.
        protein2pfam = pd.read_csv(
            protein2pfam, sep="\t", header=None, names=["protein", "pfam"]
        )

        merged_markers = merge_pfam_protein_ms_protein2pfam(
            pfam_markers, protein_markers, protein2pfam
        )
    else:
        merged_markers = merge_pfam_protein_ms_pair_count(
            pfam_markers, protein_markers, essential_pair
        )
    # merged_markers = merge_pfam_protein_ms_protein2pfam(pfam_markers, protein_markers, protein2pfam)

    return merged_markers


def compute_overlap_map(pfam_filtered, protein_filtered):
    Overlap_map = {}
    # use proteins as query, as proteins are longer than domains, which may overlap with multiple domains.
    # use target_name as the protein id.
    for genome, df in protein_filtered.groupby("genome"):
        if genome not in Overlap_map:
            Overlap_map[genome] = {}
        for _, protein_row in df.iterrows():
            start_pos = protein_row["start_pos"]
            end_pos = protein_row["end_pos"]
            # we need the hit's query name.
            protein_name = protein_row[
                "query_name"
            ]  # this will be the key in the overlap map.

            # find overlapping domains
            overlapping_domains = pfam_filtered[
                (pfam_filtered["genome"] == genome)
                & (pfam_filtered["start_pos"] < end_pos)
                & (pfam_filtered["end_pos"] > start_pos)
            ]
            domain_list = overlapping_domains["query_accession"].tolist()
            if protein_name not in Overlap_map:
                Overlap_map[genome][protein_name] = list(set(domain_list))
            else:
                Overlap_map[genome][protein_name].extend(domain_list)
                Overlap_map[genome][protein_name] = list(set(Overlap_map[protein_name]))
    """        
    finally, this dictionary will be of the form:
    {
        'genome1': {
                'protein_name1': ['domain1', 'domain2'],
                'protein_name2': ['domain3']
        },
        'genome2': {
                'protein_name3': ['domain4', 'domain5']
        }
    }
    """
    return Overlap_map


def identify_conserved_pairs(overlap_map, threshold=0.9):
    """
    Aim: identify conserved pairs of domains and proteins based on the overlap map.
    Save them into a list of tuples where each tuple contains a protein and a domain.
    """

    pair_count = dict()

    count_genome = len(overlap_map)  # total number of genomes

    for genome, values in overlap_map.items():
        # values is a dictionary storing proteins and their overlapping domains
        for protein, domains in values.items():
            for domain in domains:
                pair = (protein, domain)
                if pair not in pair_count:
                    pair_count[pair] = 0
                pair_count[pair] += 1
    # Now we have a count of how many genomes each pair appears in.
    essential_pairs = []
    for pair, count in pair_count.items():
        if count / count_genome >= threshold:
            essential_pairs.append(pair)
    # breakpoint()
    return essential_pairs, pair_count


def merge_pfam_protein_ms_pair_count(domain_markers, protein_markers, essential_pairs):
    # merge proteins by finding overlapping domains in the pair_count.

    merged_markers = set(protein_markers).union(set(domain_markers))
    # iterate over the pair_count, if the protein marker is in the pair_count, we will keep the domain marker,
    # while discarding the protein marker.
    for protein in protein_markers:
        involved_pairs = [pair for pair in essential_pairs if pair[0] == protein]
        if involved_pairs:
            count = 0
            for pair in involved_pairs:
                domain = pair[1]
                if domain in domain_markers:
                    count += 1
                    # remove the protein marker.
            if count == len(involved_pairs):
                # if all domains are in domain_markers, we will remove the protein marker.
                merged_markers.discard(protein)
                merged_markers.update([pair[1] for pair in involved_pairs])
            else:
                # if not all domains are in domain_markers, we will keep the protein marker.
                # but at the same time, if there are overlapping domains, discard them.
                for pair in involved_pairs:
                    domain = pair[1]
                    if domain in merged_markers:
                        merged_markers.discard(domain)
                        # we will keep the protein marker while removing the domains.
        else:
            # if the protein marker is not in the pair_count, we will keep it.
            continue

    return merged_markers


def merge_pfam_protein_ms_protein2pfam(domain_markers, protein_markers, protein2pfam):
    """
    Aim: merge the separate pfam marker set and protein marker set into one, but list those non-overlapping markers.
    for those overlapping markers, we will choose the one with the best prevalence and single-copy status.
    """
    # use protein markers as the key, iterate over the overlap map, and find the matching pfam markers.
    # We change the strategy: merge all domains and proteins first,
    # remove the overlapping pairs by retaining only domains.

    # protein2pfam dictionary.
    # 刻舟求剑，不可取。
    protein2pfam_dict = {}
    for index, row in protein2pfam.iterrows():
        protein = row["protein"]
        pfam = row["pfam"]
        if not protein2pfam_dict.get(protein):
            protein2pfam_dict[protein] = [pfam.replace("pfam", "PF")]
        else:
            protein2pfam_dict[protein].append(pfam.replace("pfam", "PF"))
        if not protein2pfam_dict.get(pfam):
            protein2pfam_dict[pfam] = [protein]
        else:
            protein2pfam_dict[pfam].append(protein)

    merged_markers2 = set(protein_markers).union(set(domain_markers))
    for protein in protein_markers:
        # see if there are any pfams for proteins in protein markers to map to.
        if protein2pfam_dict.get(protein):
            # if there are, we will merge them into the merged_markers2.
            count = 0
            pfams = protein2pfam_dict[protein]
            for pfam in pfams:
                if pfam in domain_markers:
                    # if the pfam is in domain markers, we will keep it.
                    count += 1
                    # remove the protein.
                else:
                    # if the pfam is not in domain markers, we will add it to the merged_markers2.
                    continue
            if count == len(pfams):
                # if all pfams are in domain markers, we will remove the protein from merged_markers2.
                merged_markers2.discard(protein)
                merged_markers2.update(pfams)
            else:
                for pfam in pfams:
                    if pfam in merged_markers2:
                        merged_markers2.discard(pfam)
                        # we will keep the protein marker while removing the domains.
        else:
            continue

    return merged_markers2


def prune_cyan(input_tree: ete3.coretype.tree.Tree):
    # Prune the cyanobacterial leaves from the tree.
    removed_leaves = []
    for leaf in input_tree.get_leaves():
        if "Bacteria" in leaf.name:
            removed_leaves.append(leaf.name)
        if "Centroplasthelida" in leaf.name:
            removed_leaves.append(leaf.name)

    keep_leaves = [
        leaf for leaf in input_tree.get_leaves() if leaf.name not in removed_leaves
    ]
    input_tree.prune(keep_leaves, preserve_branch_length=True)
    return input_tree


def preprocess_branches(input_tree):
    """
    Considering that some strings may intervene with taxon key, we should preprocess the branches,
    to make every branch distinguishable in terms of their lineage and species name.
    """
    # for every leaf in the tree, firstly isolate the first 6 taxonomic levels split by "_", do not change them.
    # then, for the rest, the last part split by "_" is the genbank id, keep it, then for the rest, delimit them by space or "@".

    all_genbanks = [leaf.name.split("_")[-1] for leaf in input_tree.get_leaves()]

    for leaf in input_tree.get_leaves():
        leaf_name = leaf.name
        # preserve the first 6 taxonomic levels.
        lineage, species_genbank = (
            leaf_name.split("_", 6)[:6],
            leaf_name.split("_", 6)[-1],
        )
        # The pattern of genbank id is: one or more alphabetical characters followed by multiple digits, alphabets and digits maybe separated by "_".
        genbank_match = re.findall(r"\_([a-zA-Z]{1,2}_\d{5,6})", species_genbank)
        if genbank_match:
            # remove the first underscore from the genbank match.
            genbank = genbank_match[-1].lstrip("_")
            # then species may be the left characters.
            species = species_genbank.replace(genbank, "").strip().split("_")
            # concatenate the strings from species by " ", and genbank by "@".
            leaf.name = "_".join(lineage) + "_" + " ".join(species) + "@" + genbank
        else:
            genbank_match = re.findall(r"\_([a-zA-Z]{1,2}\d{5,6})", species_genbank)
            genbank = (
                genbank_match[-1].lstrip("_") if genbank_match else "unknown_genbank"
            )
            # then species may be the left characters.
            species = species_genbank.replace(genbank, "").strip().split("_")
            # concatenate the strings from species by " ", and genbank by "@".
            leaf.name = "_".join(lineage) + "_" + " ".join(species) + "@" + genbank
    print("Preprocessing branches done.")

    # remove the leaves: AY286123, LC704715, LC716139 and LC716140.
    remove_empty_annotation_gs = []
    for leaf in input_tree.get_leaves():
        if "AY286123" in leaf.name:
            remove_empty_annotation_gs.append(leaf.name)
        elif "LC704715" in leaf.name:
            remove_empty_annotation_gs.append(leaf.name)
        elif "LC716139" in leaf.name:
            remove_empty_annotation_gs.append(leaf.name)
        elif "LC716140" in leaf.name:
            remove_empty_annotation_gs.append(leaf.name)

    keep_leaves = [
        leaf
        for leaf in input_tree.get_leaves()
        if leaf.name not in remove_empty_annotation_gs
    ]
    input_tree.prune(keep_leaves, preserve_branch_length=True)
    return input_tree


def find_single_copy_genes(
    genbank_ids, hmm_gene_table: tuple, single_copy_threshold=0.97
):

    branch_leaf_num = len(genbank_ids)

    pfam_table, protein_table = hmm_gene_table
    pfam_df = pd.read_csv(pfam_table, sep="\t", header=0)
    protein_df = pd.read_csv(protein_table, sep="\t", header=0)

    pfam_df["genome"] = pfam_df["genome"].apply(
        lambda x: x.split(".")[0]
    )  # Extract genbank id.
    protein_df["genome"] = protein_df["genome"].apply(
        lambda x: x.split(".")[0]
    )  # Extract genbank id.

    # Filter the pfam_df and protein_df by genbank_ids.
    pfam_df = pfam_df[pfam_df["genome"].isin(genbank_ids)]
    protein_df = protein_df[protein_df["genome"].isin(genbank_ids)]

    # now, groupby each genome, and find the single-copy genes/pfams in these tables.
    single_copy_pfam = {}
    single_copy_protein = {}

    pfam_markers = set()
    protein_markers = set()
    for genome in genbank_ids:
        genome_pfam = pfam_df[pfam_df["genome"] == genome]
        genome_protein = protein_df[protein_df["genome"] == genome]
        # construct count value table for pfam accession and protein query name.
        pfam_count = genome_pfam[PFAM_ACCESSION].value_counts()
        protein_count = genome_protein[PROTEIN_ACCESSION].value_counts()

        for gene, count in pfam_count.items():
            if count == 1:  # Only want the single-copy pfam.
                if gene not in single_copy_pfam:
                    single_copy_pfam[gene] = []
                single_copy_pfam[gene].append(genome)

        for gene, count in protein_count.items():
            if count == 1:  # Only want the single-copy protein.
                if gene not in single_copy_protein:
                    single_copy_protein[gene] = []
                single_copy_protein[gene].append(genome)

    # Now, we need to filter the single-copy genes/pfams by the threshold.
    for gene, genomes_list in single_copy_pfam.items():
        if len(genomes_list) / branch_leaf_num >= single_copy_threshold:
            # If the number of genomes with this gene is greater than the threshold, we keep it.
            pfam_markers.add(gene)
    for gene, genomes_list in single_copy_protein.items():
        if len(genomes_list) / branch_leaf_num >= single_copy_threshold:
            # If the number of genomes with this gene is greater than the threshold, we keep it.
            protein_markers.add(gene)

    return (pfam_markers, protein_markers)


def count_matching_leaves(node, taxon):
    return len([leaf for leaf in node.get_leaves() if taxon in leaf.name])


def find_node(tree, taxon_name, matching_leaves=None):
    """
    Find the node containing all leaves matching the taxon name.
    Some nodes may contain not only leaves matching, but also endosymbiotic lineages.
    """
    # iterate over all nodes, for those nodes with leaves containing those matching taxon name (and may contain other complex endosymbiotic lineages),
    # we will return the node.
    if taxon_name == "all":
        # return root.
        return tree.get_tree_root()
    if taxon_name == "Stramenopiles":
        # find the node that is the latest common ancestor of all Stramenopiles.
        node = tree.get_common_ancestor(matching_leaves)
        return node if node else None
    candidate_nodes = [
        node
        for node in tree.traverse()
        if count_matching_leaves(node, taxon_name) >= len(matching_leaves)
        and all(leaf.name in node.get_leaf_names() for leaf in matching_leaves)
    ]
    # return the candidate node with minimum number of leaves matching the taxon name.
    return (
        candidate_nodes[argmin([len(node.get_leaves()) for node in candidate_nodes])]
        if candidate_nodes
        else None
    )


def transform_ms2str(marker_set):
    # marker_set is a list of sets of markers, we need to transform it to a string like this: [set(["marker1", "marker2"]), set(["marker3"])]
    if not marker_set:
        return ""
    marker_set_str = "[" + ",".join([f"set({list(ms)})" for ms in marker_set]) + "]"
    return marker_set_str


def Calculate_nodes_marker_set(
    hmm_pfam_proteins,
    tree=None,
    outfile=None,
    single_copy_threshold=0.97,
    overlap_map=None,
    protein2pfam=None,
    merge_method="universal",
):
    # For each of the node, get their descendants, calculate marker sets.
    hmm_pfam = hmm_pfam_proteins[0]
    hmm_protein = hmm_pfam_proteins[1]
    # we need to get the marker sets combining both PFAM and proteins.
    index = 0
    for node in tqdm.tqdm(tree.traverse()):
        if not node.is_leaf():
            Node_obj = TreeNode(node)
            protein_marker_set, pfam_marker_set, merged_marker_set = (
                Node_obj.node_marker_sets(
                    hmm_gene_table=(hmm_pfam, hmm_protein),
                    single_copy_threshold=single_copy_threshold,
                    overlap_map=overlap_map,
                    protein2pfam=protein2pfam,
                    merge_method=merge_method,
                )
            )
            index += 1
            # Save the marker set into the node object, so that the newick-formatted tree has each node calculated with marker sets.
            protein_marker_str = str("@".join(protein_marker_set))
            merged_marker_str = str("@".join(merged_marker_set))
            node.merged_marker_set = merged_marker_str
            node.protein_marker_set = protein_marker_str
            node.pfam_marker_set = str("@".join(pfam_marker_set))
        else:
            continue
    tree.write(
        outfile=outfile,
        features=["protein_marker_set", "merged_marker_set", "pfam_marker_set"],
        format=1,
    )
    return


def Calculate_lineage_marker_set(
    taxon,
    output,
    tree=None,
    hmm_pfam_proteins=None,
    endosymbiosis_dict=None,
    single_copy_threshold=0.97,
    overlap_map=None,
    protein2pfam=None,
    merge_method="universal",
):
    """
    Calculate marker sets for a lineage in the tree.
    The lineage is defined as a monophyletic group in the tree.
    """
    hmm_pfam = hmm_pfam_proteins[0]
    hmm_protein = hmm_pfam_proteins[1]

    taxon_df = pd.read_csv(taxon, sep=",", header=None)

    taxon_df.columns = ["level", "taxon"]

    protein_marker_set_dict = {}
    merged_marker_set_dict = {}
    pfam_marker_set_dict = {}
    phy2acc = {}
    phy2acc_merged = {}

    tree_copy = tree.copy()  # make a copy of the tree, so that we can prune it later.

    if not os.path.exists(output):
        os.makedirs(output)

    for index, row in taxon_df.iterrows():
        taxon_name = row["taxon"]

        # try splitting it: "all" for all organisms, then others are calculated.
        if taxon_name == "all":
            # For "all", we will not prune the tree, but calculate the marker sets for all leaves in the tree.
            print("Calculating marker sets for all plastids in the tree.")
            taxon_name = "all"
            matching_leaves = tree_copy.get_leaves()
            genbank_ids = [leaf.name.split("@")[-1] for leaf in matching_leaves]
        else:
            taxon_name = (
                taxon_name.split(";")[-1] if ";" in taxon_name else taxon_name
            )  # Get the last part of the taxon name.
            matching_leaves = [
                leaf for leaf in tree_copy.get_leaves() if taxon_name in leaf.name
            ]
            genbank_ids = [leaf.name.split("@")[-1] for leaf in matching_leaves]
        # Find the node that contains all the leaves matching the taxon.
        node_of_taxon = find_node(
            tree_copy, taxon_name, matching_leaves=matching_leaves
        )

        # calculate both marker sets.
        protein_marker_set, pfam_marker_set, merged_marker_set = TreeNode(
            node_of_taxon
        ).taxon_marker_sets(
            taxon_name=taxon_name,
            picked_taxon=matching_leaves,
            hmm_gene_table=(hmm_pfam, hmm_protein),
            single_copy_threshold=single_copy_threshold,
            overlap_map=overlap_map,
            protein2pfam=protein2pfam,
            merge_method=merge_method,
        )

        print("taxon", taxon_name, protein_marker_set, "\n", merged_marker_set)
        print(
            "number of proteins:",
            len(protein_marker_set),
            "number of merged markers:",
            len(merged_marker_set),
        )

        # refresh tree_copy.
        tree_copy = tree.copy()  # reset the tree copy for the next iteration.

        colocalized_protein_ms = Colocalize_marker_sets(
            genbank_ids, protein_marker_set, hmm_pfam_proteins
        )
        # colocalized_merged_ms = Colocalize_marker_sets(genbank_ids, merged_marker_set, hmm_pfam_proteins)
        colocalized_pfam_ms = Colocalize_marker_sets(
            genbank_ids, pfam_marker_set, hmm_pfam_proteins
        )
        if taxon_name not in protein_marker_set_dict:
            protein_marker_set_dict[taxon_name] = dict()
            # merged_marker_set_dict[taxon_name] = dict()
            pfam_marker_set_dict[taxon_name] = dict()
            # save for proteins first.
        protein_marker_set_dict[taxon_name]["rank"] = row["level"]
        # merged_marker_set_dict[taxon_name]['rank'] = row['level']
        pfam_marker_set_dict[taxon_name]["rank"] = row["level"]

        if taxon_name == "all":
            # If the taxon is "all", we will not save the level, but just save the number of genomes and markers.
            protein_marker_set_dict[taxon_name]["level"] = "algae"
            # merged_marker_set_dict[taxon_name]['level'] = "algae"
            pfam_marker_set_dict[taxon_name]["level"] = "algae"
        elif taxon_name in [
            "Bacillariophyta",
            "Dictyochophyceae",
            "Eustigmatophyceae",
            "Pelagophyceae",
            "Peridiniales",
            "Phaeophyceae",
            "Raphidophyceae",
            "Synurophyceae",
        ]:
            protein_marker_set_dict[taxon_name][
                "level"
            ] = f"algae;Stramenopiles;{taxon_name}"
        else:
            protein_marker_set_dict[taxon_name]["level"] = f"algae;{taxon_name}"
            # merged_marker_set_dict[taxon_name]['level'] = f"algae;{taxon_name}"
            pfam_marker_set_dict[taxon_name]["level"] = f"algae;{taxon_name}"

        protein_marker_set_dict[taxon_name]["num_genome"] = len(matching_leaves)
        protein_marker_set_dict[taxon_name]["num_markers"] = len(protein_marker_set)
        protein_marker_set_dict[taxon_name]["num_ms"] = len(colocalized_protein_ms)
        protein_marker_set_dict[taxon_name]["markers"] = transform_ms2str(
            colocalized_protein_ms
        )

        # save for merged markers.
        # merged_marker_set_dict[taxon_name]['level'] = f"algae;{taxon_name}"
        # merged_marker_set_dict[taxon_name]["num_genome"] = len(matching_leaves)
        # merged_marker_set_dict[taxon_name]["num_markers"] = len(merged_marker_set)
        # merged_marker_set_dict[taxon_name]["num_ms"] = len(colocalized_merged_ms)
        # merged_marker_set_dict[taxon_name]['markers'] = transform_ms2str(colocalized_merged_ms)

        pfam_marker_set_dict[taxon_name]["rank"] = row["level"]
        pfam_marker_set_dict[taxon_name]["level"] = f"algae;{taxon_name}"
        pfam_marker_set_dict[taxon_name]["num_genome"] = len(matching_leaves)
        pfam_marker_set_dict[taxon_name]["num_markers"] = len(pfam_marker_set)
        pfam_marker_set_dict[taxon_name]["num_ms"] = len(colocalized_pfam_ms)
        pfam_marker_set_dict[taxon_name]["markers"] = transform_ms2str(
            colocalized_pfam_ms
        )

        phy2acc[f"p__{taxon_name}"] = list(protein_marker_set)
        # phy2acc_merged[f"p__{taxon_name}"] = list(merged_marker_set)
        # breakpoint()
    # save in terms of binny's data structure, and individual marker sets in checkm's format.
    # save the marker sets in terms of dictionary: for CompleteBin's Phy2ACC.pkl.
    phy2acc_path = os.path.join(output, "Phy2ACC_prot.pkl")
    # phy2acc_merged_path = os.path.join(output, "Phy2ACC_merged.pkl")
    with open(phy2acc_path, "wb") as f:
        pickle.dump(phy2acc, f)
    # with open(phy2acc_merged_path, 'wb') as f:
    #     pickle.dump(phy2acc_merged, f)

    taxon_marker_set = pd.DataFrame.from_dict(protein_marker_set_dict, orient="index")
    # reset index.
    taxon_marker_set = taxon_marker_set.reset_index()
    taxon_marker_set = taxon_marker_set.sort_values("level")

    new_columns = [
        "rank",
        "index",
        "level",
        "num_genome",
        "num_markers",
        "num_ms",
        "markers",
    ]
    taxon_marker_set = taxon_marker_set[new_columns]
    taxon_marker_set.to_csv(
        os.path.join(output, "taxon_marker_sets_lineage_sorted.tsv"),
        sep="\t",
        header=False,
        index=False,
    )

    # taxon_marker_set_merged = pd.DataFrame.from_dict(merged_marker_set_dict, orient='index')
    # taxon_marker_set_merged = taxon_marker_set_merged.reset_index()
    # taxon_marker_set_merged = taxon_marker_set_merged.sort_values("level")
    # taxon_marker_set_merged.to_csv(os.path.join(output, "taxon_marker_sets_lineage_sorted.merged.tsv"), sep="\t", header=False, index=False)

    taxon_marker_set_pfam = pd.DataFrame.from_dict(pfam_marker_set_dict, orient="index")
    # reset index.
    taxon_marker_set_pfam = taxon_marker_set_pfam.reset_index()
    taxon_marker_set_pfam = taxon_marker_set_pfam.sort_values("level")
    taxon_marker_set_pfam.to_csv(
        os.path.join(output, "taxon_marker_sets_lineage_sorted.pfam.tsv"),
        sep="\t",
        header=False,
        index=False,
    )
    return


def Calculate_stramenopiles_marker_set(tree=None, output=None, hmm_pfam_proteins=None):
    if not os.path.exists(output):
        os.makedirs(output)
    work_tree = Tree(
        tree, format=1
    )  # this is the output tree from the node-wise marker set calculation.
    # nodes = [node for node in work_tree.traverse() if not node.is_leaf()]
    # Because we did not define Stramenopiles in their leaf names, we can only find the node which is the latest common ancestor
    # of Bacillariophyta, Phaeophyceae, Dictyochophyceae, Pelagophyceae, Synurophyceae, Eustigmatophyceae, Raphidophyceae,
    # Xanthophyceae.
    SAR_taxa = [
        "Bacillariophyta",
        "Phaeophyceae",
        "Dictyochophyceae",
        "Pelagophyceae",
        "Synurophyceae",
        "Eustigmatophyceae",
        "Raphidophyceae",
        "Xanthophyceae",
    ]
    # Find the node that contains all the leaves matching the taxon.
    all_matching_leaves = [
        leaf
        for leaf in work_tree.get_leaves()
        if any(taxon in leaf.name for taxon in SAR_taxa)
    ]

    SAR_node = find_node(
        work_tree, "Stramenopiles", matching_leaves=all_matching_leaves
    )
    protein_marker_set = SAR_node.protein_marker_set.split("@")
    pfam_marker_set = SAR_node.pfam_marker_set.split("@")
    # merged_marker_set = SAR_node.merged_marker_set.split("@")

    # Now we have the marker sets, we can colocalize them.
    genbank_ids = [leaf.name.split("@")[-1] for leaf in all_matching_leaves]
    colocalized_protein_ms = Colocalize_marker_sets(
        genbank_ids, protein_marker_set, hmm_pfam_proteins
    )
    # colocalized_merged_ms = Colocalize_marker_sets(genbank_ids, merged_marker_set, hmm_pfam_proteins)
    # breakpoint()
    colocalized_pfam_ms = Colocalize_marker_sets(
        genbank_ids, pfam_marker_set, hmm_pfam_proteins
    )

    SAR_marker_dict = {
        "rank": "phylum",
        "name": "Stramenopiles",
        "level": "algae;Stramenopiles",
        "num_genome": len(all_matching_leaves),
        "num_markers": len(protein_marker_set),
        "num_ms": len(colocalized_protein_ms),
        "markers": transform_ms2str(colocalized_protein_ms),
    }

    # SAR_merged_dict = {
    #     "rank": "phylum",
    #     "name": "Stramenopiles",
    #     "level": "algae;Stramenopiles",
    #     "num_genome": len(all_matching_leaves),
    #     "num_markers": len(merged_marker_set),
    #     "num_ms": len(colocalized_merged_ms),
    #     "markers": transform_ms2str(colocalized_merged_ms)
    # }

    SAR_pfam_dict = {
        "rank": "phylum",
        "name": "Stramenopiles",
        "level": "algae;Stramenopiles",
        "num_genome": len(all_matching_leaves),
        "num_markers": len(pfam_marker_set),
        "num_ms": len(colocalized_pfam_ms),
        "markers": transform_ms2str(colocalized_pfam_ms),
    }

    protein_df = pd.DataFrame.from_dict(
        {"Stramenopiles": SAR_marker_dict}, orient="index"
    )
    # merged_df = pd.DataFrame.from_dict({"Stramenopiles": SAR_merged_dict}, orient='index')
    pfam_df = pd.DataFrame.from_dict({"Stramenopiles": SAR_pfam_dict}, orient="index")

    protein_df = protein_df.reset_index()
    # merged_df = merged_df.reset_index()
    pfam_df = pfam_df.reset_index()
    new_columns = [
        "rank",
        "index",
        "level",
        "num_genome",
        "num_markers",
        "num_ms",
        "markers",
    ]
    protein_df = protein_df[new_columns]
    protein_df.to_csv(
        os.path.join(output, "Stramenopiles_marker_set.tsv"),
        sep="\t",
        header=False,
        index=False,
    )
    # merged_df.to_csv(os.path.join(output, "Stramenopiles_marker_set_merged.tsv"), sep="\t", header=False, index=False)
    pfam_df.to_csv(
        os.path.join(output, "Stramenopiles_marker_set_pfam.tsv"),
        sep="\t",
        header=False,
        index=False,
    )
    return
