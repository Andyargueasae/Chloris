# Utility for calculating marker set of an algal lineage.
from Bio import SeqIO, Phylo
import pickle as pkl
import re
import ete3
from ete3 import Tree
from collections import Counter
import ete3.coretype
import ete3.coretype.tree
import networkx as nx
import os

# with open("/home/student.unimelb.edu.au/yuhtong/andy/data/species_genome_effective_dict.pkl", "rb") as f:
#     species_genome_dict = pkl.load(f)

ALL = "all"

# These genomes should be excluded from the analysis, due to misannotations and classifications.


class TreeNode(ete3.coretype.tree.TreeNode):
    def __init__(self, node):
        super().__init__()
        self.node = node

    def taxonomy_group_pooling(self):
        # This function of the class pools all leaves' taxonomic level name together.
        # s # Leave objects list.
        node_network = nx.Graph()
        for i in self.node.get_leaves():
            lineage = i.name.split("_")[:-1]
            for j in range(len(lineage)):
                # For the first group, just add it to graph, count only occur once.
                if lineage[j] not in node_network.nodes and j == 0:
                    node_network.add_node(lineage[j])
                elif lineage[j] in node_network.nodes and j == 0:
                    pass
                # But for lower taxonomic ranking, the first presence should be added to graph, and count to 1.
                elif j >= 1 and lineage[j] not in node_network.nodes:
                    node_network.add_node(lineage[j])
                    node_network.add_edge(lineage[j - 1], lineage[j])
                    node_network.nodes[lineage[j]]["count"] = 1
                # Then, for taxonomic group inside one species' lineage and has already been added to graph, simply accumulate counts.
                else:
                    # Note, here we forgot some codes.
                    # If j >= 1 and lineages already in the graph, how we gonna do?
                    # First, the node is in the graph, then we add count by 1.
                    node_network.nodes[lineage[j]]["count"] += 1

        return node_network

    def node_marker_sets(self, species_genome_dict):
        # breakpoint()
        leafnames = [leaf.name for leaf in self.node.get_leaves()]
        sp_gen_dict = species_genome_dict.copy()
        mod_dict, _ = unify_species_keys(sp_gen_dict)
        list_of_genomes = []
        num_genomes = len(leafnames)
        threshold = (
            0.97 * num_genomes
        )  # Marker gene conserved in over 97% of the genomes under a node.
        gene_counts = Counter()
        # Iterate over genomes to save gene contents from species-genome dict.
        for leafname in leafnames:
            gene_content_dict = iterative_match_species_to_key(leafname, mod_dict)
            list_of_genomes.append(gene_content_dict["gene content"]["gene contents"])
        # Now, take a look by counting the gene contents for each genome, count those genes present only once.
        for genome_content in list_of_genomes:
            genome_gene_count = Counter(genome_content)
            for gene, count in genome_gene_count.items():
                if count == 1:  # single-copy only wanted.
                    gene_counts[gene] += 1
        # Here I completely trust the annotations from NCBI, so I just count their presence.
        # IF THE count is over the threshold, then we add it to the marker set.
        marker_sets = [
            gene for gene, count in gene_counts.items() if count >= threshold
        ]
        return marker_sets


def find_conserved_markers(descendants, species_genome_dict, pthreshold=0.97):
    # You've got all descendants, not find the gene set that gets conserved.
    # We shall change a bit.
    sp_gen_dict = species_genome_dict.copy()
    mod_dict, key_map = unify_species_keys(sp_gen_dict)
    list_of_genes = []
    num_genomes = len(descendants)
    threshold = pthreshold * num_genomes
    taxon_key_gene_counter = Counter()
    for i in descendants:
        gene_content_dict = iterative_match_species_to_key(i, mod_dict)
        list_of_genes.append(gene_content_dict["gene content"]["gene contents"])

    # breakpoint()
    for genome in list_of_genes:
        genome_gene_counter = Counter(genome)
        for gene, count in genome_gene_counter.items():
            if count == 1:
                taxon_key_gene_counter[gene] += 1

    conserved_marker_sets = [
        gene for gene, count in taxon_key_gene_counter.items() if count >= threshold
    ]
    # print(conserved_marker_sets)
    return conserved_marker_sets  # a set.


def get_all_descendants(taxon_key, tree):
    # Find all matching taxon.
    tree_obj = Tree(tree)
    leaf_names = [leaf.name for leaf in tree_obj.iter_leaves()]
    genomes = set([leaf for leaf in leaf_names if taxon_key in leaf])
    # A more sophisticated way: try to search nodes and see which nodes have majorly the taxon key (> 90%)
    for node in tree_obj.traverse():
        node_leaves = set(node.get_leaves())

        # How to determine A node that contains all taxon_key bearing genomes, and also cover something other than having taxon_key?
        # Or we add annotations to the tree first, so that every monophyletic group shall get one thing they wants.
        if genomes.issubset(node_leaves) and len(genomes) > 0.95 * len(node_leaves):
            node_descendants = node_leaves
        else:
            node_descendants = genomes
    return node_descendants


def taxon_indexing(leafname):
    # Input, a string with the type: phylum_class_order_family_genus_genus@species.
    lineage, _ = leafname.split("_")[:-1], leafname[-1]
    # Just like the encoding in NLP, we add each taxon group an index.
    index_taxon_dict = {}
    lineage_lst = lineage.split("_")
    for i in range(len(lineage_lst)):
        index_taxon_dict[i] = lineage_lst[i]
    return index_taxon_dict  # Returns the dict, and the length of lineage list, to make sure we don't lose anything.


def get_universal_markers(tree, species_genome_dict):
    # In this case, all genomes are included, so we add all genomes (i.e.: all leaves) to genomes included.
    tree_obj = Tree(tree)
    leaf_names = [leaf.name for leaf in tree_obj.iter_leaves()]
    sp_gen_dict = species_genome_dict.copy()
    mod_dict, _ = unify_species_keys(sp_gen_dict)
    list_of_genes = []

    # I want to find single-copy genes conserved in > 97% of the genomes.
    num_genomes = len(leaf_names)
    threshold = 0.97 * num_genomes
    gene_counts = Counter()
    for leafname in leaf_names:
        gene_content_dict = iterative_match_species_to_key(leafname, mod_dict)
        list_of_genes.append(gene_content_dict["gene content"]["gene contents"])

    for genome_content in list_of_genes:
        genome_gene_count = Counter(genome_content)
        for gene, count in genome_gene_count.items():
            if count == 1:
                gene_counts[gene] += 1
    # Here I completely trust the annotations from NCBI, so I just count their presence.
    marker_sets = [gene for gene, count in gene_counts.items() if count >= threshold]
    return {"marker set": marker_sets}


def unify_species_keys(species_genome_dict: dict):
    # Remove every non-alphanumerical characters.
    mod_dict = dict()
    key_map = dict()
    for key, value in species_genome_dict.items():
        new_key = re.sub(r"[^a-zA-Z0-9]", "", key)
        mod_dict[new_key] = value
        key_map[new_key] = key
    return mod_dict, key_map


def iterative_match_species_to_key(leaf_name: str, species_genome_dict: dict):
    # The leafname is commonly lineage_concatenated.
    # Again, we use the same idea: by adding each character from right end to left.
    # we assay whether this new character is present in the key bag of the species_genome_dict.
    # Once matched, we return a new dict storing all these items.

    # first, we identify anything inside the leaf name that are non-alphanumerical after splitted by "_".
    leaf_name_elements = leaf_name.split("_")
    leaf_name_2 = [re.sub(r"[^a-zA-Z0-9]", "", elem) for elem in leaf_name_elements]
    i = 0
    while i < len(leaf_name_2):
        leaf_try = "".join(leaf_name_2[i:])
        if leaf_try in species_genome_dict.keys():
            species = leaf_try
            gene_content = species_genome_dict[leaf_try]
            break
        else:
            i += 1
    # breakpoint()
    return {"species": species, "gene content": gene_content}


# Pooling taxa into a networkx tree?
def tree_taxa_pool(tree):
    # The whole tree's taxon will be here.
    leaves = tree.get_leaves()
    network = nx.Graph()
    for i in leaves:
        lineage = i.name.split("_")[:-1]
        for j in range(len(lineage)):
            if lineage[j] not in network.nodes and j == 0:
                network.add_node(lineage[j])
                # No counting for these.
            elif j >= 1 and lineage[j] not in network.nodes:
                network.add_node(lineage[j])
                network.add_edge(lineage[j], lineage[j - 1])
                network[lineage[j]]["count"] = 1
            else:
                network[lineage[j]]["count"] += 1
    # An unrooted network, each supergroup is isolated from each other.
    return network


# More functions.
def Calculate_node_marker_set(species_genome_dict, tree_obj=None, outfile=None):
    # For each of the node, get their descendants, calculate marker sets.
    tree = Tree(tree_obj)
    # We only want to have internal nodes.
    # Somehow: how to transform nodes to exact taxon? That can only be done by a statistical model.
    nodes = [node for node in tree.traverse() if not node.is_leaf()]
    index = 0
    for node in tree.traverse():
        if not node.is_leaf():
            Node_obj = TreeNode(node)
            marker_set = Node_obj.node_marker_sets(species_genome_dict)
            taxon_network = Node_obj.taxonomy_group_pooling()
            node_identifier = f"id{index}_{hex(id(node))}"
            print(
                f"marker set is: \n {marker_set}\n taxa stored in {hex(id(taxon_network))}\n node: {node_identifier}"
            )
            index += 1
            # save_marker_sets_taxon(marker_set, i, node_identifier)
            # Save the marker set into the node object, so that the newick-formatted tree has each node calculated with marker sets.
            node.marker_set = str("@".join(marker_set))
            # Once these things all get done, save these files into the directory: marker_sets.
        else:
            continue
    tree.write(
        outfile=outfile, features=["marker_set"], format=1
    )  # that is already backuped.
    return


def save_marker_sets_taxon(
    marker_set: list, taxon_node: ete3.coretype.tree.TreeNode, node_identifier: str
):
    # Here saves all marker sets and taxon_networks for each node.
    # First, create the directory under /data.
    Home_dir = os.path.join("/data", "nodes_marker_sets")
    if not os.path.exists(os.path.join("/data", "nodes_marker_sets")):
        os.mkdir(os.path.join("/data", "nodes_marker_sets"))
        # Then for each node a directory will be built to save both marker sets and network graph.
    Node_dir = os.path.join(os.path.join(Home_dir, node_identifier))
    if not os.path.exists(Node_dir):
        os.mkdir(Node_dir)
    else:
        os.rmdir(Node_dir)
        os.mkdir(Node_dir)
    UNIVERSAL_MARKER_SET_FILE = "marker_set.pkl"
    UNIVERSAL_NETWORK_FILE = "node.pkl"
    with open(os.path.join(Node_dir, UNIVERSAL_MARKER_SET_FILE), "wb") as ms:
        pkl.dump(marker_set, ms)
    with open(os.path.join(Node_dir, UNIVERSAL_NETWORK_FILE), "wb") as nf:
        pkl.dump(taxon_node, nf)
    return
