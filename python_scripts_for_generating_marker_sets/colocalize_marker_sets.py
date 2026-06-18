from collections import defaultdict
from Bio import SeqIO
import os
import networkx as nx


def parse_genbank(genbank_file):
    """
    Parses a GenBank file and returns a list of gene information.
    Each gene will be represented as a dictionary with the gene name and its location.
    """
    genes = []

    # Parse the GenBank file using SeqIO
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "gene":
                # Extract gene name and location (start and end positions)
                gene_name = feature.qualifiers.get("gene", [""])[0]
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand  # 1 for forward, -1 for reverse

                if gene_name:
                    genes.append(
                        {
                            "name": gene_name,
                            "start": start,
                            "end": end,
                            "strand": strand,
                        }
                    )
    # This is summarized information, containing all genes in the genome, need selection.
    return genes


def extract_marker_positions(genes, marker_genes):
    """
    Extracts the positions of marker genes from the parsed gene list.
    """
    return {
        gene["name"]: (gene["start"], gene["end"])
        for gene in genes
        if gene["name"] in marker_genes
    }


def check_colocalization(marker_positions, max_distance=5000):
    """
    Checks if marker genes are colocalized within a certain genomic distance.
    """
    if len(marker_positions) < 2:
        return (
            False  # Colocalization is not relevant if less than 2 markers are present
        )

    # Sort positions by start location
    positions = sorted(marker_positions.values(), key=lambda x: x[0])

    for i in range(len(positions) - 1):
        start1, end1 = positions[i]
        start2, end2 = positions[i + 1]

        # Check the distance between the end of the first gene and the start of the next gene
        if start2 - end1 > max_distance:
            return False

    return True


# Does it have to be true that we need to minimize max distance for plastid genomes, like reducing it to 3500/3000/2000?
# The key is that: length do impact the colocalizations: if too long, everything will be linked.
def build_marker_gene_set_graph(marker_positions, max_distance=5000):
    G = nx.Graph()

    sorted_marker_genes = sorted(marker_positions.items(), key=lambda x: x[1][0])

    for i in range(len(sorted_marker_genes)):
        gene1, (start1, end1) = sorted_marker_genes[i]
        for j in range(i + 1, len(sorted_marker_genes)):
            gene2, (start2, end2) = sorted_marker_genes[j]

            # Colocalized within 5000 bps.
            if start2 - end1 <= max_distance:
                G.add_edge(gene1, gene2)
            else:
                break

    # Return the graph, one graph for each genome.
    return G


def find_colocalized_clusters(marker_positions, max_distance=5000):

    G = build_marker_gene_set_graph(marker_positions, max_distance)
    marker_clusters = list(nx.connected_components(G))

    return marker_clusters


def process_genomes_colocalization(genomes_marker_positions, max_distance=5000):
    """
    Processes multiple genomes and tracks the frequency of colocalized gene pairs.

    genomes_marker_positions: A list where each entry corresponds to a genome, containing a dictionary
    of gene positions for conserved marker genes. Format: [{gene_name: (start, end)}, ...]
    """
    colocalization_counter = defaultdict(int)
    total_genomes = len(genomes_marker_positions)

    # Process each genome individually
    for marker_positions in genomes_marker_positions:
        # Get clusters of colocalized genes for the current genome
        clusters = find_colocalized_clusters(marker_positions, max_distance)

        # For each cluster, record gene colocalization
        for cluster in clusters:
            cluster = sorted(list(cluster))  # Sort to ensure consistent ordering
            # Increment colocalization count for each pair in the cluster
            for i in range(len(cluster)):
                for j in range(i + 1, len(cluster)):
                    pair = (cluster[i], cluster[j])
                    colocalization_counter[pair] += 1

    return colocalization_counter, total_genomes


# Based on CheckM.
def define_lineage_specific_marker_sets(
    colocalization_counter, total_genomes, threshold=0.95
):
    """
    Defines lineage-specific marker sets by identifying sets of genes that are colocalized in >95% of genomes.
    """
    # Find pairs that colocalize in >95% of the genomes
    threshold_count = total_genomes * threshold

    # Build the final graph for colocalized gene sets
    G = nx.Graph()

    # Add edges for gene pairs that meet the colocalization threshold
    for pair, count in colocalization_counter.items():
        if count >= threshold_count:
            G.add_edge(pair[0], pair[1])

    # Find connected components (colocalized gene sets) in the final graph
    marker_sets = list(nx.connected_components(G))

    return marker_sets


def map_species_to_locus_tags(selected_genomes: list, species_to_locus_tags: dict):
    """
    Maps species names to locus tags for each gene in the GenBank files.
    """
    loci = []
    for genome in selected_genomes:
        lineage, genome_species = genome.split("@")[0], genome.split("@")[1:]
        genus = lineage.split("_")[-1]
        loci.append(species_to_locus_tags[genus + "@" + "@".join(genome_species)])

    return loci


def map_locus_to_genbank(locus_list: list, genbank_files: list):
    """
    Now the locus tags should be mapped to genbank files.
    """
    matched_genbank = []
    for i in locus_list:
        for genbank_file in genbank_files:
            if i in genbank_file:
                matched_genbank.append(genbank_file)
    return matched_genbank


def Colocalize_marker_sets(
    picked_leaves: list, marker_set: list, locus_to_species: dict, genbank_files: list
):

    # Step 1: match the selected genomes to the locus tags.
    selected_genomes = map_species_to_locus_tags(picked_leaves, locus_to_species)
    matched_genbank = map_locus_to_genbank(selected_genomes, genbank_files)

    # Step 2: parse the GenBank file and extract marker gene positions.
    genomes_marker_positions = []
    for genbank_file in matched_genbank:
        genes = parse_genbank(genbank_file)
        marker_positions = extract_marker_positions(genes, marker_set)
        genomes_marker_positions.append(
            marker_positions
        )  # a list of dictionaries, each dictionary contains marker gene positions for one genome.

    # Step 3: using marker positions to check colocalization.
    colocalization_counter, total_genomes = process_genomes_colocalization(
        genomes_marker_positions, max_distance=5000
    )

    # Step 4: define lineage-specific marker sets.
    lineage_specific_marker_sets = define_lineage_specific_marker_sets(
        colocalization_counter, total_genomes, threshold=0.95
    )

    # Now we have the lineage-specific marker sets, write it in correct order. We want it to be a line of string.
    left_genes = set(marker_set).difference(set.union(*lineage_specific_marker_sets))

    Final_marker_set = ""
    Final_marker_set += "["
    # First, add grouped marker genes as sets into this list.
    # Add comma between sets, also, add single quotation marks for each gene.
    for i, marker_set in enumerate(lineage_specific_marker_sets):
        Final_marker_set += "set(["
        Final_marker_set += ", ".join([f"'{gene}'" for gene in marker_set]) + "])"
        if i < len(lineage_specific_marker_sets) - 1:
            Final_marker_set += ", "

    # Then, add individual genes in terms of individual set into this list.
    if left_genes:
        for gene in left_genes:
            Final_marker_set += ", set(["
            Final_marker_set += f"'{gene}'" + "])"

    Final_marker_set += "]"
    return Final_marker_set
