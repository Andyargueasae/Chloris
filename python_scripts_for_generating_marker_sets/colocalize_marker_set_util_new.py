from collections import defaultdict
from typing import List
from Bio import SeqIO
import os
import re
import networkx as nx
import pandas as pd


# 2025 July 4 - 7.
def locate_pfam_proteins(
    pfam_df, protein_df, genbank_ids: List[str], marker_set: List[str]
) -> List[dict]:
    marker_positions = list()

    # Note: marker set may contain both pfam and protein IDs, should split in prior.
    # PF: for PFAM, proteins are just the protein ID.
    pfam_pattern = re.compile(r"PF\d+\.\d+")
    # Only check pfam IDs.
    pfam_marker_set = [gene for gene in marker_set if pfam_pattern.match(gene)]
    protein_marker_set = [gene for gene in marker_set if not pfam_pattern.match(gene)]

    for genbank_id in genbank_ids:
        # each genome has a dictionary to store marker positions.
        genome_marker_positions = defaultdict(dict)
        # only one genome at a time,
        selected_genomes_pfam = pfam_df[
            (pfam_df["genome"] == genbank_id)
            & (pfam_df["query_accession"].isin(marker_set))
        ]
        selected_genomes_protein = protein_df[
            (protein_df["genome"] == genbank_id)
            & (protein_df["query_name"].isin(marker_set))
        ]

        # proteins are consensus.
        if len(protein_marker_set) > 0:
            # for each row (one hit), extract the gene name and its start and end positions.
            for gene, df in selected_genomes_protein.groupby("query_name"):
                if len(df) == 1:
                    # exactly single hit.
                    row = df.iloc[0]
                    gene_name = row["query_name"]
                    start = row["start_pos"]
                    end = row["end_pos"]
                    genome_marker_positions[gene_name] = (start, end)
                else:
                    # more than one hit, then in the marker set we just make it an individual set in marker set,
                    # skip calculating gene positions.
                    continue

        if len(pfam_marker_set) > 0:
            for gene, df in selected_genomes_pfam.groupby("query_accession"):
                if len(df) == 1:
                    # exactly single hit.
                    row = df.iloc[0]
                    gene_name = row["query_accession"]
                    start = row["start_pos"]
                    end = row["end_pos"]
                    genome_marker_positions[gene_name] = (start, end)
                else:
                    # more than one hit, then in the marker set we just make it an individual set in marker set,
                    # skip calculating gene positions.
                    continue
        # if genbank_id in ["AY286123", "LC704715", "LC716139", "LC716140"]:
        #     breakpoint()
        # Add the genome's marker positions to the list.
        marker_positions.append(genome_marker_positions)
    return marker_positions


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


def Colocalize_marker_sets(genbank_ids: list, marker_set: list, hmm_gene_table: tuple):

    # Step 1: read the HMM gene table and extract hmm_pfam and hmm_protein.
    hmm_pfam, hmm_protein = hmm_gene_table
    pfam_df = pd.read_csv(hmm_pfam, sep="\t", header=0)
    protein_df = pd.read_csv(hmm_protein, sep="\t", header=0)
    pfam_df["genome"] = pfam_df["genome"].apply(
        lambda x: x.split(".")[0]
    )  # Extract genbank id.
    protein_df["genome"] = protein_df["genome"].apply(lambda x: x.split(".")[0])

    # genbank ids are already in arguments.
    # Step 2: parse the tables and extract marker positions for both protein and pfams in the marker set.
    genomes_marker_positions = locate_pfam_proteins(
        pfam_df, protein_df, genbank_ids, marker_set
    )  # a list of dictionaries.
    # Note: some genes/pfam domains may be multi-copy in some genomes, so we make them individual sets in the marker set.
    # Step 3: using marker positions to check colocalization.
    colocalization_counter, total_genomes = process_genomes_colocalization(
        genomes_marker_positions, max_distance=5000
    )

    # Step 4: define lineage-specific marker sets.
    lineage_specific_marker_sets = define_lineage_specific_marker_sets(
        colocalization_counter, total_genomes, threshold=0.95
    )

    print(
        f"Total marker sets: {len(marker_set)}, colocalized marker sets: {len(lineage_specific_marker_sets)}"
    )
    # Now we have the lineage-specific marker sets, write it in correct order. We want it to be a line of string.
    try:
        left_genes = set(marker_set).difference(
            set.union(*lineage_specific_marker_sets)
        )
    except:
        # If no colocalized marker sets, then left_genes will be the whole marker set.
        left_genes = set(marker_set)
    Final_marker_set = [marker_set for marker_set in lineage_specific_marker_sets]
    Final_marker_set.extend(
        list([{left_gene} for left_gene in left_genes])
    )  # add the left genes that are not colocalized.
    return Final_marker_set
