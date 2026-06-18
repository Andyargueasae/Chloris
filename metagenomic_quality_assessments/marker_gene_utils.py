import os
import sys
import subprocess
import pandas as pd
import networkx as nx
import pyfastx
import shutil

# Take a look at LorBin, COMEBin and CompleteBin's use of marker genes scripts for a reference.
THRESHOLD=1e-10 # for hmmsearch, the e-value threshold.

def orf_runner(contig_file: str, output_dir: str, threads: int = 10, predictor: str = "fraggenescan") -> str:
    # try run fraggenescan.
    if predictor == "fraggenescan":
        FragGeneScan_cmd = f"{shutil.which('FragGeneScan')}"
        FragGeneScan_params = f"-s {contig_file} -o {output_dir}/contigs_Fr -w 0 -t complete -p {threads}"

        cmd = subprocess.Popen(f"{FragGeneScan_cmd} {FragGeneScan_params}", shell=True)
        cmd.communicate()

        if cmd.returncode != 0:
            print("Error running FragGeneScan", file=sys.stderr)
    else:
        # use prodigal as default
        Prodigal_cmd = "prodigal"
        Prodigal_params = f"-i {contig_file} -p meta -q -m -a {output_dir}/contigs_Pr.faa"
        cmd = subprocess.Popen(f"{Prodigal_cmd} {Prodigal_params}", shell=True)
        cmd.communicate()

        if cmd.returncode != 0:
            print("Error running Prodigal", file=sys.stderr)

    return output_dir + f"/contigs_{'Fr' if predictor == 'fraggenescan' else 'Pr'}.faa"

def hmmer_runner(orf_fasta: str, output_dir: str, hmm_file: str, threads: int = 10):
    # try run hmmsearch.
    HMMER_cmd = "hmmsearch"
    outfile = os.path.join(output_dir, "hmmsearch_output.domtbl")
    HMMER_params = f"--cpu {threads} --domtblout {outfile} -E {THRESHOLD} {hmm_file} {orf_fasta}"
    if os.path.exists(outfile):
        os.remove(outfile)
    try:
        with open(f"{output_dir}/hmmsearch_output.log", "w") as log_file:
            cmd = subprocess.Popen(f"{HMMER_cmd} {HMMER_params}", shell=True, stdout=log_file, stderr=log_file)
            cmd.communicate()
            if cmd.returncode != 0:
                print("Error running hmmsearch", file=sys.stderr)
    except:
        print("Error running hmmsearch", file=sys.stderr)
    return outfile

# Three functions above has a core command to call. Now what below is to parse the output of hmmsearch.

def get_marker_hits(hmmsearch_output_file: str, contig_names: list, fasta_file: str, min_contig_len=1000,
                    predictor="fraggenescan") -> dict:
    # This is after running hmmsearch, parse the output file to get marker hits.
    sequence2markers = {}
    # turn the hmmsearch output into a table, pandas.
    cols = ['target_name', 'target_accession', 'tlen', 'query_name', 'query_accession', 'qlen',
            'E_value', 'score', 'bias', 'domE_value', 'dom_score', 'dom_bias',
            'hmm_from', 'hmm_to', 'ali_from', 'ali_to', 'env_from', 'env_to',
            'acc', 'description_of_target']
    df = pd.read_csv(hmmsearch_output_file, comment='#', sep=r"\s+", header=None,
                     usecols=[0, 3, 5, 15, 16], names=["orf", "gene", "qlen", "qstart", "qend"])
    if df.empty:
        return sequence2markers
    
    def contig_name(ell):
        if predictor == "fraggenescan":
            return ell.rsplit('_', 3)[0]
        else:
            return ell.split('_', 1)[0]
    # filter by E-value threshold
    df['contig'] = df['orf'].map(lambda x: contig_name(x))
    data = df.query('(qend - qstart) / qlen > 0.4').copy()  # filter by length of hit: at least 40% of HMM length 
    contig_name_set = set(contig_names)

    contig_len_dict = {c:len(seq) for c, seq in pyfastx.Fasta(fasta_file, build_index=False)}
    data = data[data['contig'].map(lambda c: contig_len_dict[c] >= min_contig_len)]

    data = data.drop_duplicates(['gene', 'contig'])
    marker = data['gene'].values
    contig = data['contig'].values

    for m, c in zip(marker, contig):
        if c not in sequence2markers:
            sequence2markers[c] = []
        sequence2markers[c].append(m)

    return sequence2markers

def load_marker_set(marker_set_file: str):
    lineage_dict = {0: 'domain', 1: 'phylum', 2: 'class', 3: 'order', 4: 'family', 5: 'genus', 6: 'species'}
    tms_data = nx.DiGraph()
    with open(marker_set_file, "r") as f:
        for line in f:
            line = line.strip('\n \t').split('\t')
            marker_sets = eval(line[-1])
            marker_sets = [{marker.split('.')[0] for marker in marker_set} for marker_set in
                           marker_sets]
            lineage = line[2].split(';')
            lineage = [tax + '[' + lineage_dict[tax_ind] + ']' if tax == lineage[tax_ind - 1] and len(lineage) > 1
                       else tax for tax_ind, tax in enumerate(lineage)]
            if len(lineage) == 7:
                lineage[-1] = ' '.join(lineage[-2:])

            tms_data.add_nodes_from(lineage)
            tms_data.add_edges_from([(i, lineage[index + 1]) for index, i in enumerate(lineage)
                                         if not index + 1 == len(lineage)])
            tms_data.nodes[lineage[-1]]['markers'] = int(line[4])
            tms_data.nodes[lineage[-1]]['marker_groups'] = int(line[5])
            tms_data.nodes[lineage[-1]]['marker_sets'] = marker_sets
    return tms_data

def get_marker_list_node_quality(marker_list: list, node: str, tms_data: nx.DiGraph) -> tuple:
    node_marker_sets = tms_data.nodes.data()[node]['marker_sets']

    node_marker_sets_completenesses = []
    node_marker_sets_purities = []

    for marker_set in node_marker_sets:
        marker_set_stats = get_marker_set_quality(marker_set, marker_list)
        if marker_set_stats:
            node_marker_sets_completenesses.append(marker_set_stats[0])
            if marker_set_stats[1]:
                node_marker_sets_purities.append(marker_set_stats[1])
        else:
            node_marker_sets_completenesses.append(0)

    node_marker_set_completeness = round(sum(node_marker_sets_completenesses)
                                         / len(node_marker_sets_completenesses), 3)
    if node_marker_sets_purities:
        node_marker_set_purity = round(sum(node_marker_sets_purities)
                                       / len(node_marker_sets_purities), 3)
    else:
        node_marker_set_purity = 0

    if node_marker_set_completeness > 1:
        print('Completeness of for marker set {0} is > 1 with {1} for'
                      ' marker list {2}'.format(node, node_marker_set_completeness,
                                                marker_list))
        raise Exception

    return [node_marker_set_completeness, node_marker_set_purity]

def get_marker_set_quality(marker_set: set, marker_list: list) -> tuple:

    marker_set_markers_found = [marker for marker in marker_list
                                if marker in marker_set]
    if not marker_set_markers_found:
        return [0.0, 0.0]
    node_markers_list_set = set(marker_set_markers_found)

    n_set_markers_found = marker_set.intersection(node_markers_list_set)

    marker_set_completeness = round(len(n_set_markers_found) / len(marker_set), 3)

    if marker_set_completeness > 1:
        marker_set_completeness = 1

    marker_set_marker_purities = [round(1 / marker_set_markers_found.count(marker), 3)
                                  # for marker in set(node_t2p_markers_list_set).intersection(set(marker_set_t2p_markers_found))]
                                  for marker in node_markers_list_set]
    try:
        marker_set_average_purity = round(sum(marker_set_marker_purities)
                                          / len(marker_set_marker_purities), 3)
    except ZeroDivisionError:
        print(marker_set_markers_found)
        print(node_markers_list_set)
        print(marker_set_marker_purities)
        raise ZeroDivisionError

    return [marker_set_completeness, marker_set_average_purity]

def compare_marker_set_stats(marker_set, current_best_marker_set, completeness_variability):
    if marker_set[1] >= current_best_marker_set[1] * completeness_variability:
        current_best_marker_set = marker_set
    return current_best_marker_set

def choose_marker_sets_on_quality(marker_list: list, tms_data: nx.DiGraph, max_depth_lvl: int = 4) -> dict:
    nodes = [n for n, d in tms_data.in_degree() if d == 0]
    current_node = nodes[0]
    previous_nodes = None
    best_marker_set = []
    depth_grace_count = 0
    king_lvl_stats = None
    current_depth_level = 0
    while list(tms_data[current_node]) and depth_grace_count < 2 and current_depth_level <= max_depth_lvl:
        current_level_best_marker_set = []
        if previous_nodes == nodes:
            depth_grace_count += 1
            nodes = [sub_node for node in nodes for sub_node in list(tms_data[node])]
        previous_nodes = nodes
        for index, node in enumerate(nodes):
            if isinstance(tms_data.nodes.data()[node]['marker_sets'], str):
                print('Marker set of {0} identical to higher level set {1}.'
                              ' Skipping.'.format(node,
                            tms_data.nodes.data()[node]['marker_sets'].split('_')[1]))
                if not best_marker_set:
                    best_marker_set = [node, 0.00, 0.00, 0.00, current_depth_level]
                if not current_level_best_marker_set:
                    current_level_best_marker_set = [node, 0.00, 0.00, 0.00, current_depth_level]
                continue
            node_n_markers = tms_data.nodes.data()[node]['markers']
            node_n_marker_sets = tms_data.nodes.data()[node]['marker_groups']
            node_stats = get_marker_list_node_quality(marker_list, node, tms_data)
            if not node_stats:
                continue
            node_marker_set_completeness = node_stats[0]
            node_marker_set_purity = node_stats[1]
            node_marker_set_completeness_score = round(node_marker_set_completeness
                                    * node_n_marker_sets / node_n_markers * 100, 3)
            if king_lvl_stats: # and node not in ['Bacteria', 'Archaea']:
                current_marker_set = [node, ((node_marker_set_completeness + king_lvl_stats[0]) / 2),
                                      ((node_marker_set_purity + king_lvl_stats[1]) / 2),
                                      node_marker_set_completeness_score, current_depth_level]
            else:
                current_marker_set = [node, node_marker_set_completeness, node_marker_set_purity,
                                      node_marker_set_completeness_score, current_depth_level]

            if not best_marker_set: # or (best_marker_set[0] in ['Bacteria', 'Archaea'] and current_depth_level > 0):
                best_marker_set = [node, node_marker_set_completeness, node_marker_set_purity,
                                   node_marker_set_completeness_score, current_depth_level]
            else:
                # Check if comparing same level sets, if so dont give completeness leeway
                if current_marker_set[-1] == best_marker_set[-1]:
                    completeness_variability = 1.0
                else:
                    completeness_variability = 0.975
                best_marker_set = compare_marker_set_stats(current_marker_set,
                                                              best_marker_set, completeness_variability)
            if not current_level_best_marker_set:  #  or (current_level_best_marker_set[0] in ['Bacteria', 'Archaea']
                                                   # and current_depth_level > 0):
                current_level_best_marker_set = [node, node_marker_set_completeness, node_marker_set_purity,
                                                 node_marker_set_completeness_score, current_depth_level]
            else:
                if current_marker_set[-1] == current_level_best_marker_set[-1]:
                    completeness_variability = 1.0
                else:
                    completeness_variability = 0.975
                current_level_best_marker_set = compare_marker_set_stats(current_marker_set,
                                                                         current_level_best_marker_set,
                                                                         completeness_variability)
        nodes = list(tms_data[current_level_best_marker_set[0]])

        current_node = current_level_best_marker_set[0]

        if current_level_best_marker_set[0] in ['Bacteria', 'Archaea']:
            king_lvl_stats = [node_marker_set_completeness, node_marker_set_purity]

        current_depth_level += 1

    if best_marker_set:
        return best_marker_set
    
    else:
        print('Something went wrong while choosing the best marker set. Markers:'
                      ' {0}; unique: {1}; total {2}.'.format(set(marker_list),
                                                             len(set(marker_list)), len(marker_list)))
        return ['None', 0, 0, 0]

    return

def single_contig_bins(contig_dict, fasta_file, output_dir, marker_gene_dir, minfasta=0, latent=None) -> list:
    '''
    Retrieve single-contig bins from the contig dictionary based on the quality assessed by marker sets. 
    
    :param contig_dict: Contig dictionary mapping contig names to sequences
    :param contig_list: List of contig names to consider
    :param marker_gene_dir: Directory containing marker gene files for quality assessment, lineage by lineage.
    :param minfasta: Minimum length of contig to be considered a bin
    :param latent: Latent space representation of contigs (used to filter scMAGs and leave only unselected ones)
    :return: List of single-contig bins

    # Used inside HDBSCAN, BIRCH. 
    '''
    marker_set_tsv = os.path.join(marker_gene_dir, "taxon_marker_sets_lineage_sorted.tsv")
    hmm_file = os.path.join(marker_gene_dir, "markers.hmm")

    # Step 1, given a contig file, run orf finder to find ORFs. 
    orf_fasta = orf_runner(fasta_file, output_dir=output_dir, predictor="fraggenescan")
    # Step 2, run hmmsearch to find marker genes.
    hmmsearch_output = hmmer_runner(orf_fasta, output_dir=output_dir, hmm_file=hmm_file)
    # Step 3, parse the hmmsearch output to get marker hits.
    contig_names = list(contig_dict.keys())

    sequence2markers = get_marker_hits(hmmsearch_output, contig_names, fasta_file, min_contig_len=minfasta, predictor="fraggenescan")
    # Step 4, load marker set data.
    tms_data = load_marker_set(marker_set_tsv)
    # Step 5, assess quality of each contig based on marker hits. Also take a look at latent space to filter out scMAGs.
    single_contig_bins = []

    # print something in prior.
    print("Maximum number of markers on contig:", max(len(v) for v in sequence2markers.values()))
    for contig, marker_gene in sequence2markers.items():
        marker_set = choose_marker_sets_on_quality(marker_gene, tms_data)
        taxon, com, pur = marker_set[0], marker_set[1], marker_set[2]
        com = 100*com
        pur = 100*pur
        if contig in contig_dict:
            print(f"Contig: {contig}, Taxon: {taxon}, Completeness: {com}, Purity: {pur}, Length: {len(contig_dict[contig])}")
        if com >= 90 and pur >= 92.5:
            single_contig_bins.append(contig)
    
    # filter that in the latent embeddings, simply select out those that have been filtered.
    # assume latent is a pd.dataframe with index as contig names.
    if latent is not None:
        latent = latent.drop(index=single_contig_bins)

    return single_contig_bins, latent.to_numpy()
