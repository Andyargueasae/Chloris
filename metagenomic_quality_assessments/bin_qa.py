import typer
import os
import pandas as pd
import tempfile
import shutil
import subprocess
import sys
import pyfastx
from marker_gene_utils import load_marker_set, choose_marker_sets_on_quality

THRESHOLD=1e-10

app = typer.Typer()

def orf_runner(contig_file: str, output_dir: str, threads: int = 10, predictor: str = "fraggenescan") -> str:
    # try run fraggenescan.
    basename = os.path.basename(contig_file).rsplit('.', 1)[0]

    if predictor == "fraggenescan":
        FragGeneScan_cmd = f"{shutil.which('FragGeneScan')}"
        FragGeneScan_params = f"-s {contig_file} -o {output_dir}/{basename}_Fr -w 0 -t complete -p {threads}"

        cmd = subprocess.Popen(f"{FragGeneScan_cmd} {FragGeneScan_params}", shell=True)
        cmd.communicate()

        if cmd.returncode != 0:
            print("Error running FragGeneScan", file=sys.stderr)
    else:
        # use prodigal as default
        Prodigal_cmd = "prodigal"
        Prodigal_params = f"-i {contig_file} -p meta -q -m -a {output_dir}/{basename}_Pr.faa"
        cmd = subprocess.Popen(f"{Prodigal_cmd} {Prodigal_params}", shell=True)
        cmd.communicate()

        if cmd.returncode != 0:
            print("Error running Prodigal", file=sys.stderr)

    return output_dir + f"/{basename}_{'Fr' if predictor == 'fraggenescan' else 'Pr'}.faa"

def hmmer_runner(orf_fasta: str, output_dir: str, hmm_file: str, threads: int = 10, use_tc: bool = False) -> str:
    # try run hmmsearch.
    basename = os.path.basename(orf_fasta).rsplit('.', 1)[0]
    HMMER_cmd = "hmmsearch"
    outfile = os.path.join(output_dir, f"{basename}_hmmsearch_output.domtbl")

    HMMER_params = f"--cpu {threads} --domtblout {outfile} -E {THRESHOLD} {hmm_file} {orf_fasta}"
    if use_tc:
        HMMER_params = f"--cpu {threads} --domtblout {outfile} --cut_tc {hmm_file} {orf_fasta}"
    if os.path.exists(outfile):
        os.remove(outfile)
    try:
        with open(f"{output_dir}/{basename}_hmmsearch_output.log", "w") as log_file:
            cmd = subprocess.Popen(f"{HMMER_cmd} {HMMER_params}", shell=True, stdout=log_file, stderr=log_file)
            cmd.communicate()
            if cmd.returncode != 0:
                print("Error running hmmsearch", file=sys.stderr)
    except:
        print("Error running hmmsearch", file=sys.stderr)
    return outfile

def get_marker_hits(hmmsearch_output_file: str,
                    predictor="fraggenescan", orf_align_prop: float = 0.60) -> dict:
    # This is after running hmmsearch, parse the output file to get marker hits.
    bin2markers = {}
    basename = os.path.basename(hmmsearch_output_file).rsplit('.', 1)[0]
    # turn the hmmsearch output into a table, pandas.
    cols = ['target_name', 'target_accession', 'tlen', 'query_name', 'query_accession', 'qlen',
            'E_value', 'score', 'bias', 'domE_value', 'dom_score', 'dom_bias',
            'hmm_from', 'hmm_to', 'ali_from', 'ali_to', 'env_from', 'env_to',
            'acc', 'description_of_target']
    df = pd.read_csv(hmmsearch_output_file, comment='#', sep=r"\s+", header=None,
                     usecols=[0, 3, 4, 5, 15, 16], names=["orf", "gene", "accession", "qlen", "qstart", "qend"])
    breakpoint()
    if df.empty:
        return bin2markers
    
    def contig_name(ell):
        if predictor == "fraggenescan":
            return ell.rsplit('_', 3)[0]
        else:
            return ell.split('_', 1)[0]
    # filter by E-value threshold
    df['contig'] = df['orf'].map(lambda x: contig_name(x))
    data = df.query(f'(qend - qstart) / qlen > {orf_align_prop}').copy()  # filter by length of hit: at least 40% of HMM length 
    data = data.drop_duplicates(['accession', 'contig'])
    marker = data['accession'].values
    contig = data['contig'].values
    for m, c in zip(marker, contig):
        if c not in bin2markers:
            bin2markers[c] = []
        bin2markers[c].append(m)

    # concatenate markers for each contig into a long list.
    final_markers = [m for markers in bin2markers.values() for m in markers]

    return {basename: final_markers}

@app.command()
def assess_quality(bin_dir: str, marker_gene_dir: str, output_file: str, suffix: str="fa", predictor: str="fraggenescan",
                   use_tc: bool = False, orf_align_prop: float = 0.60):
    if os.listdir(bin_dir) == []:
        print("No bins found in the specified directory.")
        return
    
    marker_set_tsv = os.path.join(marker_gene_dir, "taxon_marker_sets_lineage_sorted.tsv")
    bins = [os.path.join(bin_dir, f) for f in os.listdir(bin_dir) if os.path.isfile(os.path.join(bin_dir, f)) and (f.endswith(f'.{suffix}'))]

    tmp_dir = os.path.join(os.path.dirname(output_file), "tmp_marker_genes")
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    if os.path.exists(output_file):
        os.remove(output_file)
    
    with open(output_file, "a") as out_f:
        out_f.write("Bin\ttaxonomy\tCompleteness\tPurity\n")
    
    for bin in bins:
        orf_bin = orf_runner(bin, output_dir=tmp_dir, predictor=predictor)
        hmmsearch_output = hmmer_runner(orf_bin, output_dir=tmp_dir, 
                                        hmm_file=os.path.join(marker_gene_dir, "markers.hmm"), use_tc=use_tc)
        # parse hmmsearch output and write to output_file
        # hold on, this is the marker list for the whole bin, not per contig.
        bin2markers = get_marker_hits(hmmsearch_output, predictor=predictor, orf_align_prop=orf_align_prop)
        # use this, to estimate the quality. 
        tms_data = load_marker_set(marker_set_tsv)
        with open(output_file, "a") as out_f:
            bin_name, markers = bin2markers.popitem()
            marker_set = choose_marker_sets_on_quality(markers, tms_data)
            taxon, com, pur = marker_set[0], marker_set[1], marker_set[2]
            com = 100*com
            pur = 100*pur
            out_f.write(f"{bin_name}\t{taxon}\t{com:.2f}\t{pur:.2f}\n")
    # shutil.rmtree(tmp_dir)
    print(f"Bin quality assessment completed. Results saved to {output_file}.")
    return


if __name__ == "__main__":
    app()