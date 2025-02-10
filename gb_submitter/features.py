#TODO: make output path if doesnt exist 
#TODO: add help info
#TODO: do not reverse complement when negarnaviricota ==> see taxonomy script(?)
import os

import click
import Bio.SeqIO
import pandas as pd
import pyrodigal_gv
from Bio.SeqIO import write
from gb_submitter import utils

def calculate_coding_capacity(genes, seq_length):
        """Calculate the total coding capacity for a list of genes."""
        return sum((gene.end - gene.begin) / seq_length for gene in genes)

def find_orientation(genes):
    """
    Calculate the sum of the strand orientations for a list of genes.

    Parameters:
    genes (list): A list of gene objects, each having a 'strand' attribute.

    Returns:
    int: The sum of strand orientations across all genes, where each strand is typically represented as 1 (forward) or -1 (reverse).
    """
    return sum(gene.strand for gene in genes)
    
def extract_gene_results(genes, record_id, seq_length):
    """Extract gene prediction results for a sequence."""
    return [
        [
            record_id,
            seq_length,
            f"{record_id}_{i+1}",
            gene.begin,
            gene.end,
            gene.strand,
            gene.start_node.type,
            gene.partial_begin,
            gene.partial_end
        ]
        for i, gene in enumerate(genes)
    ]
def write_proteins(genes, record_id, dst_path, overwrite):
    """Write protein translations to a file."""
    with open(dst_path, "w" if overwrite else "a") as dst:
        genes.write_translations(
            dst,
            sequence_id=f"{record_id}",
            width=80,
            translation_table=genes[0].translation_table,
            include_stop=False,
        )
    return False  # Update overwrite flag
def write_nucleotides(sequence, output_handle, overwrite):
    """Write protein translations to a file."""
    with open(output_handle, "w" if overwrite else "a") as dst:
        write(sequence, dst, "fasta")
    return False

def select_top_structure(df):
    """
    Select the top structure for each query based on the bitscore.

    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame with columns 'query', 'bits'.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with the top structure for each query.
    """
    highest_bits_idx = df.groupby("query")["bits"].idxmax()
    # Select those rows
    result = df.loc[highest_bits_idx]
    return result

# Function to generate and save NCBI feature tables
def save_ncbi_feature_tables(df, output_dir="./"):
    """
    Generate and save NCBI feature tables for sequences in a DataFrame.

    This function creates feature table files for each unique sequence ID
    in the input DataFrame and saves them in the specified output directory.
    Each feature table file is named using the accession number of the sequence.
    The content includes feature annotations such as CDS type, product name,
    and notes about ORF predictions and alternative start codons.

    Parameters:
    df (pd.DataFrame): DataFrame containing sequence data with columns
                       ['seqid', 'accession', 'start', 'end', 'strand', 'type',
                        'Protein names', 'source', 'start_codon', 'partial_begin', 'partial_end'].
    output_dir (str): Directory path to save the feature tables. Defaults to "./".

    """
    os.makedirs(output_dir, exist_ok=True)

    for seqid, group in df.groupby("seqid"):
        accession = group["seqid"].iloc[0]
        filename = f"{output_dir}{accession}.tbl"
        with open(filename, "w") as file:
            file.write(f">Feature {accession}\n")
            for _, row in group.iterrows():
                end = row['end']
                start = row['start']
                if row['partial_end']:
                    # This should not be possible with extended pyrodigal enforcing start codon (?)
                    if row['strand'] == -1:
                        end = f"<{row['end']}"
                    # This is fine
                    else:
                        end = f">{row['end']}"
                elif row['partial_begin']:
                    if row['strand'] == -1:
                        start = f">{row['start']}"

                if row["strand"] == -1:
                    file.write(f"{end}\t{start}\t{row['type']}\n")
                else:
                    file.write(f"{start}\t{end}\t{row['type']}\n")

                if pd.notna(row["Protein names"]):
                    protein = row["Protein names"]
                else:
                    protein = "hypothetical protein"

                file.write(f"\t\t\tproduct\t{protein}\n")
                file.write(f"\t\t\tinference\tab initio prediction:{row['source']}\n")

                if row["start_codon"] != "ATG":
                    file.write(
                        f"\t\t\tnote\tAlternative start codon: {row['start_codon']}\n"
                    )
                if protein != "hypothetical protein":
                    #file.write(f"\t\t\tnote\tAnnotation database: {row['annotation_source']}\n")
                    file.write((f"\t\t\tinference\talignment:{row['aligner']}:{row['aligner_version']}:UniProtKB:{row['Uniref_entry']},BFVD:{row['model']}\n"))

        print(f"Saved: {filename}")

@click.command(help="Create feature tables for sequences.")
@click.option('-i', '--input-file', 'fasta_file', required=True, type=click.Path(exists=True), help='Input fasta file')
@click.option('-o', '--output-path', 'output_path', required=True, type=click.Path(exists=False), help='Output directory')
@click.option('-d', '--database', 'database', required=True, type=click.Path(exists=True), help='BFVD diamond database path')
#TODO: limit genetic codes?
@click.option('-g', '--translation-table', 'transl_table', required=False, type=click.Choice(list(range(1, 32))), default=1, metavar="<1-31>", help='Translation table to use')
@click.option('-t', '--threads', 'threads', required=False, default=4, type=int, help='Number of threads to use')
def features(fasta_file, output_path, database, transl_table, threads):
    records = list(Bio.SeqIO.parse(fasta_file, "fasta"))

    orf_finder = pyrodigal_gv.ViralGeneFinder()
    training_info = orf_finder.train(*(bytes(seq.seq) for seq in records), translation_table=transl_table)


    orf_finder = pyrodigal_gv.ViralGeneFinder(
        meta=False, viral_only=True, closed=True, training_info=training_info
    )

    orf_finder2 = pyrodigal_gv.ViralGeneFinder(
        meta=False, viral_only=True, closed=[True,False], training_info=training_info
    )

    prot_path = f"{output_path}/proteins.faa"
    nucl_path = f"{output_path}/reoriented_nucleotide_sequences.fna"
    results = []
    no_orf_pred = []
    overwrite = True
    overwrite_n = True
    for record in records:
        seq_length = len(record.seq)

        seq_results = []

        genes = orf_finder.find_genes(bytes(record.seq))
        coding_capacity = calculate_coding_capacity(genes, seq_length)
        orientation = find_orientation(genes)

        if coding_capacity > 0.5:
            if orientation < 0:
                record.seq = record.seq.reverse_complement()
                genes = orf_finder.find_genes(bytes(record.seq))
            seq_results = extract_gene_results(genes, record.id, seq_length)
            overwrite = write_proteins(genes, record.id, prot_path, overwrite)
            overwrite_n = write_nucleotides(record, nucl_path, overwrite_n)
        else:
            print(f"Repredicting ORFs for {record.id} because of low coding capacity.")
            genes = orf_finder2.find_genes(bytes(record.seq))
            coding_capacity = calculate_coding_capacity(genes, seq_length)
            orientation = find_orientation(genes)

            if coding_capacity > 0.5:
                if orientation < 0:
                    record.seq = record.seq.reverse_complement()
                    genes = orf_finder2.find_genes(bytes(record.seq))
                seq_results = extract_gene_results(genes, record.id, seq_length)
                overwrite_n = write_nucleotides(record, nucl_path, overwrite_n)
                overwrite = write_proteins(genes, record.id, prot_path, overwrite)
            else:
                no_orf_pred.append(record.id)
                print(
                    f"No ORF predictions with a total coding capacity over 50% for {record.id}."
                )

        results.extend(seq_results)
    
    with open(f'{output_path}/no_ORF_prediction.txt', 'w') as f:
        for line in no_orf_pred:
            f.write(f"{line}\n")


    # Create DataFrame from results
    columns = ["seqid", "seq_length", "orf", "start", "end", "strand", "start_codon", "partial_begin", "partial_end"]
    df = pd.DataFrame(results, columns=columns)

    df["seqid"] = df["seqid"].str.strip()
    df["type"] = "CDS"
    df["source"] = f"pyrodigal-gv:{pyrodigal_gv.__version__}"
    #df["source"] = f"pyrodigal-gv"
    #df["annotation_source"]=f"BFVD (https://doi.org/10.1093/nar/gkae1119)"
    df["annotation_source"]="UniProtKB"

    #Cmd = "diamond blastp "
    #Cmd += f"--db {database}/foldseek_db/bfvd.dmnd "
    #Cmd += f"--query {output_path}/proteins.faa "
    #Cmd += f"--out {output_path}/gb_sub_proteins.m8 "
    #Cmd += f"--threads {threads} "
    #Cmd += "--sensitive "
    #Cmd += "--index-chunks 1 "
    #Cmd += "--block-size 8 "
    #Cmd += "--unal 1 "
    #Cmd += "--tmpdir /dev/shm "
    #Cmd += "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"
    #utils.Exec(Cmd)
#
    #aligner = "Diamond"
    #aligner_version = utils.Exec("diamond version")
    #aligner_version = aligner_version.strip().split()[2]

    Cmd = "mmseqs easy-search "
    Cmd += f"{output_path}/proteins.faa " #input
    Cmd += f"{database}/foldseek_db" # database TODO: fix path
    Cmd += f"{output_path}/gb_sub_proteins.m8 " #output
    Cmd += "tmp " #temp directory
    Cmd += "-s 7.5 "
    Cmd += "--format-mode 4 "
    Cmd += "--format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader "
    Cmd += f"--threads {threads}"

    aligner = "MMseqs2"
    aligner_version = utils.Exec("mmseqs version")

    m8 = pd.read_csv("gb_sub_proteins.m8", sep="\t", header=None)
    m8.rename(
        {
            0: "query",
            1: "target",
            2: "pident",
            3: "len",
            4: "mismatch",
            5: "gapopen",
            6: "qstart",
            7: "qend",
            8: "sstart",
            9: "send",
            10: "evalue",
            11: "bits",
            12: "stitle",
        },
        axis=1,
        inplace=True,
    )


    m8=m8[m8["evalue"]<1e-3]

    m8['aligner'] = aligner
    m8['aligner_version'] = aligner_version

    m8_top = select_top_structure(m8)
    names_df = pd.read_csv(f"{database}/bfvd_uniprot_names.tsv", sep="\t")

    # remove all trailing strings within brackets from protein names
    names_df['Protein names'] = names_df['Protein names'].str.replace(r'[\(\[].*?[\)\]]$', '', regex=True)

    meta_df = pd.read_csv(f"{database}/bfvd_metadata.tsv", sep="\t", header=None)
    meta_df.rename(
        {
            0: "Uniref_entry",
            1: "model",
            2: "length",
            3: "avg_pLDDT",
            4: "pTM",
            5: "splitted",
        },
        axis=1,
        inplace=True,
    )

    meta_df["model"] = meta_df["model"].str.replace(".pdb", "")

    merged_df = pd.merge(
        meta_df, names_df, left_on="Uniref_entry", right_on="From", how="left"
    )

    prot_df = pd.merge(m8_top, merged_df, left_on="target", right_on="model", how="left")

    #prot_df["Protein names"]

    prot_df.to_csv(f"{output_path}/diamond_names.tsv", sep="\t", index=False)

    diamond = prot_df
    final_df = pd.merge(df, diamond, left_on="orf", right_on="query", how="left")  
    
    # Call the function to save feature tables
    save_ncbi_feature_tables(
        final_df,
        output_dir=f"{output_path}/feature_tables/",
    )

if __name__ == "__main__":
    features()