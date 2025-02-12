# TODO: add help info
import os

import Bio.SeqIO
import click
import pandas as pd
import pyrodigal_gv
import taxopy
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
            gene.partial_end,
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
    """Write nucleotide sequences to a file."""
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


def read_taxonomy_table(taxonomy):
    taxonomy_data = None
    if taxonomy:
        taxonomy_data = pd.read_csv(taxonomy, sep="\t")
    return taxonomy_data


def predict_orfs(orf_finder, seq):
    """Find genes, compute coding capacity, and determine orientation."""
    genes = orf_finder.find_genes(bytes(seq))
    coding_capacity = calculate_coding_capacity(genes, len(seq))
    orientation = find_orientation(genes)
    return genes, coding_capacity, orientation, orf_finder


def get_lineage(record_id, taxonomy_data, taxdb):
    """Retrieve the lineage of a given record from the taxonomy table."""
    record_taxonomy = taxonomy_data[taxonomy_data["contig"] == record_id]
    if record_taxonomy.empty:
        return []
    taxid_dict = record_taxonomy.set_index("contig")["taxid"].to_dict()
    lineage = taxopy.Taxon(taxid_dict[record_id], taxdb).name_lineage
    return lineage


# Function to generate and save NCBI feature tables
# Function to generate and save NCBI feature tables
def save_ncbi_feature_tables(df, output_dir="./", single_file=True):
    """
    Generate and save NCBI feature tables for sequences in a DataFrame.

    This function creates a single feature table file by default, but can
    also save separate files for each unique sequence ID when specified.

    Parameters:
    df (pd.DataFrame): DataFrame containing sequence data with columns
                       ['seqid', 'accession', 'start', 'end', 'strand', 'type',
                        'Protein names', 'source', 'start_codon', 'partial_begin', 'partial_end'].
    output_dir (str): Directory path to save the feature tables. Defaults to "./".
    single_file (bool): If True, saves all features to one file; otherwise, saves separate files.
    """

    if single_file:
        filename = os.path.join(output_dir, "featuretable.tbl")
        with open(filename, "w") as file:
            for seqid, group in df.groupby("seqid"):
                accession = group["seqid"].iloc[0]
                file.write(f">Feature {accession}\n")
                write_feature_entries(file, group)
        click.echo(f"Saved: {filename}")

    else:
        os.makedirs(os.path.join(output_dir, "feature_tables"), exist_ok=True)
        for seqid, group in df.groupby("seqid"):
            accession = group["seqid"].iloc[0]
            filename = os.path.join(output_dir, "feature_tables", f"{accession}.tbl")
            with open(filename, "w") as file:
                file.write(f">Feature {accession}\n")
                write_feature_entries(file, group)
            click.echo(f"Saved: {filename}")


def write_feature_entries(file, group):
    """Helper function to write feature entries to a file."""
    for _, row in group.iterrows():
        start, end = row["start"], row["end"]

        if row["partial_end"]:
            end = f"<{row['end']}" if row["strand"] == -1 else f">{row['end']}"
        if row["partial_begin"] and row["strand"] == -1:
            start = f">{row['start']}"

        file.write(
            f"{end}\t{start}\t{row['type']}\n"
            if row["strand"] == -1
            else f"{start}\t{end}\t{row['type']}\n"
        )

        protein = (
            row["Protein names"]
            if pd.notna(row["Protein names"])
            else "hypothetical protein"
        )
        file.write(f"\t\t\tproduct\t{protein}\n")
        file.write(f"\t\t\tinference\tab initio prediction:{row['source']}\n")

        if row["start_codon"] != "ATG":
            file.write(f"\t\t\tnote\tAlternative start codon: {row['start_codon']}\n")
        if protein != "hypothetical protein":
            file.write(
                f"\t\t\tinference\talignment:{row['aligner']}:{row['aligner_version']}:UniProtKB:{row['Uniref_entry']},BFVD:{row['model']}\n"
            )


@click.command(help="Create feature tables for sequences.")
@click.option(
    "-i",
    "--input-file",
    "fasta_file",
    required=True,
    type=click.Path(exists=True),
    help="Input fasta file",
)
@click.option(
    "-o",
    "--output-path",
    "output_path",
    required=True,
    type=click.Path(exists=False),
    help="Output directory",
)
@click.option(
    "-d",
    "--database",
    "database",
    required=True,
    type=click.Path(exists=True),
    help="BFVD diamond database path",
)
# TODO: limit genetic codes?
@click.option(
    "-g",
    "--translation-table",
    "transl_table",
    required=False,
    type=click.Choice(list(range(1, 32))),
    default=1,
    metavar="<1-31>",
    help="Translation table to use",
)
@click.option(
    "--taxonomy",
    required=False,
    type=click.Path(exists=True),
    help="Taxonomy file",
)
@click.option(
    "--separate-files",
    required=False,
    is_flag=True,
    help="Save feature tables into separate files",
)
@click.option(
    "-t",
    "--threads",
    "threads",
    required=False,
    default=4,
    type=int,
    help="Number of threads to use",
)
def features(
    fasta_file, output_path, database, transl_table, taxonomy, separate_files, threads
):
    if os.path.exists(output_path):
        click.echo(
            f"Warning: Output directory '{output_path}' already exists and may be overwritten."
        )

    os.makedirs(output_path, exist_ok=True)

    records = list(Bio.SeqIO.parse(fasta_file, "fasta"))

    # Train ORF finder
    orf_finder = pyrodigal_gv.ViralGeneFinder()
    training_info = orf_finder.train(
        *(bytes(seq.seq) for seq in records), translation_table=transl_table
    )

    # Initialize ORF finders
    orf_finder1 = pyrodigal_gv.ViralGeneFinder(
        meta=False, viral_only=True, closed=True, training_info=training_info
    )
    orf_finder2 = pyrodigal_gv.ViralGeneFinder(
        meta=False, viral_only=True, closed=[True, False], training_info=training_info
    )

    # Load taxonomy database
    taxdb = taxopy.TaxDb(
        nodes_dmp="/lustre1/scratch/337/vsc33750/ictv_db/ictv_taxdump/nodes.dmp",
        names_dmp="/lustre1/scratch/337/vsc33750/ictv_db/ictv_taxdump/names.dmp",
    )  # TODO: Set database path
    taxonomy_data = read_taxonomy_table(taxonomy)

    # Define output paths
    prot_path = f"{output_path}/proteins.faa"
    nucl_path = f"{output_path}/reoriented_nucleotide_sequences.fna"

    results, no_orf_pred = [], []
    overwrite, overwrite_n = True, True

    for record in records:
        lineage = get_lineage(record.id, taxonomy_data, taxdb) if taxonomy else []

        # Predict ORFs using orf_finder1 first
        genes, coding_capacity, orientation, chosen_orf_finder = predict_orfs(
            orf_finder1, record.seq
        )

        # If coding capacity is too low, use orf_finder2 instead
        if coding_capacity <= 0.5:
            # click.echo(f"Repredicting ORFs for {record.id} due to low coding capacity.")
            genes, coding_capacity, orientation, chosen_orf_finder = predict_orfs(
                orf_finder2, record.seq
            )

        if coding_capacity > 0.5:
            # Adjust orientation based on lineage
            if (orientation < 0 and "Negarnaviricota" not in lineage) or (
                orientation > 0 and "Negarnaviricota" in lineage
            ):
                record.seq = record.seq.reverse_complement()
                genes, _, _, _ = predict_orfs(
                    chosen_orf_finder, record.seq
                )  # Use the last used ORF finder

            results.extend(extract_gene_results(genes, record.id, len(record.seq)))
            overwrite = write_proteins(genes, record.id, prot_path, overwrite)
            overwrite_n = write_nucleotides(record, nucl_path, overwrite_n)
        else:
            no_orf_pred.append(record.id)
            # click.echo(
            #    f"No ORF predictions with start site and >50% coding capacity for {record.id}."
            # )

    with open(f"{output_path}/no_ORF_prediction.txt", "w") as f:
        for line in no_orf_pred:
            f.write(f"{line}\n")

    # Create DataFrame from results
    columns = [
        "seqid",
        "seq_length",
        "orf",
        "start",
        "end",
        "strand",
        "start_codon",
        "partial_begin",
        "partial_end",
    ]
    df = pd.DataFrame(results, columns=columns)

    df["seqid"] = df["seqid"].str.strip()
    df["type"] = "CDS"
    df["source"] = f"pyrodigal-gv:{pyrodigal_gv.__version__}"
    # df["source"] = f"pyrodigal-gv"
    # df["annotation_source"]=f"BFVD (https://doi.org/10.1093/nar/gkae1119)"
    df["annotation_source"] = "UniProtKB"

    # Cmd = "diamond blastp "
    # Cmd += f"--db {database}/foldseek_db/bfvd.dmnd "
    # Cmd += f"--query {output_path}/proteins.faa "
    # Cmd += f"--out {output_path}/gb_sub_proteins.m8 "
    # Cmd += f"--threads {threads} "
    # Cmd += "--sensitive "
    # Cmd += "--index-chunks 1 "
    # Cmd += "--block-size 8 "
    # Cmd += "--unal 1 "
    # Cmd += "--tmpdir /dev/shm "
    # Cmd += "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    # utils.Exec(Cmd)
    #
    # aligner = "Diamond"
    # aligner_version = utils.Exec("diamond version", capture=True)
    # aligner_version = aligner_version.strip().split()[2]

    Cmd = "mmseqs easy-search "
    Cmd += f"{output_path}/proteins.faa "  # input
    Cmd += f"{database} "  # database
    Cmd += f"{output_path}/gb_sub_proteins.m8 "  # output
    Cmd += "tmp "  # temp directory
    Cmd += "-s 7.5 "
    Cmd += "--format-mode 4 "
    Cmd += "--format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits "
    Cmd += f"--threads {threads}"

    utils.Exec(Cmd)

    aligner = "MMseqs2"
    aligner_version = utils.Exec("mmseqs version", capture=True).strip()

    m8 = pd.read_csv(
        f"{output_path}/gb_sub_proteins.m8", sep="\t", header=None, low_memory=False
    )
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
            8: "tstart",
            9: "tend",
            10: "evalue",
            11: "bits",
        },
        axis=1,
        inplace=True,
    )

    m8["evalue"] = pd.to_numeric(m8["evalue"], errors="coerce")
    m8 = m8[m8["evalue"] < 1e-3]

    m8["aligner"] = aligner
    m8["aligner_version"] = aligner_version

    m8_top = select_top_structure(m8)
    names_df = pd.read_csv(
        f"{database}_uniprot_names.tsv", sep="\t"
    )  # TODO find better solution

    # remove all trailing strings within brackets from protein names
    names_df["Protein names"] = names_df["Protein names"].str.replace(
        r"[\(\[].*?[\)\]]$", "", regex=True
    )

    meta_df = pd.read_csv(
        f"{database}_metadata.tsv", sep="\t", header=None
    )  # TODO find better solution
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

    prot_df = pd.merge(
        m8_top, merged_df, left_on="target", right_on="model", how="left"
    )

    # prot_df["Protein names"]

    prot_df.to_csv(
        f"{output_path}/diamond_names.tsv", sep="\t", index=False  # TODO change name
    )

    diamond = prot_df
    final_df = pd.merge(df, diamond, left_on="orf", right_on="query", how="left")

    single_file = False if separate_files else True
    # Call the function to save feature tables
    save_ncbi_feature_tables(
        final_df, output_dir=f"{output_path}", single_file=single_file
    )


if __name__ == "__main__":
    features()
