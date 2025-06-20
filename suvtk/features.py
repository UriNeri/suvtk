"""
features.py
===========

This script processes input sequences to predict open reading frames (ORFs),
aligns the predicted protein sequences against a database, and generates feature
tables for submission to GenBank.

Functions
---------
validate_translation_table(ctx, param, value)
    Validate the given translation table.

calculate_coding_capacity(genes, seq_length)
    Calculate the total coding capacity for a list of genes.

find_orientation(genes)
    Determine the orientation of genes based on strand information.

predict_orfs(orf_finder, seq)
    Predict ORFs, compute coding capacity, and determine orientation.

features(fasta_file, output_path, database, transl_table, coding_complete, taxonomy, separate_files, threads)
    Main command to create feature tables for sequences.
"""

# TODO: add count for sequences without ORF prediction
# TODO: mmseqs log to file for clarity
import os
import shutil

import needletail
import click
import polars as pl
import pyrodigal_gv
import taxopy
from Bio.SeqIO import write

from suvtk import utils

# Define the valid genetic codes based on NCBI
VALID_GENETIC_CODES = set(range(1, 7)) | set(range(9, 17)) | set(range(21, 32))


def validate_translation_table(ctx, param, value):
    """
    Validate that the given translation table is one of the valid genetic codes.

    Parameters
    ----------
    ctx : click.Context
        The Click context object. Unused.
    param : click.Parameter
        The parameter object. Unused.
    value : int
        The given translation table.

    Returns
    -------
    int
        The given translation table if it is valid.

    Raises
    ------
    click.BadParameter
        If the given translation table is not valid.
    """
    if value not in VALID_GENETIC_CODES:
        raise click.BadParameter(
            f"Invalid translation table. Must be one of {sorted(VALID_GENETIC_CODES)}."
        )
    return value


def calculate_coding_capacity(genes, seq_length):
    """
    Calculate the total coding capacity for a list of genes.

    Parameters
    ----------
    genes : list
        A list of gene objects.
    seq_length : int
        The length of the sequence.

    Returns
    -------
    float
        The total coding capacity.
    """
    return sum((gene.end - gene.begin) / seq_length for gene in genes)


def find_orientation(genes):
    """
    Calculate the sum of the strand orientations for a list of genes.
    If the sum is zero, return the orientation of the largest gene.

    Parameters
    ----------
    genes : list
        A list of gene objects, each having 'strand', 'begin', and 'end' attributes.

    Returns
    -------
    int
        The sum of strand orientations across all genes, or the orientation of the largest gene if the sum is zero.
    """
    orientation_sum = sum(gene.strand for gene in genes)

    if orientation_sum == 0 and genes:
        # Find the largest gene based on absolute length (|end - begin|)
        largest_gene = max(genes, key=lambda gene: abs(gene.end - gene.begin))
        return largest_gene.strand

    return orientation_sum


def extract_gene_results(genes, record_id, seq_length):
    """
    Extract gene prediction results for a sequence.

    Parameters
    ----------
    genes : list
        A list of gene objects.
    record_id : str
        The ID of the sequence record.
    seq_length : int
        The length of the sequence.

    Returns
    -------
    list
        A list of gene prediction results.
    """
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
    """
    Write protein translations to a file.

    Parameters
    ----------
    genes : list
        A list of gene objects.
    record_id : str
        The ID of the sequence record.
    dst_path : str
        The destination file path.
    overwrite : bool
        Whether to overwrite the file.

    Returns
    -------
    bool
        Updated overwrite flag.
    """
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
    """
    Write nucleotide sequences to a file.

    Parameters
    ----------
    sequence : dict
        The sequence record with 'id' and 'seq' keys.
    output_handle : str
        The output file path.
    overwrite : bool
        Whether to overwrite the file.

    Returns
    -------
    bool
        Updated overwrite flag.
    """
    with open(output_handle, "w" if overwrite else "a") as dst:
        dst.write(f">{sequence['id']}\n")
        # Write sequence with line breaks every 80 characters
        seq = sequence['seq']
        for i in range(0, len(seq), 80):
            dst.write(seq[i:i+80] + "\n")
    return False


def select_top_structure(df):
    """
    Select the top structure for each query based on the bitscore.

    Parameters
    ----------
    df : polars.DataFrame
        A DataFrame with columns 'query' and 'bits'.

    Returns
    -------
    polars.DataFrame
        A DataFrame with the top structure for each query.
    """
    # In polars, we use window functions to get the row with maximum bits per group
    result = df.filter(
        pl.col("bits") == pl.col("bits").max().over("query")
    )
    return result


def predict_orfs(orf_finder, seq):
    """
    Find genes, compute coding capacity, and determine orientation.

    Parameters
    ----------
    orf_finder : pyrodigal_gv.ViralGeneFinder
        The ORF finder object.
    seq : str
        The sequence to analyze.

    Returns
    -------
    tuple
        A tuple containing genes, coding capacity, orientation, and the ORF finder used.
    """
    genes = orf_finder.find_genes(bytes(seq))
    coding_capacity = calculate_coding_capacity(genes, len(seq))
    orientation = find_orientation(genes)
    return genes, coding_capacity, orientation, orf_finder


def get_lineage(record_id, taxonomy_data, taxdb):
    """
    Retrieve the lineage of a given record from the taxonomy table.

    Parameters
    ----------
    record_id : str
        The ID of the sequence record.
    taxonomy_data : polars.DataFrame
        The taxonomy data table.
    taxdb : taxopy.TaxDb
        The taxonomy database.

    Returns
    -------
    list
        The lineage of the record.
    """
    record_taxonomy = taxonomy_data.filter(pl.col("contig") == record_id)
    if record_taxonomy.height == 0:
        return []
    # taxid_dict = record_taxonomy.set_index("contig")["taxid"].to_dict()
    tax = record_taxonomy["taxonomy"].item().removesuffix(" sp.")
    if tax == "unclassified viruses":
        return []
    try:
        taxid = taxopy.taxid_from_name(tax, taxdb)
        lineage = taxopy.Taxon(taxid[0], taxdb).name_lineage
        return lineage
    except IndexError:
        click.echo(
            f"Warning: '{tax}' is not part of the official ICTV taxonomy. Its lineage can not be looked up and therefore the nucleotide reorientation could not be performed for {record_taxonomy['contig'].item()}"
        )
        return []


# Functions to generate and save NCBI feature tables
def save_ncbi_feature_tables(df, output_dir=".", single_file=True):
    """
    Generate and save NCBI feature tables for sequences in a DataFrame.

    This function creates a single feature table file by default, but can
    also save separate files for each unique sequence ID when specified.

    Parameters
    ----------
    df : polars.DataFrame
        DataFrame containing sequence data with columns
        ['seqid', 'accession', 'start', 'end', 'strand', 'type',
        'Protein names', 'source', 'start_codon', 'partial_begin', 'partial_end'].
    output_dir : str, optional
        Directory path to save the feature tables. Defaults to ".".
    single_file : bool, optional
        If True, saves all features to one file; otherwise, saves separate files.

    Returns
    -------
    None
    """

    if single_file:
        filename = os.path.join(output_dir, "featuretable.tbl")
        with open(filename, "w") as file:
            for seqid, group in df.group_by("seqid"):
                accession = group["seqid"].item(0)
                file.write(f">Feature {accession}\n")
                write_feature_entries(file, group)
        click.echo(f"Saved: {filename}")

    else:
        os.makedirs(os.path.join(output_dir, "feature_tables"), exist_ok=True)
        for seqid, group in df.group_by("seqid"):
            accession = group["seqid"].item(0)
            filename = os.path.join(output_dir, "feature_tables", f"{accession}.tbl")
            with open(filename, "w") as file:
                file.write(f">Feature {accession}\n")
                write_feature_entries(file, group)
            click.echo(f"Saved: {filename}")


def write_feature_entries(file, group):
    """
    Helper function to write feature entries to a file.

    Parameters
    ----------
    file : file-like object
        The file to write to.
    group : polars.DataFrame
        The group of feature entries to write.

    Returns
    -------
    None
    """
    for row in group.iter_rows(named=True):
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
            if row["Protein names"] is not None
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


@click.command(short_help="Create feature tables for sequences.")
@click.option(
    "-i",
    "--input",
    "fasta_file",
    required=True,
    type=click.Path(exists=True),
    help="Input fasta file",
)
@click.option(
    "-o",
    "--output",
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
    help="Path to the suvtk database folder.",
)
@click.option(
    "-g",
    "--translation-table",
    "transl_table",
    required=False,
    type=int,
    callback=validate_translation_table,
    default=1,
    metavar="",
    help="Translation table to use. Only genetic codes from https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi are allowed (1-6, 9-16, 21-31).",
)
@click.option(
    "--coding-complete",
    required=False,
    is_flag=True,
    help="Do not predict incomplete genes (no stop codon) and only keep genomes that are 'coding complete' (>50% coding capacity). [This can not be turned off for now]",  # TODO: check on pyrodigal to implement fixed start codon
)
@click.option(
    "--taxonomy",
    required=False,
    type=click.Path(exists=True),
    help="Taxonomy file to adjust sequence orientation (ssRNA- sequences will get 3' -> 5' orientation, all others 5' -> 3').",
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
    default=utils.get_available_cpus(),
    type=int,
    metavar="",
    help="Number of threads to use",
)
def features(
    fasta_file,
    output_path,
    database,
    transl_table,
    coding_complete,
    taxonomy,
    separate_files,
    threads,
):
    """
    Create feature tables for sequences from an input fasta file.

    This command processes the input sequences to predict open reading frames (ORFs),
    aligns the predicted protein sequences against a specified database with proteins and their function, and generates
    feature tables for submission to GenBank.
    """

    if os.path.exists(output_path):
        click.echo(
            f"Warning: Output directory '{output_path}' already exists and will be overwritten."
        )

    os.makedirs(output_path, exist_ok=True)

    # Parse FASTA file with needletail
    records = []
    for record in needletail.parse_fastx_file(fasta_file):
        records.append({
            'id': record.id,
            'seq': record.seq
        })

    # Train ORF finder
    orf_finder = pyrodigal_gv.ViralGeneFinder()
    training_info = orf_finder.train(
        *(bytes(seq['seq'], 'utf-8') for seq in records), translation_table=transl_table
    )

    # Initialize ORF finders
    orf_finder1 = pyrodigal_gv.ViralGeneFinder(
        meta=False, viral_only=True, closed=True, training_info=training_info
    )

    # Not possible for now because pyrodigal does not support 'closed=[True, False]' yet
    # orf_finder2 = pyrodigal_gv.ViralGeneFinder(
    #    meta=False, viral_only=True, closed=[True, False], training_info=training_info
    # )

    # Load taxonomy database
    if taxonomy:
        taxdb = taxopy.TaxDb(
            nodes_dmp=os.path.join(database, "nodes.dmp"),
            names_dmp=os.path.join(database, "names.dmp"),
        )  # TODO: Set better database path?
        taxonomy_data = utils.safe_read_csv(taxonomy, sep="\t")

    # Define output paths
    prot_path = os.path.join(output_path, "proteins.faa")
    nucl_path = os.path.join(output_path, "reoriented_nucleotide_sequences.fna")

    results, no_orf_pred = [], []
    overwrite, overwrite_n = True, True

    for record in records:
        lineage = get_lineage(record['id'], taxonomy_data, taxdb) if taxonomy else []

        # Predict ORFs using orf_finder1 first
        genes, coding_capacity, orientation, chosen_orf_finder = predict_orfs(
            orf_finder1, record['seq']
        )

        # Commented out because it is not possible to use orf_finder2 yet
        # If coding capacity is too low, use orf_finder2 instead
        # if coding_capacity < 0.5 and not coding_complete:
        #    # click.echo(f"Repredicting ORFs for {record['id']} due to low coding capacity.")
        #    genes, coding_capacity, orientation, chosen_orf_finder = predict_orfs(
        #        orf_finder2, record['seq']
        #    )

        if coding_capacity >= 0.5:
            # Adjust orientation based on lineage
            if (orientation < 0 and "Negarnaviricota" not in lineage) or (
                orientation > 0 and "Negarnaviricota" in lineage
            ):
                # Reverse complement the sequence
                record['seq'] = needletail.reverse_complement(record['seq'])
                genes, _, _, _ = predict_orfs(
                    chosen_orf_finder, record['seq']
                )  # Use the last used ORF finder

            results.extend(extract_gene_results(genes, record['id'], len(record['seq'])))
            overwrite = write_proteins(genes, record['id'], prot_path, overwrite)
            overwrite_n = write_nucleotides(record, nucl_path, overwrite_n)
        else:
            no_orf_pred.append(record['id'])
            # click.echo(
            #    f"No ORF predictions with start site and >50% coding capacity for {record.id}."
            # )

    with open(os.path.join(output_path, "no_ORF_prediction.txt"), "w") as f:
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
    df = pl.DataFrame(results, schema=columns)

    feat_pred = f"pyrodigal-gv"
    feat_pred_version = f"{pyrodigal_gv.__version__}"

    df = df.with_columns([
        pl.col("seqid").str.strip_chars(),
        pl.lit("CDS").alias("type"),
        pl.lit(f"{feat_pred}:{feat_pred_version}").alias("source"),
        pl.lit("UniProtKB").alias("annotation_source")
    ])

    # Cmd = "diamond blastp "
    # Cmd += f"--db {database}/foldseek_db/bfvd.dmnd "
    # Cmd += f"--query {output_path}/proteins.faa "
    # Cmd += f"--out {output_path}/alignment.m8 "
    # Cmd += "--threads {threads} "
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

    m8_path = os.path.join(output_path, "alignment.m8")

    Cmd = "mmseqs easy-search "
    Cmd += f"{prot_path} "  # input
    Cmd += os.path.join(database, "bfvd") + " "  # database
    Cmd += f"{m8_path} "  # output
    Cmd += "tmp "  # temp directory
    Cmd += "-s 7.5 "
    Cmd += "--format-mode 0 "
    Cmd += "--format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits "
    Cmd += f"--threads {threads}"
    utils.Exec(Cmd)

    shutil.rmtree("tmp")

    aligner = "MMseqs2"
    aligner_version = utils.Exec("mmseqs version", capture=True).strip()

    m8 = utils.safe_read_csv(
        m8_path,
        separator="\t",
        has_header=False,
    )
    
    # Rename columns in polars
    column_mapping = {
        "column_1": "query",
        "column_2": "target", 
        "column_3": "pident",
        "column_4": "len",
        "column_5": "mismatch",
        "column_6": "gapopen",
        "column_7": "qstart",
        "column_8": "qend",
        "column_9": "tstart", 
        "column_10": "tend",
        "column_11": "evalue",
        "column_12": "bits",
    }
    m8 = m8.rename(column_mapping)

    m8 = m8.filter(pl.col("evalue") < 1e-3)

    m8 = m8.with_columns([
        pl.lit(aligner).alias("aligner"),
        pl.lit(aligner_version).alias("aligner_version")
    ])

    m8_top = select_top_structure(m8)
    names_df = utils.safe_read_csv(
        os.path.join(database, "bfvd_uniprot_names.tsv"), separator="\t"
    )  # TODO find better solution?

    # remove all trailing strings within brackets from protein names
    names_df = names_df.with_columns(
        pl.col("Protein names").str.replace_all(r"[\(\[].*?[\)\]]$", "")
    )

    meta_df = utils.safe_read_csv(
        os.path.join(database, "bfvd_metadata.tsv"), separator="\t", has_header=False
    )  # TODO find better solution?
    
    # Rename columns in polars
    meta_column_mapping = {
        "column_1": "Uniref_entry",
        "column_2": "model",
        "column_3": "length",
        "column_4": "avg_pLDDT",
        "column_5": "pTM",
        "column_6": "splitted",
    }
    meta_df = meta_df.rename(meta_column_mapping)

    meta_df = meta_df.with_columns(
        pl.col("model").str.replace(".pdb", "")
    )

    merged_df = meta_df.join(
        names_df, left_on="Uniref_entry", right_on="From", how="left"
    )

    prot_df = m8_top.join(
        merged_df, left_on="target", right_on="model", how="left"
    )

    # prot_df["Protein names"]

    prot_df.write_csv(
        os.path.join(output_path, "tophit_info.tsv"),
        separator="\t",
        include_header=True
    )

    diamond = prot_df
    final_df = df.join(diamond, left_on="orf", right_on="query", how="left")

    single_file = False if separate_files else True
    # Call the function to save feature tables
    save_ncbi_feature_tables(
        final_df, output_dir=f"{output_path}", single_file=single_file
    )

    with open(os.path.join(output_path, "miuvig_features.tsv"), "w") as file:
        file.write(f"MIUVIG_parameter\tvalue\n")
        file.write(
            f"feat_pred\t{feat_pred};{feat_pred_version};-g {transl_table}, default otherwise\n"
        )
        file.write(
            f"ref_db\tBFVD;2023_02;https://bfvd.steineggerlab.workers.dev\n"
        )  # TODO: read DB version from version.txt or something
        file.write(
            f"sim_search_meth\t{aligner};{aligner_version};-s 7.5, default otherwise\n"
        )


if __name__ == "__main__":
    features()
