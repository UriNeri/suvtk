import click
import pandas as pd
from gb_submitter import utils

@click.command(help="Assign virus taxonomy to sequences.")
@click.option('-i', '--input-file', 'fasta_file', required=True, type=click.Path(exists=True), help='Input fasta file')
@click.option('-o', '--output-path', 'output_path', required=True, type=click.Path(exists=False), help='Output directory')
@click.option('-d', '--database', 'database', required=True, type=click.Path(exists=True), help='ICTV MMseqs database path')
@click.option('-s', '--identity', 'seqid', required=False, default=0.3, type=float, help='Minimum sequence identity for hits to be considered')
@click.option('-t', '--threads', 'threads', required=False, default=4, type=int, help='Number of threads to use')
def taxonomy(fasta_file, database, output_path, seqid, threads):
    # Add RAM restrictions?
    # Add error handling
    Cmd = 'mmseqs easy-taxonomy '
    Cmd += f'{fasta_file} '
    Cmd += f'{database} '
    Cmd += f'{output_path}/taxresults '
    Cmd += 'tmp '
    Cmd += '--blacklist "" --tax-lineage 1 '
    Cmd += f'--min-seq-id {seqid} '
    Cmd += f'--threads {threads}'
    utils.Exec(Cmd)

    taxonomy = pd.read_csv(f"{output_path}/taxresults_lca.tsv", sep="\t", header=None)

    tax_names = []
    for index, row in taxonomy.iterrows():
        if row[2] == "no rank":
            print(f"No taxonomy for {row[0]}")
            last_known = "unknown"
        elif row[2] == "species": # Fix issue when species contains sp.
            last_known = row[8].split(";")[-2].replace("g_", "")
            last_known += " sp."
        else:
            last_known = row[3].strip()
            last_known += " sp."
        tax_names.append([row[0], last_known])

    tax_df = pd.DataFrame(tax_names, columns=['contig', 'taxonomy'])
    tax_df.to_csv(f'{output_path}/taxonomy.tsv', sep='\t', index=False)        

if __name__ == "__main__":
    taxonomy()