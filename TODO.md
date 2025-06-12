- **File:** download_database.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/download_database.py
  - TODO: If the archive in zenodo is in .zip format, it's content will be listed in the web view. Just a nice thing to have, if the compression ratio isn't detrimental.
  - Context: """
  - Scope: def find_tar_file(files: list) -> dict:
  - Modified: 2025-06-12 15:45:25.572044

- **File:** features.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/features.py
  - TODO: add count for sequences without ORF prediction
  - Context: 
  - Scope: main
  - Modified: 2025-06-12 14:14:06.052994

- **File:** features.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/features.py
  - TODO: mmseqs log to file for clarity
  - Context: # TODO: add count for sequences without ORF prediction
  - Scope: main
  - Modified: 2025-06-12 14:14:06.052994

- **File:** features.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/features.py
  - TODO: check on pyrodigal to implement fixed start codon
  - Context: is_flag=True,
  - Scope: def write_feature_entries(file, group):
  - Modified: 2025-06-12 14:14:06.052994

- **File:** features.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/features.py
  - TODO: Set better database path?
  - Context: names_dmp=os.path.join(database, "names.dmp"),
  - Scope: def features(
  - Modified: 2025-06-12 14:14:06.052994

- **File:** features.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/features.py
  - TODO: find better solution?
  - Context: os.path.join(database, "bfvd_uniprot_names.tsv"), separator="\t"
  - Scope: def features(
  - Modified: 2025-06-12 14:14:06.052994

- **File:** features.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/features.py
  - TODO: find better solution?
  - Context: os.path.join(database, "bfvd_metadata.tsv"), separator="\t", has_header=False
  - Scope: def features(
  - Modified: 2025-06-12 14:14:06.052994

- **File:** features.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/features.py
  - TODO: read DB version from version.txt or something
  - Context: f"ref_db\tBFVD;2023_02;https://bfvd.steineggerlab.workers.dev\n"
  - Scope: def features(
  - Modified: 2025-06-12 14:14:06.052994

- **File:** gbk2tbl.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/gbk2tbl.py
  - TODO: generally improve script
  - Context: # adapted from https://github.com/wanyuac/BINF_toolkit/blob/master/gbk2tbl.py
  - Scope: main
  - Modified: 2025-06-12 16:28:41.342243

- **File:** gbk2tbl.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/gbk2tbl.py
  - TODO: figure out a way to read gbk without darn biopython
  - Context: )  # read a GenBank file from the standard input and convert it into a list of SeqRecord objects
  - Scope: def gbk2tbl(input, mincontigsize, prefix):
  - Modified: 2025-06-12 16:28:41.342243

- **File:** table2asn.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/table2asn.py
  - TODO: add date correction option
  - Context: 
  - Scope: main
  - Modified: 2025-06-11 14:11:15.045690

- **File:** table2asn.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/table2asn.py
  - TODO: check/display stats file for errors (https://www.ncbi.nlm.nih.gov/genbank/validation/#BioSourceMissing) -> OKish
  - Context: # TODO: add date correction option
  - Scope: main
  - Modified: 2025-06-11 14:11:15.045690

- **File:** table2asn.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/table2asn.py
  - TODO: add missing required miuvig params from src file? -> comments.py?:  env_broad_scale, env_local_scale, env_medium, investigation_type, project_name, seq_meth
  - Context: # TODO: check/display stats file for errors (https://www.ncbi.nlm.nih.gov/genbank/validation/#BioSourceMissing) -> OKish
  - Scope: main
  - Modified: 2025-06-11 14:11:15.045690

- **File:** table2asn.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/table2asn.py
  - TODO: add option to make genbank file -> made by default
  - Context: # TODO: add missing required miuvig params from src file? -> comments.py?:  env_broad_scale, env_local_scale, env_medium, investigation_type, project_name, seq_meth
  - Scope: main
  - Modified: 2025-06-11 14:11:15.045690

- **File:** table2asn.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/table2asn.py
  - TODO: add check for required columns (src, comments)
  - Context: # TODO: add option to make genbank file -> made by default
  - Scope: main
  - Modified: 2025-06-11 14:11:15.045690

- **File:** taxonomy.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/taxonomy.py
  - TODO: save MIUVIG file with pred_genome_type and pred_genome_struc
  - Context: 
  - Scope: main
  - Modified: 2025-06-12 16:28:51.495373

- **File:** taxonomy.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/taxonomy.py
  - TODO: Change to genomad taxonomy
  - Context: # TODO: save MIUVIG file with pred_genome_type and pred_genome_struc
  - Scope: main
  - Modified: 2025-06-12 16:28:51.495373

- **File:** taxonomy.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/taxonomy.py
  - TODO: mmseqs overwrite tmp file (mmseqs fails when command was previously aborted) ||| why not use --remove-tmp-files ? shouldn't unless you intend to try and reuse the tmp stuff from a failed run (--force-reuse)
  - Context: # TODO: Change to genomad taxonomy
  - Scope: main
  - Modified: 2025-06-12 16:28:51.495373

- **File:** taxonomy.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/taxonomy.py
  - TODO: add memory size check and then manually set the  --split INT and --split-memory-limit BYTE so that at leas 200mb of the system availalbe memory are not used.
  - Context: # TODO: mmseqs overwrite tmp file (mmseqs fails when command was previously aborted) ||| why not use --remove-tmp-files ? shouldn't unless you intend to try and reuse the tmp stuff from a failed run (--force-reuse)
  - Scope: main
  - Modified: 2025-06-12 16:28:51.495373

- **File:** taxonomy.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/taxonomy.py
  - TODO: Add RAM restrictions?
  - Context: 
  - Scope: def taxonomy(fasta_file, database, output_path, seqid, threads):
  - Modified: 2025-06-12 16:28:51.495373

- **File:** taxonomy.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/taxonomy.py
  - TODO: Add error handling
  - Context: # TODO Add RAM restrictions?
  - Scope: def taxonomy(fasta_file, database, output_path, seqid, threads):
  - Modified: 2025-06-12 16:28:51.495373

- **File:** taxonomy.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/taxonomy.py
  - TODO: removing tmp before running mmseqs might be dangerous
  - Context: # TODO Add error handling
  - Scope: def taxonomy(fasta_file, database, output_path, seqid, threads):
  - Modified: 2025-06-12 16:28:51.495373

- **File:** taxonomy.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/taxonomy.py
  - TODO: check best cutoff
  - Context: elif (
  - Scope: def taxonomy(fasta_file, database, output_path, seqid, threads):
  - Modified: 2025-06-12 16:28:51.495373

- **File:** virus_info.py
  - Path: /home/neri/Documents/GitHub/suvtk/suvtk/virus_info.py
  - TODO: Rename to genome_info?
  - Context: 
  - Scope: main
  - Modified: 2025-06-12 14:19:22.879980

