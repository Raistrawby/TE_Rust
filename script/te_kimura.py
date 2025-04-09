#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
import multiprocessing
import pybedtools
import subprocess
import shlex
import shutil
from Bio import AlignIO, SeqIO
from math import log, sqrt
from functools import partial
from time import time
import random
import string
import re

def file_check(repeat_library, in_gff, genome, out_gff, temp_dir):
    # Check that all input files exist.
    for f in [repeat_library, in_gff, genome]:
        if not os.path.exists(f):
            sys.exit(f'Error: File {f} not found.')
    # Index genome if necessary.
    if not os.path.exists(genome + ".fai"):
        print("Indexing genome using samtools faidx...")
        subprocess.run(["samtools", "faidx", genome],
                       stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    # Create temporary directories if they don't exist.
    for d in [temp_dir, os.path.join(temp_dir, "qseqs"), os.path.join(temp_dir, "split_library")]:
        if not os.path.exists(d):
            os.makedirs(d)

def splitter(in_seq, temp_dir):
    """Split the TE consensus FASTA into individual files for each consensus sequence."""
    with open(in_seq, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # Use the part before any '#' as the consensus name and lower-case it.
            consensus_name = record.name.split(sep="#")[0].lower()
            out_file = os.path.join(temp_dir, "split_library", f"{consensus_name}.fasta")
            SeqIO.write(record, out_file, "fasta-2line")

def parse_gff(in_gff):
    """
    Parse the TE annotation GFF file.
    
    This function adapts to your format. It reads all columns and then extracts a new
    column 'consensus' by parsing the metadata. The consensus is taken from the field that starts 
    with 'Target=', from which the first token (before any whitespace) is used.
    """
    gff = pd.read_csv(in_gff, sep="\t", header=None,
                      names=['seqnames', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'metadata'])
    
    def extract_consensus(meta):
        # Split metadata by semicolon.
        fields = meta.split(";")
        for field in fields:
            field = field.strip()
            if field.startswith("Target="):
                # Remove "Target=" and split on whitespace; take the first token.
                return field.replace("Target=", "").split()[0].lower()
        # Fallback to NA if no Target is found.
        return "NA"
    
    gff['consensus'] = gff['metadata'].apply(extract_consensus)
    return gff

def file_name_generator():
    """Generate a random temporary file name."""
    return ''.join(random.sample(string.ascii_letters, 12)) + '.tmp'

def Kimura80(qseq, sseq):
    """
    Calculate the Kimura 2-parameter distance between two aligned sequences.
    """
    transitions = {"AG", "GA", "CT", "TC"}
    transversions = {"AC", "CA", "AT", "TA", "GC", "CG", "GT", "TG"}
    matches = {"AA", "GG", "CC", "TT"}
    m, ts, tv = 0, 0, 0
    for a, b in zip(qseq, sseq):
        pair = a + b
        if pair in matches:
            m += 1
        elif pair in transitions:
            ts += 1
        elif pair in transversions:
            tv += 1
    aln_len = m + ts + tv
    if aln_len == 0:
        return "NA"
    p = ts / aln_len
    q = tv / aln_len
    try:
        kimura_dist = -0.5 * log((1 - 2 * p - q) * (1 - 2 * q)**0.5)
        return str(round(kimura_dist, 4))
    except ValueError:
        return "NA"

# Reorder outer_func so that the first parameter is the GFF chunk.
def outer_func(gff_chunk, genome_path, temp_dir, timeoutSeconds):
    """
    Process a chunk of the GFF: extract the corresponding genome sequence, align it to the TE consensus 
    (using the 'matcher' tool), and compute the Kimura distance.
    """
    tmp_file = os.path.join(temp_dir, file_name_generator())
    failed_file = os.path.join(temp_dir, "failed_" + os.path.basename(tmp_file))
    with open(tmp_file, 'w') as tmp_out:
        # Write header: use all columns except the first (seqnames) plus a Kimura column.
        header = "\t".join(gff_chunk.columns[1:].tolist() + ["Kimura"]) + "\n"
        tmp_out.write(header)
        for idx, row in gff_chunk.iterrows():
            seqname = row['seqnames']
            # Convert to 0-based start for BEDtools.
            start = str(row['start'] - 1)
            end = str(row['end'])
            strand = row['strand']
            consensus = row['consensus']
            bed_str = f"{seqname} {start} {end} . . {strand}"
            query_path = os.path.join(temp_dir, "qseqs", str(idx))
            try:
                a = pybedtools.BedTool(bed_str, from_string=True)
                a.sequence(fi=genome_path, fo=query_path, s=True)
            except Exception:
                with open(failed_file, "a") as ff:
                    ff.write(f"{seqname}:{start}-{end}_{strand}_{consensus}\n")
            if os.path.exists(query_path) and os.path.getsize(query_path) > 0:
                # Look for the consensus FASTA file.
                subject_path = os.path.join(temp_dir, "split_library", f"{consensus}.fasta")
                command = f"matcher {query_path} {subject_path} -outfile {query_path}.matcher -aformat fasta"
                test_command = shlex.split(command)
                try:
                    subprocess.run(
                        test_command,
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.STDOUT,
                        timeout=timeoutSeconds
                    )
                except subprocess.TimeoutExpired:
                    with open(failed_file, "a") as ff:
                        ff.write(f"Timeout expired for {seqname}:{start}-{end}_{strand}_{consensus}\n")
                except Exception as e:
                    with open(failed_file, "a") as ff:
                        ff.write(f"Error processing {seqname}:{start}-{end}_{strand}_{consensus}: {str(e)}\n")
                    
            if not os.path.exists(query_path + ".matcher") or os.path.getsize(query_path + ".matcher") == 0:
                Kdist = "NA"
                if os.path.exists(query_path):
                    os.remove(query_path)
                if os.path.exists(query_path + ".matcher"):
                    os.remove(query_path + ".matcher")
            else:
                aln = list(SeqIO.parse(query_path + ".matcher", 'fasta'))
                ref_seq = str(aln[0].seq).upper()
                gen_seq = str(aln[1].seq).upper()
                if len(ref_seq) == len(gen_seq):
                    Kdist = Kimura80(ref_seq, gen_seq)
                else:
                    Kdist = "NA"
                os.remove(query_path + ".matcher")
                os.remove(query_path)
            # Write the original row data (except the first column) plus the Kimura distance.
            row_values = "\t".join(str(x) for x in row[1:].tolist())
            tmp_out.write(f"{row_values}\t{Kdist}\n")
    return tmp_file

def tmp_out_parser(file_list):
    """Combine temporary output files and update the metadata with the Kimura values."""
    combined = pd.DataFrame()
    for file in file_list:
        temp_df = pd.read_csv(file, sep="\t")
        combined = pd.concat([combined, temp_df], ignore_index=True)
    # Append the Kimura value to the metadata field.
    combined['metadata'] = combined['metadata'] + ";KIMURA80=" + combined['Kimura'].astype(str)
    # Change sorting to use 'source' instead of 'seqnames'
    combined = combined.sort_values(by=['source', 'start']).reset_index(drop=True)
    return combined

def main():
    parser = argparse.ArgumentParser(
        description="Estimate Kimura distance between TE copies and TE consensus sequences.")
    parser.add_argument('-l', '--repeat_library', type=str, required=True,
                        help="Path to TE consensus sequences FASTA")
    parser.add_argument('-i', '--in_gff', type=str, required=True,
                        help="Path to TE annotation GFF file")
    parser.add_argument('-g', '--genome', type=str, required=True,
                        help="Path to genome assembly FASTA")
    parser.add_argument('-o', '--out_gff', type=str, required=True,
                        help="Output GFF file with Kimura distances")
    parser.add_argument('-t', '--temp_dir', type=str, default='tmp/',
                        help='Temporary directory')
    parser.add_argument('-s', '--timeoutSeconds', type=int, default=1200,
                        help="Timeout (in seconds) for matcher tool execution")
    parser.add_argument('-p', '--parallel', type=int, default=multiprocessing.cpu_count(),
                        help="Number of parallel processes to run")
    args = parser.parse_args()

    start_time = time()
    file_check(args.repeat_library, args.in_gff, args.genome, args.out_gff, args.temp_dir)
    
    # Parse GFF file and split the repeat library.
    gff_data = parse_gff(args.in_gff)
    splitter(args.repeat_library, args.temp_dir)
    
    # Divide the GFF data into chunks for parallel processing.
    chunk_size = len(gff_data) // args.parallel if args.parallel > 0 else len(gff_data)
    chunks = [gff_data.iloc[i:i+chunk_size] for i in range(0, len(gff_data), chunk_size)]
    
    # Set temporary directory for pybedtools.
    pybed_temp = os.path.join(args.temp_dir, "pybedtools")
    if not os.path.exists(pybed_temp):
        os.mkdir(pybed_temp)
    pybedtools.set_tempdir(pybed_temp)
    
    # Process in parallel.
    with multiprocessing.Pool(args.parallel) as pool:
        file_list = pool.map(
            partial(outer_func,
                    genome_path=args.genome,
                    temp_dir=args.temp_dir,
                    timeoutSeconds=args.timeoutSeconds),
            chunks
        )
    
    # Combine the output files from all chunks.
    combined_df = tmp_out_parser(file_list)
    combined_df.to_csv(args.out_gff, sep="\t", index=False)
    print(f"Output written to {args.out_gff}")
    run_time = time() - start_time
    print(f"Total run time: {run_time:.2f} seconds")
    
    # Clean up temporary directories.
    shutil.rmtree(os.path.join(args.temp_dir, "split_library"), ignore_errors=True)

if __name__ == "__main__":
    main()