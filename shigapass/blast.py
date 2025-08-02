"""Module to handle blast requests"""

import os
import subprocess
import logging

from .utils import dict_to_str, write_out

logger = logging.getLogger(__name__)

class Blast:
    """Class for Blast arguments"""
    def __init__(self, db, threads):
        self.db = db
        self.threads = threads

    def get_hits(self, outdir, input_filename, db_marker):
        gene_counts = {}
        filtered_records_fpath = os.path.join(outdir, input_filename, f"{db_marker}_allrecords.txt")
        with open(filtered_records_fpath, "r") as fin:
            for line in fin:
                parts = line.strip().split("\t")
                gene = parts[1].split("_")[0]
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
        sorted_gene_counts = dict(sorted(gene_counts.items(), key=lambda item: item[1], reverse=True))
        outfpath = os.path.join(outdir, input_filename, f"{db_marker}_hits.txt")
        content = dict_to_str(sorted_gene_counts)
        write_out(content, outfpath)
        logger.info("The gene hits output was written to: %s" % outfpath)
        return sorted_gene_counts

    def coverage_hits(self, outdir, input_filename, db_marker):
        hits_fpath = os.path.join(outdir, input_filename, f"{db_marker}_hits.txt")
        rfb_hits_count_fpath = os.path.join(self.db, "RFB_hits_count.csv")
        hitscov_fpath = os.path.join(outdir, input_filename, f"{db_marker}_hitscoverage.txt")

        # Read reference hit counts
        ref_counts = {}
        with open(rfb_hits_count_fpath, "r") as fin2:
            for line2 in fin2:
                parts2 = line2.rstrip().split(";")
                if len(parts2) >= 2:
                    ref_counts[parts2[0]] = parts2[1]

        hits_coverage = 0
        with open(hits_fpath, "r") as fin1, open(hitscov_fpath, "w") as fout:
            for line in fin1:
                parts1 = line.rstrip().split(";")
                if len(parts1) >= 2 and parts1[0] in ref_counts:
                    hits_coverage = float(parts1[1]) / float(ref_counts[parts1[0]]) * 100
                    fout.write(f"{line.rstrip()};{ref_counts[parts1[0]]};{hits_coverage:.2f}\n")
                    break

        logger.info("The hits coverage was written to: %s" % hitscov_fpath)
        return hits_coverage

    def run_blastn(self, outdir, input_filename, input_fasta_file, db_fasta, db_marker, identity, coverage):
        #blast_cmd = f"blastn -db {self.db}/{db_fasta} -query {input_fasta_file} -out {os.path.join(outdir, input_filename, f"{db_marker}_blastout.txt")} -num_threads {self.threads} -num_alignments 10000 -outfmt 6 -word_size 11 -dust no"
        blast_cmd = ["blastn", "-query", input_fasta_file, "-db", os.path.join(self.db, db_fasta),
                     "-out", os.path.join(outdir, input_filename, f"{db_marker}_blastout.txt"),
                     "-num_threads", str(self.threads), "-num_alignments", "10000", "-outfmt", "6",
                     "-word_size", "11", "-dust", "no"]
        subprocess.run(blast_cmd, check=True)
        blast_outfpath = os.path.join(outdir, input_filename, f"{db_marker}_blastout.txt")
        filtered_blast_outfpath = os.path.join(outdir, input_filename, f"{db_marker}_allrecords.txt")
        self.filter_blast(blast_outfpath, filtered_blast_outfpath, identity, coverage)

    def filter_blast(self, blast_outfpath, filtered_blast_outfpath, identity, coverage):
        with open(blast_outfpath, "r") as blastfile, \
            open(filtered_blast_outfpath, "w") as fout:
            for line in blastfile:
                parts = line.rstrip().split("\t")
                length = float(parts[3])
                try:
                    q_len = float(parts[1].split("_")[-1])
                except ValueError:
                    q_len = float(parts[1].split(":")[-1])
                if float(parts[2]) >= identity and (length / q_len) * 100 >= coverage:
                    fout.write(line)
        logger.info("The blast filtered output was written to: %s" % filtered_blast_outfpath)

    def search(self, outdir, input_filename, input_fasta_file, db_array):
        search_parameters = db_array[0]
        self.run_blastn(outdir, input_filename, input_fasta_file, search_parameters["db_fasta"],
            search_parameters["db_marker"], search_parameters["identity"],
            search_parameters["coverage"])
        sorted_gene_counts = self.get_hits(outdir, input_filename, search_parameters["db_marker"])
        hits_coverage = self.coverage_hits(outdir, input_filename, search_parameters["db_marker"])
        return sorted_gene_counts, hits_coverage

    def make_blast_db(self):
        for root, dirs, files in os.walk(self.db):
            for filename in files:
                if filename.endswith(".fasta"):
                    makeblastdb_cmd = ["makeblastdb", "-dbtype", "nucl", "-in", os.path.join(root, filename)]
                    subprocess.run(makeblastdb_cmd, check=True)
