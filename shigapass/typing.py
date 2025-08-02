"""Module to handle blast requests"""

import re
import os
import shutil
import tempfile
import logging
import pandas as pd

from .config import db_config_array, db_alt_config_array, db_mlst_config_array, flex_score_array
from .utils import write_out, read_file, get_val_from_file, check_ascending

logger = logging.getLogger(__name__)

class Typing:
    """Class for Blast arguments"""
    def __init__(self, blast):
        self.blast = blast

    def get_db_array(self, db_marker, rfb_filter=None):
        """Retrieve db array from the config"""
        return [db for db in db_config_array if db.get("db_marker") == db_marker and db.get("rfb") == rfb_filter]

    def determine_ipah(self, outdir, input_filename, input_fasta_file):
        """Determine the ipah type"""
        logger.info("Determing ipaH  for: %s" % input_filename)
        search_parameters = self.get_db_array("ipaH")
        sorted_gene_counts, hits_coverage = self.blast.search(outdir, input_filename, input_fasta_file, search_parameters)
        if sorted_gene_counts:
            ipah = "ipaH+"
            ipah_hits = sorted_gene_counts.get("ipaH")
            ipah_coverage = hits_coverage
        else:
            ipah = "ipaH-"
            ipah_hits = 0
            ipah_coverage = 0
        return ipah, ipah_hits, ipah_coverage

    def determine_alt_rfb(self, outdir, input_filename, input_fasta_file, rfb, rfb_hits, rfb_coverage):
        """Determine alternative rfb types"""
        for db_array in db_alt_config_array:
            if rfb == db_array["rfb"]:
                search_parameters = [db_array]
                sorted_gene_counts, hits_coverage = self.blast.search(outdir, input_filename, input_fasta_file, search_parameters)
                if sorted_gene_counts:
                    logger.info("Rfb has changed from %s to %s" % (rfb, db_array["new_rfb"]))
                    rfb = db_array["new_rfb"]
                    if rfb != "C6":
                        rfb_hits = sorted_gene_counts.get(rfb)
                        rfb_coverage = hits_coverage
                else:
                    rfb = db_array["rfb"]
                    logger.info("Rfb has remained the same: %s" % rfb)
                break
        return rfb, rfb_hits, rfb_coverage

    def determine_rfb(self, outdir, input_filename, input_fasta_file):
        """Determine the rfb type"""
        logger.info("Determining rfb type for: %s" % input_filename)
        search_parameters = self.get_db_array("rfb")
        sorted_gene_counts, hits_coverage = self.blast.search(outdir, input_filename, input_fasta_file, search_parameters)
        rfb, rfb_hits, rfb_coverage, flexserotype, comments = [""] * 5
        if sorted_gene_counts:
            rfb = list(sorted_gene_counts.keys())[0]
            rfb_hits = sorted_gene_counts[rfb]
            rfb_coverage = hits_coverage
            rfb, rfb_hits, rfb_coverage = self.determine_alt_rfb(outdir, input_filename, input_fasta_file, rfb, rfb_hits, rfb_coverage)
            if rfb == "A3b":
                rfb_hits_fpath = os.path.join(outdir, input_filename, "rfb_hits.txt")
                rfb = get_val_from_file(rfb_hits_fpath, 0, 0)[0]
                if rfb == "A3a":
                    rfb = "A3"
                    logger.info("Rfb has changed to %s; hits detected are unique for %s" % (rfb, rfb))
                else:
                    rfb = "A3/A16"
                    logger.info("Hits detected are common with A3 and Aprov97-10607")
            elif rfb == "A3a":
                rfb = "A3"
                logger.info("Rfb has changed to %s; hits detected are unique for %s" % (rfb, rfb))
            elif rfb == "C1":
                search_parameters = self.get_db_array("additionalrfb", "C1")
                sorted_gene_counts, hits_coverage = self.blast.search(outdir, input_filename, input_fasta_file, search_parameters)
                add_rfb_hits_fpath = os.path.join(outdir, input_filename, "additionalrfb_hits.txt")
                try:
                    galF_hit = get_val_from_file(add_rfb_hits_fpath, 1, 0)[0]
                    if galF_hit:
                        rfb = "C1"
                        logger.info("Rfb has remained %s" % rfb)
                    else:
                        rfb = "C20"
                        logger.info("Rfb has changed to %s" % rfb)
                except (IndexError, TypeError):
                    rfb = "C20"
                    logger.info("Rfb has changed to %s" % rfb)
            elif rfb == "B1-5":
                logger.info("Determining phage and plasmid encoded O-antigen modification genes")
                search_parameters = self.get_db_array("POAC")
                sorted_gene_counts, hits_coverage = self.blast.search(outdir, input_filename, input_fasta_file, search_parameters)
                flexserotype = self.get_poac(outdir, input_filename)
            elif not rfb:
                # Search for D serogroup
                search_parameters = [{"db_fasta": "RFB/RFB_serogroup_D_150-mers.fasta", "db_marker": "additionalrfb", "identity": 100, "coverage": 100, "rfb": None}]
                sorted_gene_counts, hits_coverage = self.blast.search(outdir, input_filename, input_fasta_file, search_parameters)
                add_rfb_hits_fpath = os.path.join(outdir, input_filename, "additionalrfb_hitscoverage.txt")
                try:
                    rfb_hits_val = get_val_from_file(add_rfb_hits_fpath, 3, 1)[0]
                    rfb_coverage_val = get_val_from_file(add_rfb_hits_fpath, 3, 3)[0]
                    if rfb_hits_val and int(rfb_hits_val) >= 30:
                        rfb = "D"
                        rfb_hits = rfb_hits_val
                        rfb_coverage = rfb_coverage_val
                        logger.info("Rfb has changed to %s" % rfb)
                    else:
                        rfb = "none"
                        rfb_hits = 0
                        rfb_coverage = 0
                except (IndexError, TypeError, ValueError):
                    rfb = "none"
                    rfb_hits = 0
                    rfb_coverage = 0
            else:
                logger.info("Rfb blast completed")
            comments = self.check_multiple_rfbs(outdir, input_filename, "rfb", rfb)
        else:
            rfb = "none"
        return rfb, rfb_hits, rfb_coverage, flexserotype, comments

    def mlst_search(self, search_dict):
        """Search mlst db for allele matches"""
        mlst_db_fpath = os.path.join(self.blast.db, "MLST/ST_profiles.txt")
        mlst_db_df = pd.read_csv(mlst_db_fpath, delimiter="\t")

        # Convert the DataFrame to a list of dictionaries
        mlst_db_list = mlst_db_df.to_dict(orient="records")
        for st_dict in mlst_db_list:
            if all(st_dict.get(key) == int(value) for key, value in search_dict.items() if key in st_dict):
                return st_dict['ST']
        return "~"

    def determine_mlst(self, outdir, input_filename, input_fasta_file):
        """Determine MLST based on BLAST results"""
        logger.info("Determining MLST type for: %s" % input_filename)
        mlst_loci_dict = {}
        st_loci_missing = False

        for mlst_array in db_mlst_config_array:
            gene = mlst_array["db_marker"]
            self.blast.run_blastn(outdir, input_filename, input_fasta_file,
                             mlst_array["db_fasta"], gene,
                             mlst_array["identity"], mlst_array["coverage"])
            mlst_blastfpath = os.path.join(outdir, input_filename, f"{gene}_allrecords.txt")

            try:
                mlst_blast = read_file(mlst_blastfpath).strip()
                if mlst_blast:
                    # Extract allele number from blast result
                    first_line = mlst_blast.split("\n")[0]
                    allele_info = first_line.split("\t")[1]
                    # Handle different allele naming formats
                    if ":" in allele_info:
                        st = allele_info.split(":")[0].split("-")[1]
                    else:
                        st = allele_info.split("-")[1]
                    mlst_loci_dict[gene] = st
                else:
                    mlst_loci_dict[gene] = "0"
                    logger.info("No MLST hit for: %s" % gene)
                    st_loci_missing = True
            except (IndexError, ValueError) as e:
                mlst_loci_dict[gene] = "0"
                logger.info("Error parsing MLST for %s: %s" % (gene, e))
                st_loci_missing = True

            output = f"{gene}:{mlst_loci_dict[gene]}\n"
            mlst_loci_outfpath = os.path.join(outdir, input_filename, "mlst_alleles.txt")
            write_out(output, mlst_loci_outfpath, "a")

        mlst = self.mlst_search(mlst_loci_dict) if not st_loci_missing else "none"
        mlst_outfpath = os.path.join(outdir, input_filename, "mlst_ST.txt")
        write_out(f"ST{mlst}", mlst_outfpath)
        return f"ST{mlst}"

    def determine_flic(self, outdir, input_filename, input_fasta_file):
        """Determine fliC type"""
        search_parameters = self.get_db_array("flic")
        _, _ = self.blast.search(outdir, input_filename, input_fasta_file, search_parameters)
        flic_records_fpath = os.path.join(outdir, input_filename, "flic_allrecords.txt")

        try:
            flic_hits = get_val_from_file(flic_records_fpath, 11, 1, "\t")
            if flic_hits and flic_hits[0]:
                flic_top_hit = flic_hits[0]
                flic = flic_top_hit.split("_")[1]
            else:
                # Try with reduced stringency
                blast_outfpath = os.path.join(outdir, input_filename, f"flic_blastout.txt")
                filtered_blast_outfpath = os.path.join(outdir, input_filename, f"flic_allrecords.txt")
                self.blast.filter_blast(blast_outfpath, filtered_blast_outfpath, 98, 45)
                flic_hits = get_val_from_file(flic_records_fpath, 11, 1, "\t")
                if flic_hits and flic_hits[0]:
                    flic_top_hit = flic_hits[0]
                    flic = flic_top_hit.split("_")[1]
                else:
                    flic = "none"
        except (IndexError, AttributeError):
            flic = "none"

        return flic

    def determine_crispr(self, outdir, input_filename, input_fasta_file):
        """Determine CRISPR type"""
        logger.info("Determining CRISPR type for: %s" % input_filename)
        search_parameters = self.get_db_array("crispr")
        _, _ = self.blast.search(outdir, input_filename, input_fasta_file, search_parameters)
        crispr_allrecords_fpath = os.path.join(outdir, input_filename, "crispr_allrecords.txt")
        crispr_records_fpath = os.path.join(outdir, input_filename, "crispr_records.txt")

        crispr = "none"
        try:
            ascending = check_ascending(crispr_allrecords_fpath, 8, 9)
            crispr_hits = get_val_from_file(crispr_allrecords_fpath, 11, 1, "\t", ascending, 3)
            if crispr_hits and crispr_hits != [""]:
                crispr = ";".join([crispr_hit.split("_")[1] for crispr_hit in crispr_hits if crispr_hit])
                write_out(crispr, crispr_records_fpath)
        except Exception as e:
            logger.info("Error determining CRISPR for %s: %s" % (input_filename, e))
            crispr = "none"
        return crispr

    def determine_serotype(self, rfb, rfb_hits, rfb_coverage, flexserotype, comments, mlst, flic, crispr):
        """Combine the results to infer the serotype"""
        serotype = f"{rfb},{rfb_hits},{rfb_coverage},{flexserotype},{comments},{mlst},{flic},{crispr}"
        return serotype

    def check_multiple_rfbs(self, outdir, input_filename, db_marker, rfb):
        comments = ""
        hitscov_fpath = os.path.join(outdir, input_filename, f"{db_marker}_hitscoverage.txt")
        content = read_file(hitscov_fpath)
        # Count non-empty lines
        lines = [line for line in content.strip().split('\n') if line.strip()]
        num_lines = len(lines)

        if num_lines >= 3:
            comments = f"More than one rfb is detected: {num_lines}"
            logger.info("Multiple rfbs")
        elif num_lines == 2:
            if rfb == "AprovBEDP02-5104" or rfb == "A16" or  rfb == "A3":
                comments = ""
            else:
                comments = f"More than one rfb is detected: {num_lines}"
        return comments

    def get_poac(self, outdir, input_filename):
        # Read and modify the file content
        poac_hits_fpath = os.path.join(outdir, input_filename, "POAC_hits.txt")
        content = read_file(poac_hits_fpath)

        # Perform the replace using re
        replacements = {
            "gtrX": "32", "gtrII": "4", "gtrIC": "2", "gtrIV": "8", "gtrV": "16",
            "gtrI": "1", "oac1b": "128", "oac": "64", "optII": "256"
        }
        for pattern, repl in replacements.items():
            content = re.sub(pattern, repl, content)
        write_out(content, poac_hits_fpath)
        # Get score
        lines = content.rstrip().split("\n")
        score = sum(int(line.split(";")[0]) for line in lines)
        score_outfpath = os.path.join(outdir, input_filename, "score.txt")
        write_out(f"score={score}\n", score_outfpath)
        logger.info("This is the score: %s" % score)

        # Find the corresponding flexserotype
        flexserotype = "Unknown"
        for flex_array in flex_score_array:
            if score == flex_array["score"]:
                flexserotype = flex_array["flex"]
                break
        logger.info("This is the flexserotype: %s" % flexserotype)

        # Sort, and concatenate phages
        phages = sorted(set(line.split(";")[0] for line in lines))
        phages_str = ";".join(phages)

        # Write the summary to a CSV file
        summary_outfpath = os.path.join(outdir, "ShigaPass_Flex_summary.csv")
        write_out(f"{input_filename};{phages_str};{flexserotype}\n", summary_outfpath, "a")
        return flexserotype

    def run(self, list_file, outdir, keep_files):
        """Main analysis pipeline"""
        # Initialize summary CSV file
        summary_fpath = os.path.join(outdir, "ShigaPass_summary.csv")
        header = "Name;rfb;rfb_hits,(%);MLST;fliC;CRISPR;ipaH;Predicted_Serotype;Predicted_FlexSerotype;Comments\n"
        write_out(header, summary_fpath)

        # Initialize Flex summary if needed
        flex_summary_fpath = os.path.join(outdir, "ShigaPass_Flex_summary.csv")
        if os.path.exists(flex_summary_fpath):
            os.remove(flex_summary_fpath)

        with open(list_file, "r") as f:
            input_files = f.read().splitlines()

            for filepath in input_files:
                logger.info(f"Processing: {filepath}")
                input_filename = os.path.splitext(os.path.basename(filepath))[0]
                tempdir = tempfile.mkdtemp(prefix="ShigaPass_")

                # Copy and parse FASTA file (replace underscores with tildes like bash version)
                input_fasta_file = os.path.join(tempdir, input_filename + "_parsed.fasta")
                with open(filepath, 'r') as fin, open(input_fasta_file, 'w') as fout:
                    for line in fin:
                        if line.startswith('>'):
                            fout.write(line.replace('_', '~'))
                        else:
                            fout.write(line)

                new_outdir = os.path.join(outdir, input_filename)
                if not os.path.exists(new_outdir):
                    os.makedirs(new_outdir)

                # Run analysis pipeline
                ipah, ipah_hits, ipah_coverage = self.determine_ipah(outdir, input_filename, input_fasta_file)

                if ipah == "ipaH-":
                    # Not Shigella/EIEC
                    serotype = "Not Shigella/EIEC"
                    rfb = "ND"
                    mlst = "ND"
                    flic = "ND"
                    crispr = "ND"
                    rfb_hits = "ND"
                    rfb_coverage = "0"
                    flexserotype = ""
                    comments = ""
                else:
                    # Is Shigella/EIEC - run full analysis
                    rfb, rfb_hits, rfb_coverage, flexserotype, comments = self.determine_rfb(outdir, input_filename, input_fasta_file)
                    mlst = self.determine_mlst(outdir, input_filename, input_fasta_file)
                    flic = self.determine_flic(outdir, input_filename, input_fasta_file)
                    crispr = self.determine_crispr(outdir, input_filename, input_fasta_file)

                    # Determine serotype based on results
                    if rfb == "none":
                        serotype = "none"
                    elif rfb == "B1-5" and flexserotype:
                        serotype = f"S. flexneri {flexserotype}"
                    elif rfb in ["A1", "A2", "A3", "A16", "AprovBEDP02-5104"]:
                        serotype = f"S. dysenteriae {rfb}"
                    elif rfb.startswith("B"):
                        serotype = f"S. flexneri {rfb}"  
                    elif rfb.startswith("C"):
                        serotype = f"S. boydii {rfb}"
                    elif rfb == "D":
                        serotype = "S. sonnei"
                    else:
                        serotype = rfb

                # Write summary line
                if not flexserotype:
                    flexserotype = ""

                summary_line = f"{input_filename};{rfb};{rfb_hits},({rfb_coverage});{mlst};{flic};{crispr};{ipah};{serotype};{flexserotype};{comments}\n"
                write_out(summary_line, summary_fpath, "a")

                logger.info(f"Completed analysis for: {input_filename}")

                # Cleanup
                shutil.rmtree(tempdir)
                if not keep_files:
                    shutil.rmtree(new_outdir)

        logger.info(f"Analysis complete. Results written to: {summary_fpath}")
