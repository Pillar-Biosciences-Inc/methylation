#!/usr/bin/env python3

__version__ = "2023.04.11"

import argparse
from datetime import datetime
import fnmatch
import glob
import gzip
import json
import logging
import multiprocessing as mp
import os
import pandas as pd
from pandas.errors import EmptyDataError
import re
import shutil
import signal
import subprocess
import sys
import time

# from subprocess import Popen
# import shlex

__log_path__ = '/tmp/meth.log'

logging.basicConfig(format='[P:%(process)d|%(asctime)s|%(funcName)s:%(lineno)d|%(levelname)s] %(message)s',
                    datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG, encoding='utf-8',
                    filename=__log_path__, filemode="w")
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter(fmt=f"[P:%(process)d|%(asctime)s|%(funcName)s:%(lineno)d|%(levelname)s] %(message)s",
                              datefmt="%m-%d-%Y %I:%M:%S %p")
console.setFormatter(formatter)
logging.getLogger().addHandler(console)


def run_command(command_str: str = None, shell=True):
    if command_str is None:
        raise ValueError("Support.run_command() was called without any command to execute.")
    try:
        logging.debug(f"Attempting to run: {command_str}")
        output = subprocess.check_output(command_str, encoding="utf-8", shell=shell,
                                         stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as ex:
        logging.error(f"Encountered an error executing the command: ")
        logging.error(command_str)
        logging.error(f"Error details:")
        logging.error(f"Exit code={ex.returncode}")
        logging.error(f"Error message={ex.output}")
        logging.error(f"STDERR={ex.stderr}")
        raise subprocess.CalledProcessError
    logging.debug("Command executed without raising any exceptions")
    return output


class Methylation:
    def __init__(self, options):
        # Input arguments
        self.input_dir = None
        self.output_dir = None
        self.threads = None
        self.bed_file = None
        self.genome_path = None
        self.name = None

        # Arguments created during runtime
        self.paired_fastq_files = {}
        self.merged_fastq_files = {}
        self.trimmed_fastq_files = {}
        self.unmapped_fastq_files = {}
        self.bismark_bam_files = {}
        self.bismark_sorted_bam_files = {}
        self.merged_fastq_folder = None
        self.trimmed_fastq_folder = None
        self.bismark_out_folder = None
        self.bismark_sorted_out_folder = None
        self.report_folder = None

        # Output pandas object
        self.bismark_cov_pd = None
        self.bismark_call_pd = None
        self.bismark_bam_pd = None
        self.fq_stats_pd = None

        if self._arg_parse(options):
            self.analysis_stack = ["get_paired_fastq_files",
                                   "merge_paired_fastq",
                                   "parallel_trim_galore",
                                   "parallel_bismark",
                                   "parallel_bismark_report_extractor",
                                   "summary",
                                   "move_files_for_staging"]

    def _arg_parse(self, options):
        logging.info("Parsing arguments...")
        self.input_dir = options.INPUT_DIR
        self.threads = options.THREAD
        self.bed_file = options.BED
        self.genome_path = options.GENOME

        if options.OUTPUT_DIR is None:
            self.output_dir = os.path.join(options.INPUT_DIR, "out_" + datetime.now().strftime("%Y%m%d-%H%M%S"))
        else:
            self.output_dir = os.path.abspath(options.OUTPUT_DIR)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        if len(options.NAME):
            self.name = options.NAME + "_" + datetime.now().strftime("%Y%m%d-%H%M%S")
        else:
            self.name = f"Pillar_Bismark_{__version__}_" + datetime.now().strftime("%Y%m%d-%H%M%S")

        if not os.path.isdir(self.input_dir) or not os.path.exists(self.input_dir):
            # try to split the input in case of galaxy
            logging.error(f"input is not a valid directory {self.input_dir}")
            raise IOError

        if self.bed_file is not None:
            if not os.path.exists(self.bed_file):
                logging.warning(f"The bed file specified {self.bed_file} doesn't exist")

        return True

    def go(self):
        # self.summarize_run()
        while True:
            step = self.analysis_stack[0]
            self.analysis_stack = self.analysis_stack[1:]
            function = getattr(self, step)
            function()
            if len(self.analysis_stack) == 0:
                break

    def get_paired_fastq_files(self):
        """
        search all fastq files in the folder
        get all the pairs by the same name
        separate each sample by lane number again
        produce a nested dictionary structure
        :return:
        :rtype:
        """
        files = self._recursive_glob(self.input_dir, ['*.fastq.gz', '*.fq.gz'])
        logging.info(f"{len(files)} files are found")
        logging.debug(f"Files are:\n{files}")

        if len(files) == 0:
            logging.error(f"No files found!")
            raise IOError

        for mfile in files:
            if not os.path.exists(mfile):
                logging.error(f"{mfile} doesn't exist")
                raise IOError

        # dictionary to hold the info by each sample
        try:
            for file in files:
                pattern = re.compile(r'(.*)_(S\d{1,5})_(L\d{3})_(R\d)_(\d{3})(.*)')
                file_pattern = pattern.search(file)
                file_path = file_pattern.group(0)
                read_number = file_pattern.group(4)
                sample_name = os.path.basename(file_pattern.group(1))

                if sample_name not in self.paired_fastq_files:
                    self.paired_fastq_files[sample_name] = {read_number: [file_path, ]}
                else:
                    if read_number not in self.paired_fastq_files[sample_name]:
                        self.paired_fastq_files[sample_name][read_number] = [file_path, ]
                    else:
                        self.paired_fastq_files[sample_name][read_number].append(file_path)

        except AttributeError:
            logging.error(f"{files} are not in the Illumina fastq format")
            raise IOError

        return True

    @staticmethod
    def _recursive_glob(tree_root, patterns):
        """
        this is a recursive function search all sub-directories of the input treeroot for matched pattern files
        :param tree_root: an input folder
        :type tree_root: file path
        :param patterns: a list of unix expression for fnmatch file selection
        :type patterns: list
        :return: a list of all matched file paths in abs form
        :rtype: list
        """
        if not os.path.isdir(tree_root):
            logging.error(f"input folder {tree_root} is invalid")
            raise IOError
        results = []
        for pattern in patterns:
            for base, dirs, files in os.walk(tree_root):
                good_files = fnmatch.filter(files, pattern)
                results.extend(os.path.join(base, f) for f in good_files)
        if results:
            idx = len(results) - 1
            for filename in reversed(results):
                basename = os.path.basename(filename)
                if len(basename) == 0 or basename[0] == "." or basename[0] == "~" or \
                        basename.startswith("Undetermined"):
                    del results[idx]
                idx -= 1
        return results

    def merge_paired_fastq(self):
        logging.info("+" * 10 + " Concatenating FASTQ files ... " + "+" * 10)
        self.merged_fastq_folder = os.path.join(self.output_dir, "merged_fastqs")
        if not os.path.exists(self.merged_fastq_folder):
            os.makedirs(self.merged_fastq_folder)

        for sample_name, paired_reads in self.paired_fastq_files.items():
            if sample_name not in self.merged_fastq_files:
                self.merged_fastq_files[sample_name] = {"R1": None, "R2": None}
            for read_number, paired_fastq_files in paired_reads.items():
                output_name = os.path.join(self.merged_fastq_folder,
                                           sample_name + "_S1_L001_" + read_number + "_001.fastq.gz")
                cmd = "cat " + " ".join(sorted(paired_fastq_files)) + " > " + output_name
                run_command(cmd)
                self.merged_fastq_files[sample_name][read_number] = output_name

        return True

    @staticmethod
    def _trim_galore(paired_fastq_files, this_output_dir, sample_name):
        logging.info("✂" * 10 + f" Running FASTQ trimming: {sample_name} " + "✂" * 10)
        logging.debug(f"Trimming files: {paired_fastq_files}")
        paired_fastq_file_list = [paired_fastq_files["R1"], paired_fastq_files["R2"]]

        cmd = f"trim_galore --paired -q 30 -o {this_output_dir} " + " ".join(paired_fastq_file_list)
        run_command(cmd)

        return True

    def parallel_trim_galore(self):
        self.trimmed_fastq_folder = os.path.join(self.output_dir, "trimmed_fastq")
        if not os.path.exists(self.trimmed_fastq_folder):
            os.makedirs(self.trimmed_fastq_folder)

        logging.info(f"Trimming {len(self.merged_fastq_files)} samples")
        logging.debug("Processing the following files: " + json.dumps(self.merged_fastq_files, indent=4))
        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        pool = mp.Pool(self.threads, maxtasksperchild=1)
        signal.signal(signal.SIGINT, original_sigint_handler)
        tracking = []
        rets = []
        for sample_name, paired_fastq_files in self.merged_fastq_files.items():
            worker = pool.apply_async(self._trim_galore,
                                      args=(paired_fastq_files, self.trimmed_fastq_folder, sample_name))
            tracking.append(worker)
        try:
            for worker in tracking:
                ret = worker.get()
                rets.append(ret)
        except KeyboardInterrupt:
            pool.terminate()
        finally:
            pool.close()
            pool.join()
            del pool

        # Set the new paths
        for sample_name, paired_fastq_files in self.merged_fastq_files.items():
            r1_file = os.path.join(self.trimmed_fastq_folder, os.path.basename(paired_fastq_files["R1"]))
            r2_file = os.path.join(self.trimmed_fastq_folder, os.path.basename(paired_fastq_files["R2"]))

            # Replace the extension to match trim_galore's behavior
            new_r1_file = re.sub(pattern=r"\.f(ast)?q", repl="_val_1.fq", string=r1_file)
            new_r2_file = re.sub(pattern=r"\.f(ast)?q", repl="_val_2.fq", string=r2_file)
            # Rename files for consistency
            os.rename(src=new_r1_file, dst=r1_file)
            os.rename(src=new_r2_file, dst=r2_file)

            self.trimmed_fastq_files[sample_name] = {"R1": r1_file, "R2": r2_file}

        return True

    @staticmethod
    def _bismark(paired_trimmed_fastq_files, genome, this_output_dir):
        cmd = f"bismark --genome {genome} --non_directional -1 {paired_trimmed_fastq_files['R1']}" + \
              f" -2 {paired_trimmed_fastq_files['R2']} -o {this_output_dir} --temp_dir " + \
              os.path.join(this_output_dir, "bismark_temp --dovetail --un -N 1 --local")
        logging.debug(f"Running {cmd}")
        run_command(cmd)
        return True

    def _sort_index_bam(self):
        cmd = f'find {self.bismark_out_folder} -name \"*.bam\" | grep -v _.bam | sed \'s/.bam//\' | ' \
              f'xargs -I F basename F | xargs -P {self.threads} -I FILE sh -c "samtools sort -o ' \
              f'{self.bismark_sorted_out_folder}/FILE.sorted.bam {self.bismark_out_folder}/FILE.bam; ' \
              f'samtools index {self.bismark_sorted_out_folder}/FILE.sorted.bam"'
        run_command(cmd)
        return True

    def parallel_bismark(self):
        logging.info("䷉" * 10 + " Running Bismark... " + "䷉" * 10)
        # run bismark
        # parallel thread for each sample process.
        self.bismark_out_folder = os.path.join(self.output_dir, "bismark")
        self.bismark_sorted_out_folder = os.path.join(self.output_dir, "bismark_sorted")

        if not os.path.exists(self.bismark_out_folder):
            os.makedirs(self.bismark_out_folder)
        if not os.path.exists(self.bismark_sorted_out_folder):
            os.makedirs(self.bismark_sorted_out_folder)

        logging.debug(f"File list: {self.trimmed_fastq_files}")
        try:
            for sample_name, paired_trimmed_fastq_files in self.trimmed_fastq_files.items():
                logging.info(f"Mapping sample {sample_name}")
                logging.debug(f"Sample {sample_name}, files: {json.dumps(paired_trimmed_fastq_files, indent=4)}")
                self._bismark(paired_trimmed_fastq_files, self.genome_path, self.bismark_out_folder)

                r1_file = os.path.join(self.bismark_out_folder,
                                       os.path.basename(self.trimmed_fastq_files[sample_name]["R1"]))
                r2_file = os.path.join(self.bismark_out_folder,
                                       os.path.basename(self.trimmed_fastq_files[sample_name]["R2"]))
                self.unmapped_fastq_files[sample_name] = {"R1": r1_file + "_unmapped_reads_1.fq.gz",
                                                          "R2": r2_file + "_unmapped_reads_2.fq.gz"}

                bam_file = os.path.join(self.bismark_out_folder,
                                        os.path.basename(self.trimmed_fastq_files[sample_name]["R1"]))
                bam_file = re.sub(pattern="R1_001.*", repl="R1_001_bismark_bt2_pe.bam", string=bam_file)
                self.bismark_bam_files[sample_name] = bam_file

                bam_file = os.path.join(self.bismark_sorted_out_folder,
                                        os.path.basename(self.trimmed_fastq_files[sample_name]["R1"]))
                bam_file = re.sub(pattern="R1_001.*", repl="R1_001_bismark_bt2_pe.sorted.bam", string=bam_file)
                self.bismark_sorted_bam_files[sample_name] = bam_file

        except Exception as ex:
            logging.error(f'Exception raised: {ex}')
            raise RuntimeError

        # sort and index all the bams
        self._sort_index_bam()

    @staticmethod
    def _bismark_methylation_extractor(bam_file, this_output_dir):
        cmd = f"bismark_methylation_extractor -p --bedGraph --report -o {this_output_dir} {bam_file}"
        run_command(cmd)
        return True

    def parallel_bismark_report_extractor(self):
        logging.info("🔍" * 10 + f" Running Bismark extractor on {self.bismark_out_folder}... " + "🔍" * 10)
        self.report_folder = os.path.join(self.output_dir, "report")

        bam_file_list = glob.glob(os.path.join(self.bismark_out_folder, "*.bam"))
        logging.info(f"Found {len(bam_file_list)} BAM files")
        logging.debug(f"File list: {bam_file_list}")

        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
        pool = mp.Pool(self.threads, maxtasksperchild=1)
        signal.signal(signal.SIGINT, original_sigint_handler)
        tracking = []
        rets = []
        for bam_file in bam_file_list:
            worker = pool.apply_async(self._bismark_methylation_extractor, args=(bam_file, self.report_folder))
            tracking.append(worker)
        try:
            for worker in tracking:
                ret = worker.get()
                rets.append(ret)
        except KeyboardInterrupt:
            pool.terminate()
        finally:
            pool.close()
            pool.join()
            del pool

    def _bedtool_coverage(self):
        logging.info(f"Calculating coverage of BAM files in {self.bismark_sorted_out_folder}")
        cmd = f'find {self.bismark_sorted_out_folder} -name \"*.bam\" | grep -v _.bam | sed \'s/.bam//\' | ' \
              f'parallel -j ' + str(self.threads) + ' bedtools coverage -a ' + self.bed_file + ' -b {}.bam' + \
              ' \'>\' \'{}\'_coverage_bp.txt'
        run_command(cmd)

    @staticmethod
    def _get_this_fq_stats(fq_file):
        if fq_file.endswith(".gz"):
            fh = gzip.open(fq_file, "rb")
        else:
            fh = open(fq_file, "r")

        line_counter = 0
        len_counter = 0
        for _ in fh:
            read_seq = fh.readline().rstrip()
            line_counter += 1
            len_counter += len(read_seq)
            fh.readline()
            fh.readline()

        fh.close()

        return [line_counter, len_counter]

    def _get_all_fq_stats(self):
        fq_stats = {"Sample": [],
                    "Pre-QC Reads": [], "Pre-QC Avg Read Len": [],
                    "Post-QC Reads": [], "Post-QC Avg Read Len": [],
                    "Unmapped Reads": [], "Unmapped Avg Read Len": [],
                    "Mapped Reads": [], "Mapping rate (%)": []}

        # stats expected to be identical for R1 and R2
        read_number = "R1"
        for sample in self.merged_fastq_files:
            fq_stats["Sample"].append(sample)

            line_counter, len_counter = self._get_this_fq_stats(self.merged_fastq_files[sample][read_number])
            fq_stats["Pre-QC Reads"].append(line_counter)
            if line_counter == 0:
                line_counter = 1
            fq_stats["Pre-QC Avg Read Len"].append(len_counter / line_counter)

            line_counter, len_counter = self._get_this_fq_stats(self.trimmed_fastq_files[sample][read_number])
            fq_stats["Post-QC Reads"].append(line_counter)
            if line_counter == 0:
                line_counter = 1
            fq_stats["Post-QC Avg Read Len"].append(len_counter / line_counter)

            line_counter, len_counter = self._get_this_fq_stats(self.unmapped_fastq_files[sample][read_number])
            fq_stats["Unmapped Reads"].append(line_counter)
            if line_counter == 0:
                line_counter = 1
            fq_stats["Unmapped Avg Read Len"].append(len_counter / line_counter)

            # Mapped from BAM file
            cmd = f"samtools view -c {self.bismark_sorted_bam_files[sample]}"
            fq_stats["Mapped Reads"].append(int(run_command(cmd).rstrip()) / 2)
            denom = fq_stats["Post-QC Reads"][-1]
            if denom == 0:
                denom = 1
            fq_stats["Mapping rate (%)"].append(fq_stats["Mapped Reads"][-1] * 100 / denom)

        self.fq_stats_pd = pd.DataFrame.from_dict(fq_stats)
        return True

    def _bam_stats(self):
        logging.info(f"Calculating stats of BAM files in {self.bismark_sorted_out_folder}")
        cmd = f'find {self.bismark_sorted_out_folder} -name \"*.bam\" | grep -v _.bam | sed \"s/.bam//\" | ' \
              f'parallel -j {self.threads} samtools flagstat ' \
              '{}.bam \">\" \"{}\"_stats.txt'
        run_command(cmd)

    def _collect_all_cov(self):
        logging.info("Collecting all BED coverage files")
        all_coverage_files = glob.glob(os.path.join(self.bismark_sorted_out_folder, "*.sorted_coverage_bp.txt"))
        if len(all_coverage_files) == 0:
            logging.warning(f"No coverage files found! Can't proceed")
            return False
        cov_df_list = []
        file_ct = 0
        for cov_file in all_coverage_files:
            pattern = re.compile(r'(.*)_(S\d{1,5})_(L\d{3})_(R\d)_(\d{3})(.*)')
            file_pattern = pattern.search(cov_file)
            sample_name = os.path.basename(file_pattern.group(1))
            try:
                df = pd.read_csv(cov_file, sep='\t', header=None)
                df.columns = ["chrom", "start", "end", "amplicon_name", f"{sample_name}_coverage", "roi_length",
                              f"{sample_name}_roi_covered_length", f"{sample_name}_roi_covered_ratio"]
            except EmptyDataError:
                df = pd.DataFrame(columns=["chrom", "start", "end", "amplicon_name", f"{sample_name}_coverage",
                                           "roi_length", f"{sample_name}_roi_covered_length",
                                           f"{sample_name}_roi_covered_ratio"])
            if file_ct == 0:
                df = df[["chrom", "start", "end", "amplicon_name", f"{sample_name}_coverage", ]]
            else:
                df = df[[sample_name + "_coverage"]]
            file_ct += 1
            cov_df_list.append(df)
        self.bismark_cov_pd = pd.concat(cov_df_list, axis=1)
        return True

    def _collect_bam_stats(self):
        logging.info("Collecting BAM stats files")
        all_stats_files = glob.glob(os.path.join(self.bismark_sorted_out_folder, "*_stats.txt"))
        stats_results = []
        for stat_file in all_stats_files:
            pattern = re.compile(r'(.*)_(S\d{1,5})_(L\d{3})_(R\d)_(\d{3})(.*)')
            file_pattern = pattern.search(stat_file)
            sample_name = os.path.basename(file_pattern.group(1))
            pname = f"{sample_name}_passing_reads"
            fname = f"{sample_name}_failing_reads"
            dict_to_store = {f"Metric": [], pname: [], fname: []}
            with open(stat_file, "r") as f:
                header = re.sub(pattern=r"\sin total.*", repl="", string=f.readline().rstrip())
                pass_stat, fail_stat = re.split(pattern=r"\s*\+\s*", string=header, maxsplit=1)
                dict_to_store["Metric"].append("BAM QC reads")
                dict_to_store[pname].append(pass_stat)
                dict_to_store[fname].append(fail_stat)
                for line in f:
                    pass_stat, fail_stat = re.split(pattern=r"\s*\+\s*", string=line.rstrip(), maxsplit=1)
                    fail_stat, metric = re.split(pattern=r"\s+", string=fail_stat, maxsplit=1)
                    dict_to_store["Metric"].append(metric)
                    dict_to_store[pname].append(pass_stat)
                    dict_to_store[fname].append(fail_stat)

            # df = pd.read_csv(stat_file, sep='\t', header=None)
            # df.columns = [sample_name]
            df = pd.DataFrame.from_dict(dict_to_store)
            df = df.set_index('Metric')
            stats_results.append(df)

        self.bismark_bam_pd = pd.concat(stats_results, axis=1)

    def _collect_bismark_call(self):
        logging.info("Collect Bismark call files")
        all_call_files = glob.glob(os.path.join(self.report_folder, "*.bismark.cov.gz"))
        call_results = []
        for call_file in all_call_files:
            pattern = re.compile(r'(.*)_(S\d{1,5})_(L\d{3})_(R\d)_(\d{3})(.*)')
            file_pattern = pattern.search(call_file)
            sample_name = os.path.basename(file_pattern.group(1))
            df = pd.read_csv(call_file, compression="gzip", sep='\t', header=None)
            df.columns = ["chrom", "start", "end", "methylation", "methylated_reads", "unmethylated_reads"]
            df["Sample_ID"] = sample_name
            call_results.append(df)

        self.bismark_call_pd = pd.concat(call_results)
        return True

    def summary(self):
        # output all dfs to an Excel summary
        output_path = os.path.join(self.output_dir, self.name + '.xlsx')
        writer = pd.ExcelWriter(output_path, engine='xlsxwriter')

        # get the coverage info
        if self.bed_file is not None:
            self._bedtool_coverage()
            # aggregate coverage data and methylation call data
            # bed files first
            if self._collect_all_cov():
                self.bismark_cov_pd.to_excel(writer, sheet_name="coverage", index=False)
        else:
            logging.warning("no bed file provided for coverage")

        # get FASTQ stats:
        logging.info(f"Calculating FASTQ stats")
        self._get_all_fq_stats()
        logging.info("Writing FASTQ stats to the excel file")
        self.fq_stats_pd.to_excel(writer, sheet_name="fastq_stats", index=False)

        # get general BAM stats:
        logging.info(f"Calculating BAM stats")
        self._bam_stats()
        logging.info(f"Collating all the stats")
        self._collect_bam_stats()
        logging.info("Writing BAM stats to the excel file")
        self.bismark_bam_pd.to_excel(writer, sheet_name="bam_stats", index=True)

        logging.info(f"Aggregating all the methylation information")
        self._collect_bismark_call()
        logging.info("Writing call stats to the excel file")
        self.bismark_call_pd.to_excel(writer, sheet_name="variants", index=False)

        writer.close()

    def move_files_for_staging(self):
        logging.info(f"Moving sorted BAM/BAI files to output directory {self.output_dir}")
        for sample in self.bismark_sorted_bam_files:
            bam_file = self.bismark_sorted_bam_files[sample]
            logging.debug(f"Moving {bam_file} to {self.output_dir}")
            shutil.move(bam_file, os.path.join(self.output_dir, os.path.basename(bam_file)))

            bai_file = bam_file + ".bai"
            logging.debug(f"Moving {bai_file} to {self.output_dir}")
            shutil.move(bai_file, os.path.join(self.output_dir, os.path.basename(bai_file)))

            r1_un = self.unmapped_fastq_files[sample]["R1"]
            r2_un = self.unmapped_fastq_files[sample]["R2"]
            logging.debug(f"Moving {r1_un} and {r2_un} to {self.output_dir}")
            shutil.move(r1_un, os.path.join(self.output_dir, os.path.basename(r1_un)))
            shutil.move(r2_un, os.path.join(self.output_dir, os.path.basename(r2_un)))


if __name__ == "__main__":

    logging.info(f"Methylation v.{__version__}")
    start_time = time.time()

    parser = argparse.ArgumentParser(prog='Bismark')
    parser.add_argument('-i', dest='INPUT_DIR', required=True,
                        help='REQUIRED! input directory or files containing Fastq files separated by comma.')
    parser.add_argument('-o', dest='OUTPUT_DIR', default=None,
                        help='REQUIRED! output directory, Galaxy instance defaults to service definition')
    parser.add_argument("-t", dest='THREAD', default=1,
                        type=int, help="Number of threads used for parallelism, defaults to 1")
    parser.add_argument("-b", dest='BED', default=None,
                        help="Bed file input")
    parser.add_argument("-g", dest='GENOME', default="/pillar/genomes/hg19/",
                        help="Genome location")
    parser.add_argument("-n", dest='NAME', default="", help='run_name')
    parser.add_argument("-v", action="version", version='Pillar %(prog)s {version}'.format(version=__version__))
    args, unknown_arguments = parser.parse_known_args()

    if len(unknown_arguments) > 0:
        print("User specified unknown arguments: " + " ".join([str(x) for x in unknown_arguments]))
        print("Please see the correct usage below:")
        parser.print_help()
        raise ValueError

    logging.debug(f"Input arguments: {sys.argv}")

    try:
        meth_head = Methylation(args)
        meth_head.go()
        end_time = time.time()
        logging.info("🕙" * 10 + f" Total runtime {round(end_time - start_time, 2)}s " + "🕞" * 10)
    except Exception as e:
        # copy the log folder to the output
        logging.error(e)
        shutil.copy(__log_path__, os.path.join(meth_head.output_dir, "log-" + meth_head.name + ".txt"))
        sys.exit(1)
    else:
        # copy the log folder to the output
        shutil.copy(__log_path__, os.path.join(meth_head.output_dir, "log-" + meth_head.name + ".txt"))
        sys.exit(0)

