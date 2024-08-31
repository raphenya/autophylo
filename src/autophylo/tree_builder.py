#!/usr/bin/env python

import configparser
import os
import sys
import glob
import json
import logging
import gzip
import bz2
import csv
import math
import re
import subprocess
from statistics import mean
from Bio import SeqIO
from joblib import Parallel, delayed
import warnings


# ==============================================================================
# LOGGING CONFIG
# ==============================================================================
level = logging.WARNING
logger = logging.getLogger(__name__)
logger.setLevel(level)

# detailed log
# formatter = logging.Formatter('%(levelname)s %(asctime)s : \
# (%(filename)s::%(funcName)s::%(lineno)d) : %(message)s')
# basic log
formatter = logging.Formatter('%(levelname)s %(asctime)s : %(message)s')

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

logger.addHandler(stream_handler)

# ==============================================================================
# CONFIG VARIABLES
# ==============================================================================

global MUSCLE, TRIMAL, FastTree, RAxML_PTHREADS, RAxML_HYBRID, GBLOCKS, USEARCH, BLASTDB, TAXIDS


def load_config(config_file_path):
    config_path = os.path.abspath(config_file_path)
    logger.info("config_path %s", config_path)

    if os.path.isfile(config_path) == False:
        logger.error(
            f"Missing 'config' file in the working directory '{os.getcwd()}' ")
        sys.exit(1)

    config = configparser.ConfigParser()
    config.read(config_path)

    global MUSCLE, TRIMAL, FastTree, RAxML_PTHREADS, RAxML_HYBRID, GBLOCKS, USEARCH, BLASTDB, TAXIDS

    MUSCLE = config["ALIGNMENT"]["MUSCLE"]
    TRIMAL = config["ALIGNMENT"]["TRIMAL"]
    FastTree = config["TREE"]["FastTree"]
    RAxML_PTHREADS = config["TREE"]["RAxML_PTHREADS"]
    RAxML_HYBRID = config["TREE"]["RAxML_HYBRID"]
    GBLOCKS = config["CLUSTERING"]["GBLOCKS"]
    USEARCH = config["CLUSTERING"]["USEARCH"]
    BLASTDB = config["DATABASES"]["BLASTDB"]
    TAXIDS = config["DATABASES"]["TAXIDS"]


def is_fasta(input_sequence, extension=""):
    """Checks for valid fasta format."""
    if extension == "":
        with open(input_sequence, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            check_record(fasta)
            return True
    elif extension in ["gz", "bz2"]:
        if extension == "gz":
            with gzip.open(input_sequence, "rt") as handle:
                fasta = SeqIO.parse(handle, "fasta")
                check_record(fasta)
        else:
            with bz2.open(input_sequence, "rt") as handle:
                fasta = SeqIO.parse(handle, "fasta")
                check_record(fasta)
        return True
    else:
        return False


def check_record(fasta, input_type="protein"):
    # check each record in the file
    for record in fasta:
        if any(record.id) == False or any(record.seq) == False:
            return False
        if input_type == "protein":
            return is_protein(record.seq)
        else:
            return False


def is_protein(sequence):
    amino_acids_dict = {
        # common symbols between protein and dna codes
        'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0, 'U': 0,
        # other amino acids
        'R': 0, 'D': 0, 'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0,
        'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0,
        'W': 0, 'Y': 0, 'V': 0, 'X': 0, 'Z': 0, 'J': 0, 'B': 0
    }
    count = 0
    for amino_acid in sequence:
        try:
            amino_acids_dict[amino_acid.upper()] += 1
        except Exception as e:
            sys.stderr.write(f"invalid protein fasta due to: {e}")
            return False

    for a in amino_acids_dict.keys():
        if a not in 'ATGCNU':
            count = count + amino_acids_dict[a]

    if count == 0:
        sys.stderr.write(f"invalid protein fasta: {amino_acids_dict}")
        return False

    # sys.stdout.write("valid protein fasta: {}".format(amino_acids_dict))
    return True


def sampling_ncbi_blastp(input_sequence, output_file, sample_size, threads,
                         skip_taxids, percent_positive_scoring,
                         percent_identity, skip_blast,
                         specific_taxa):
    # blast sequence against RefSeq
    # /usr/bin/get_species_taxids.sh get specific taxon

    # BLASTDBv5
    # https://ftp.ncbi.nlm.nih.gov/blast/db/v5/blastdbv5.pdf
    #

    # https://www.biostars.org/p/107167/
    """
    Note: comment out the following line in get_species_taxids.sh script
    export PATH=/bin:/usr/bin:/am/ncbiapdata/bin:$HOME/edirect:$PATH

    # check database version

    blastdbcmd -info -db nr > info.nr
    blastdbcmd -info -db db > info.db

    Genera:
        cat ../check_taxons/simplefied_taxon_calls.tsv | cut -f7 | sort -h |
        uniq > genera.txt
    Species:
        cat ../check_taxons/simplefied_taxon_calls.tsv | cut -f8 | sort -h |
        uniq > species.txt

    Phyla:
        cat ../check_taxons/simplefied_taxon_calls.tsv | cut -f3 | sort -h |
        uniq > phyla.txt

        - note add Fusobacteria as it was not in the samples
        - failed for https://lpsn.dsmz.de/phylum/desulfobacterota

    Class:
        cat ../check_taxons/simplefied_taxon_calls.tsv | cut -f4 | sort -h |
        uniq > class.txt

    Note: renamed some taxa for searching

    """
    # '''

    output_path, output_filename = os.path.split(
        os.path.abspath(input_sequence))

    TAXIDS = os.path.join(output_path, f"{output_filename}.taxids")

    logger.info("NEW TAXIDS path is: %s", TAXIDS)

    if skip_taxids is False:
        if not os.path.exists(TAXIDS):
            logger.info(
                "download taxids for selected phyla found in the human gut (see phyla.txt file)")
            # make dir
            logger.info("create directory at %s", TAXIDS)
            os.makedirs(TAXIDS)
            get_ids = "get_ids.sh"
            logger.info("specific_taxa: %s", specific_taxa)
            if specific_taxa == "":
                os.system(f"""
				 for taxa in \
				    "Actinobacteriota" \
				    "Bacteroidota" \
				    "Desulfobacterota" \
				    "Firmicutes" \
				    "Proteobacteria" \
				    "Synergistota" \
				    "Verrucomicrobiota" \
				    "Fusobacteria"
				do
				esearch -db taxonomy -query $taxa < /dev/null | efetch -format docsum | xtract -pattern DocumentSummary -element TaxId | while read id
					do
					echo "get_species_taxids.sh -t $id > {TAXIDS}/$id.ids" >> {get_ids}
					done
					sleep 5
				done
				""")
            else:
                os.system(f"""
				 for taxa in "{specific_taxa}"
				do
				esearch -db taxonomy -query $taxa < /dev/null | efetch -format docsum | xtract -pattern DocumentSummary -element TaxId | while read id
					do
					echo "get_species_taxids.sh -t $id > {TAXIDS}/$id.ids" >> {get_ids}
					done
					sleep 5
				done
				""")

            # change file to executable
            logger.info("change file to executable %s", get_ids)
            os.system(f"chmod +x {get_ids}")
            # move script to TAXIDS directory
            logger.info("move script to TAXIDS directory: mv {script_file} {taxids_directory}".format(
                script_file=get_ids, taxids_directory=TAXIDS))
            os.system("mv {script_file} {taxids_directory}/{script_file}".format(
                script_file=get_ids, taxids_directory=TAXIDS))
            # get all taxids? use parallel?
            logger.info("get all taxids: {taxids_directory}/{script_file}".format(
                taxids_directory=TAXIDS, script_file=get_ids))
            os.system(
                "{taxids_directory}/{script_file}".format(taxids_directory=TAXIDS, script_file=get_ids))
            # remove script
            logger.info("remove script: rm {taxids_directory}/{script_file}".format(
                taxids_directory=TAXIDS, script_file=get_ids))
            os.system("rm {taxids_directory}/{script_file}".format(
                taxids_directory=TAXIDS, script_file=get_ids))
        else:
            logger.warning("directory: %s exists", TAXIDS)
            logger.info("skipping taxids downloads ...")

    if skip_blast is False:
        # prepare jobs
        logger.info("TAXIDS : %s", TAXIDS)
        jobs = []
        taxon_files = glob.glob(os.path.join(TAXIDS, "*"))
        logger.info("path: %s", os.path.join(TAXIDS, "*"))
        logger.info("taxon_files : %s", len(taxon_files))
        for taxa_file in taxon_files:
            taxa_file_name = os.path.basename(taxa_file)
            # add filter by identity
            cmd = "blastp -query {input_file} -db {db} -out {input_file}.{taxa_file_name}.tblastnout -seg 'no' -outfmt '7 std positive ppos qcovs' -num_alignments {number_of_hits} -num_threads {num_threads} -taxidlist {taxids}".format(
                input_file=input_sequence,
                number_of_hits=sample_size,
                num_threads=threads,
                db=BLASTDB,
                taxids=taxa_file,
                taxa_file_name=taxa_file_name
            )
            logger.info(cmd)
            jobs.append(cmd)

        n_jobs = len(jobs)
        logger.info("number of jobs: %s", n_jobs)
        completed_jobs = []
        try:
            logger.info("Running jobs in parallel ...")
            completed_jobs = Parallel(n_jobs=n_jobs, verbose=100)(
                delayed(job)(jobs[i]) for i in range(0, n_jobs))
        except Exception as e:
            logger.debug(e)

        # result = {}
        # logger.info("jobs completed: ")
        # print(completed_jobs)

    # '''
    # when running in parallel completed_jobs is a dict of dicts, this for loop merge them all in one
    # for d in completed_jobs:
    # 	result.update(d)
    # print(json.dumps(result,indent=2))

    # loop over the blast results and download sequences using blastdbcmd command, use percent_identity or
    # percent_positive_scoring for filtering
    logger.info("download fasta for each hit ...")
    parse_blast_results(input_sequence, output_file,
                        percent_positive_scoring, percent_identity)


def parse_blast_results(input_sequence, output_file, percent_positive_scoring, percent_identity):
    logger.info("output_file : %s", output_file)
    options = "{ print $2, $3, $14 }"
    options_sequence = '{print ">"$1,$2,$3"\\n"$4}'
    os.system("""cat {input_file}.*.ids.tblastnout | grep -v '#' | awk '{options}' | while read subject_accession percent_identity percent_positive_scoring
	do
		if [ $(echo "$percent_positive_scoring >= {positive_scoring}" | bc) -eq 1 ]; then
			if [ $(echo "$percent_identity >= {identity}" | bc) -eq 1 ]; then
	   			echo "download for $subject_accession $percent_identity $percent_positive_scoring ...";
	   			blastdbcmd -target_only -db {db} -dbtype prot -entry $subject_accession -outfmt '%a|%T|%t|%s' | awk -F '|' '{options_sequence}' >> {output_file}

	   		fi
	   	fi
	done
	""".format(
        input_file=input_sequence,
        output_file=output_file,
        options=options,
        options_sequence=options_sequence,
        db=BLASTDB,
        positive_scoring=percent_positive_scoring,
        identity=percent_identity
    )
    )

    # parallel --no-notice --progress -j+0 'blastdbcmd -target_only -db'


def job(i):
    logger.info("running job: %s", i)
    os.system(i)


def get_average_length(input_sequence):
    sequences = (fasta.seq for fasta in SeqIO.parse(
        open(input_sequence), "fasta"))
    return int(math.ceil(mean(len(seq) for seq in sequences)))


def get_phyla_taxid(path, taxid):
    phyla_taxid = 0
    files = glob.glob(os.path.join(path, "*"))
    for f in files:
        if os.path.isfile(f) and os.path.splitext(os.path.basename(f))[1][1:].strip() in ["ids"]:
            temp = os.path.splitext(os.path.basename(f))
            phyla_taxid = temp[0]
            with open(f, "r") as handle:
                alllines = handle.readlines()
                for lines in alllines:
                    if re.search("^" + lines.strip() + "$", taxid):
                        return phyla_taxid
    return phyla_taxid


def get_phyla_taxid_alt(path):
    phyla_taxid = 0
    results = {}
    files = glob.glob(os.path.join(path, "*"))
    for f in files:
        if os.path.isfile(f) and os.path.splitext(os.path.basename(f))[1][1:].strip() in ["ids"]:
            temp = os.path.splitext(os.path.basename(f))
            phyla_taxid = temp[0]
            with open(f, "r") as handle:
                lines = handle.readlines()
                for line in lines:
                    results[line.strip()] = phyla_taxid
    return results


def clean_fasta(putative_sequence, input_sequence, minimum_length, maximum_length):
    # use this pull species names
    # esearch -db protein -query "WP_115653772.1" | efetch -format docsum | xtract -pattern DocumentSummary -element Organism
    """
    clean fasta: remove duplicates, sort, and fix headers
    """
    # make copy of input file
    c = f"cat {input_sequence} > {input_sequence}.chars.removed"
    logger.debug(c)
    logger.debug(
        f"putative_sequence: {putative_sequence},input_sequence: {input_sequence},minimum_length: {minimum_length},maximum_length: {maximum_length}")
    os.system(c)

    # remove special characters from headers
    os.system(
        f"sed 's/[/()\//g;\]//g' {input_sequence}.chars.removed > {input_sequence}.special_chars.removed")

    # rename duplicates?
    os.system(
        f"seqkit rename --by-name -w 0 < {input_sequence}.special_chars.removed > {input_sequence}.dups_renamed")
    # sort sequences
    os.system(
        f"seqkit sort --by-length --reverse {input_sequence}.dups_renamed > {input_sequence}.sorted")

    # remove duplicates
    os.system(
        f"seqkit rmdup -s -w 0 < {input_sequence}.sorted > {input_sequence}.dups.removed")

    min_len = -1
    max_len = -1

    putative_average = get_average_length(putative_sequence)

    if minimum_length:
        min_len = putative_average - math.floor(putative_average * 0.4)
    if maximum_length:
        max_len = putative_average + math.ceil(putative_average * 0.4)

    logger.info("minimum_length : %s", min_len)
    logger.info("maximum_length : %s", max_len)
    # filter by specified length based on putative sequence
    os.system(
        f"seqkit seq --min-len {min_len} --max-len {max_len} {input_sequence}.dups.removed > {input_sequence}.dups.removed.length")

    # remove sequence with partials in the headers
    info = []
    start = "["
    end = "]"

    with open("{input_file}.dups.removed.clean".format(input_file=input_sequence), 'w') as fout:
        # with open ("{input_file}.dups.removed".format(input_file=input_sequence)) as f1:
        with open("{input_file}.dups.removed.length".format(input_file=input_sequence)) as f1:
            for record in SeqIO.parse(f1, "fasta"):
                r1 = record.description.split(record.id + " ")
                s = record.description
                organism = (s[s.find(start) + len(start):s.rfind(end)])
                r2 = r1[-1].split(" [" + organism + "]")
                description = (r2[0])
                label1 = record.id + "_" + "_".join(organism.split(" "))
                label2 = record.id + "_" + "_".join(description.split(" "))
                label3 = record.id + "_" + \
                    "_".join(organism.split(" ")) + "_" + \
                    "_".join(description.split(" "))
                info.append([record.id, description,
                            organism, label1, label2, label3])
                if "partial" not in record.description:
                    fout.write(">%s\n" % record.id)
                    fout.write("%s\n" % record.seq)

    with open("{input_file}.dups.removed.clean.tsv".format(input_file=input_sequence), 'w') as tab_out:
        writer = csv.writer(tab_out, delimiter='\t', dialect='excel')
        writer.writerow(["id", "description", "organism",
                        "label1", "label2", "label3"])
        for i in info:
            writer.writerow(i)

    output_path, output_filename = os.path.split(
        os.path.abspath(putative_sequence))
    TAXIDS = os.path.join(output_path, "{}.taxids".format(output_filename))
    logger.info("NEW TAXIDS path is: {}".format(TAXIDS))

    res = get_phyla_taxid_alt(TAXIDS)
    # print(json.dumps(res,indent=2))
    # write tabular with phyla taxids
    with open("{input_file}.dups.removed.clean.phyla.tsv".format(input_file=input_sequence), 'w') as tab_out:
        with open("{input_file}.dups.removed.clean.tsv".format(input_file=input_sequence), 'r') as tab_in:
            reader = csv.reader(tab_in, delimiter='\t', quotechar='|')
            writer = csv.writer(tab_out, delimiter='\t', dialect='excel')
            writer.writerow(["id", "description", "organism",
                            "label1", "label2", "label3", "phyla_taxid"])
            for row in reader:
                if row[0] != "id":
                    taxid = row[1].split(" ")[0]
                    try:
                        writer.writerow(
                            [row[0], row[1], row[2], row[3], row[4], row[5], res[taxid]])
                    except Exception as e:
                        print(row, e)


def cluster_samples(samples_fasta, cluster_size, percent_identity):
    """
    ./usearch11.0.667_i86linux32 -cluster_fast cluster_60_70.fasta -id 0.9 -centroids nr.fasta -uc clusters.uc -sort length -clusters clusters/c_
    """
    logger.info("CLUSTER_SIZE: {}".format(cluster_size))
    logger.info("PERCENT_IDENTITY: {}".format(percent_identity))

    # make directory for clusters
    clusters_directory = "{samples_fasta}.clusters.directory".format(
        samples_fasta=samples_fasta)
    if not os.path.exists(clusters_directory):
        os.makedirs(clusters_directory)

    # cluster sequences
    os.system("{usearch} -cluster_fast {samples_fasta} -id {percent_identity} -centroids {samples_fasta}.centroids.fasta\
		-uc {samples_fasta}.clusters.uc -sort length -clusters {clusters_directory}/c_".format(
        usearch=USEARCH, samples_fasta=samples_fasta,
        clusters_directory=clusters_directory,
        percent_identity=percent_identity
    )
    )

    with open("{input_file}.usearch.subsample.fasta".format(input_file=samples_fasta), 'w') as fout:
        # pick X number of sequences from each cluster
        files = glob.glob(os.path.join(clusters_directory, "*"))
        for fsa in files:
            count = 0
            for record in SeqIO.parse(fsa, 'fasta'):
                fout.write(">%s\n" % record.id)
                fout.write("%s\n" % record.seq)
                count = count + 1
                if count == cluster_size:
                    count = 0
                    break


def check_required_programs():
    """ check if the following programs are installed:
     - blast+
     - trimal
     - muscle
     - seqkit
    """
    pass


def add_line(line):
    sys.stderr.write(
        "============================================================\n")
    sys.stderr.write(line + "\n")
    sys.stderr.write(
        "============================================================\n")


def remove_similar_sequences(file, keep):
    all_values = {}
    cutoff = 0.95
    string_filter = "## Identity for most similar pair-wise sequences matrix"
    remove_sequences = []

    with open(file, "r") as handle:
        lines = handle.readlines()
        start = False
        for index, line in enumerate(lines):
            if string_filter in line:
                start = True
            elif start:
                if string_filter not in line:
                    values = line.strip().split("\t")
                    if len(values) == 3:
                        if values[1] in all_values.keys():
                            all_values[values[1]].append(values[0].strip())
                            all_values[values[1]].append(values[2])
                        else:
                            all_values[values[1]] = [
                                values[0].strip(), values[2]]

    # print(json.dumps(all_values,indent=2))
    representative = []
    for key in all_values:
        # keep one sequence from these
        if float(key) >= cutoff:
            if keep in all_values[key]:
                # only use putative
                for i in all_values[key]:
                    if keep != i:
                        remove_sequences.append(i)
            else:
                representative.append(all_values[key][0])
                for i in all_values[key]:
                    if i not in representative:
                        remove_sequences.append(i)

    return remove_sequences


def clean():
    files = glob.glob(os.path.join(os.getcwd(), "*"))
    for f in files:
        if os.path.isfile(f) and "dups" in os.path.basename(f) and os.path.splitext(os.path.basename(f))[1][1:].strip() not in ["tsv", "aln", "tree", "html", "1"]:
            print("Remove file: {}".format(f))


def main(args):
    if args.debug:
        logger.setLevel(10)
        # logger.info(json.dumps(args,indent=2))
    # load config
    add_line("LOAD CONFIG FILE")
    load_config(args.config_file)
    logger.info('%s', MUSCLE)
    logger.info('%s', TRIMAL)
    logger.info('%s', FastTree)
    logger.info('%s', RAxML_PTHREADS)
    logger.info('%s', RAxML_HYBRID)
    logger.info('%s', GBLOCKS)
    logger.info('%s', USEARCH)
    logger.info('%s', BLASTDB)
    logger.info('%s', TAXIDS)

    # check require programs
    add_line("CHECK REQUIRED PROGRAMS")
    check_required_programs()
    add_line("QC INPUTS")
    # qc inputs
    if (is_fasta(args.input_sequence)) is False:
        sys.stderr.write("Invalid FASTA\n")
        logger.error("Invalid FASTA")
        sys.exit(1)

    add_line("SAMPLING: NCBI nr databases")

    if args.skip_sampling is False:
        logger.info("sampling ...")
        sampling_ncbi_blastp(args.input_sequence, args.output_file, args.sample_size, args.threads, args.skip_taxids,
                             args.percent_positive_scoring, args.percent_identity, args.skip_blast, args.specific_taxa)
    else:
        logger.info("skipping sampling ncbi...")
        os.system("echo ''> {}".format(args.output_file))
        parse_blast_results(args.input_sequence, args.output_file,
                            args.percent_positive_scoring, args.percent_identity)

        # qc sampling
    add_line("CLEANING SEQUENCES")
    logger.info("cleaning ...")
    clean_fasta(args.input_sequence, "{input_file}".format(
        input_file=args.output_file), args.minimum_length, args.maximum_length)

    add_line("CLUSTERING")
    if args.skip_usearch is False:
        logger.info("cluster using USEARCH ...")
        # cluster sample sequences and pick representative sets
        cluster_samples("{input_file}.dups.removed.clean".format(
            input_file=args.output_file), args.usearch_cluster_size, args.usearch_percent_identity)

    else:
        logger.info("skipping USEARCH ...")
        os.system("cat {input_file}.dups.removed.clean > {input_file}.dups.removed.clean.usearch.subsample.fasta".format(
            input_file=args.output_file))

        # create a muscle profile before adding putatives

        # add putative
    if args.skip_sampling is True:
        os.system("cat {putative} > {input_file}.dups.removed.clean.usearch.subsample.fasta".format(
            putative=args.input_sequence, input_file=args.output_file))
    else:
        os.system("cat {putative} >> {input_file}.dups.removed.clean.usearch.subsample.fasta".format(
            putative=args.input_sequence, input_file=args.output_file))

        # multiple sequence alignment
        # qc alignment
        # align sequences using muscle
    add_line("ALIGN SEQUENCES")
    logger.info("align using MUSCLE ...")
    # measure entropy for columns
    #  use -align for # sequences < 500 and -super5 for # sequences > 500
    num_sequences = len([1 for line in open("{input_file}.dups.removed.clean.usearch.subsample.fasta".format(
        input_file=args.output_file)) if line.startswith(">")])

    algorithm = "-super5"

    if int(num_sequences) < 500:
        algorithm = "-align"

    logger.info("number of sequences {num_sequences}; align using {algorithm} ...".format(
        num_sequences=num_sequences, algorithm=algorithm))
    muscle_cmd = "{muscle} {algorithm} {input_file}.dups.removed.clean.usearch.subsample.fasta -output {input_file}.dups.removed.clean.usearch.subsample.fasta.muscle.aln -threads {threads} -log {input_file}.dups.removed.clean.usearch.subsample.fasta.muscle.log".format(
        muscle=MUSCLE,
        algorithm=algorithm,
        input_file=args.output_file,
        threads=args.threads
    )
    logger.info(muscle_cmd)
    os.system(muscle_cmd)

    # trim using trimai
    add_line("ALIGNMENT TRIMMING")

    logger.info("trim using TRIMAI (auto) ...")
    os.system("{trimal} -in {input_file}.dups.removed.clean.usearch.subsample.fasta.muscle.aln -out {input_file}.dups.removed.clean.usearch.subsample.fasta.muscle.aln.trimmed.trimal.auto \
		-htmlout {input_file}.dups.removed.clean.usearch.subsample.fasta.muscle.aln.trimmed.trimal.auto.html -fasta -sident -automated1 > {input_file}.pair-wise.txt".format(
        trimal=TRIMAL,
        input_file=args.output_file
    )
    )

    # use pair-wise information from trimal to filter out similar sequences
    sequences_to_remove = []

    if args.remove_similar_sequences is True:
        logger.info("filtering using pair-wise information from trimal ...")
        header_putative = ""
        for record in SeqIO.parse(args.input_sequence, 'fasta'):
            header_putative = record.id
            logger.info("header_putative: ", header_putative)
        sequences_to_remove = remove_similar_sequences(
            f"{args.output_file}.pair-wise.txt", header_putative)
        logger.info("sequences_to_remove: ", sequences_to_remove)

    with open(f"{args.output_file}.dups.removed.clean.usearch.subsample.fasta.muscle.aln.trimmed.trimal.auto.1", 'w') as fout1:
        for record in SeqIO.parse(f"{args.output_file}.dups.removed.clean.usearch.subsample.fasta.muscle.aln.trimmed.trimal.auto", 'fasta'):
            if record.id not in sequences_to_remove:
                fout1.write(">%s\n" % record.id)
                fout1.write("%s\n" % record.seq)

        # trim with Gblocks
        # os.system("{Gblocks} {input_file}.dups.removed.clean.usearch.subsample.fasta.muscle.aln -t=p".format(
        # 	Gblocks=GBLOCKS,
        # 	input_file=args.input_sequence
        # 	)
        # )

        # alignment trimming
        # qc trimming
        # NJ is used for development
    add_line("NEIGHBOUR JOINING (NJ): FastTree")
    # 1st pass NJ phylogenetic tree
    # qc NJ tree
    logger.info("generate NJ using FastTree ...")

    # generate NJ using fasttree
    os.system(f"{FastTree} -log {args.output_file}.dups.removed.clean.usearch.subsample.fasta.muscle.aln.trimmed.trimal.fastree.log \
		-fastest -nj {args.output_file}.dups.removed.clean.usearch.subsample.fasta.muscle.aln.trimmed.trimal.auto.1 \
		> {args.output_file}.dups.removed.clean.usearch.subsample.fasta.muscle.aln.trimmed.trimal.fastree.tree")

    # final phylogenetic tree
    add_line("MAXIMUM LIKELIHOOD (ML) TREE: RAxML")
    if args.skip_maximum_likelihood_tree is False:
        # plot trees using RAxML
        output_path, output_filename = os.path.split(
            os.path.abspath(args.output_file))
        # the directory needs to exist
        os.system("mkdir {output_path}/{output_filename}.results".format(
            output_path=output_path, output_filename=output_filename))
        cmd = "{RAxML} \
			-s {input_file}.dups.removed.clean.usearch.subsample.fasta.muscle.aln.trimmed.trimal.auto.1 -w {output_path}/{output_filename}.results -n {output_filename}.{prefix} \
			-m {model} -T {threads} -f a -p 12345 -x 12345 -N autoMRE".format(
            RAxML=(RAxML_PTHREADS if args.RAxML ==
                   "PTHREADS" else RAxML_HYBRID),
            prefix=(args.specific_taxa if args.specific_taxa != "" else "All"),
            output_path=output_path,
            output_filename=output_filename,
            input_file=args.output_file,
            threads=args.threads,
            model=args.model
        )
        add_line(cmd)
        os.system(cmd)
    else:
        logger.info("skipping ML tree ...")

    if args.clean:
        add_line("CLEAN TEMPORARY FILES")
        clean()

    logger.info("Done.")
