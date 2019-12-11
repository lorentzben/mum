
# script to run mothur following Kelly lab protocol or user-submitted paramters
# Ben Lorentz 11.19.19
# Input : Directory of fastq files or tar.gz files, design.txt
# Output: FILENAME.axes.csv, FILENAME.final.subs.opti_mcc.groups.summary, FILENAME.loadings, design.txt, version.DATE
# Required: mothur-py, mothur (v1.39.5), python 3, os.sys, readTableNumpy.py
import os
import glob
import csv
import readTableNumpy
import mothur_batch
import argparse
import datetime
import logging
from Bio import SeqIO


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# Logging handler which catches EVERYTHING
file_logger = logging.FileHandler('mum.log')
file_logger.setLevel(logging.DEBUG)
# Logging handler which logs less
console_logger = logging.StreamHandler()
console_logger.setLevel(logging.ERROR)

# Formats the logs so they are pretty
logFormatter = '%(asctime)s- %(name)s - %(lineno)s - %(levelname)s - %(message)s'
formatter = logging.Formatter(logFormatter)
file_logger.setFormatter(formatter)
console_logger.setFormatter(formatter)

# adds handlers to logger
logger.addHandler(file_logger)
logger.addHandler(console_logger)


try:
    logger.info("checking if mothur_py is installed")
    from mothur_py import Mothur
except ImportError:
    logging.critical("install mothur-py using pip")

# Checks to make sure mothur is installed


def setup():
    try:
        logger.info("Making a mothur instace exists")
        m = Mothur()
    except:
        logging.critical("install mothur-py using pip")
    design_oligos_ex()
    check_rdp()
    check_silva()
    samples_ex()
    check_sample_qual()


# Checks if a design file and primers exists
def design_oligos_ex():
    if os.path.exists("design.txt"):
        logger.info("design.txt exists")temp_qual_score = 0
            "Make sure a design file named design.txt is in this directory")
        exit(1)
    if os.path.exists("oligos.txt"):
        logger.info("oligos.txt exists")
    else:
        logger.critical(
            "Make sure the primers used in PCR are in oligos.txt in this directory")
        exit(1)

# Checks to see if .fastq files exist, if yes if they are div by 2, if no check fastq.tar.gz exist

def check_sample_qual():
    where_am_i = os.getcwd()
    my_fastqs = []
    #Looks through directory and adds a location for all .fastq files found 
    for file in os.listdir(where_am_i):
        if file.endswith(".fastq"):
            my_fastqs = os.path.join(where_am_i , file)
    #will hold tuples of seq name and average qual score
    average_qual_scores = []
    # creates a list of quality scores and appends it to a list 
    for read in my_fastqs:
        for record in SeqIO.parse(read, "fastq"):
            quals = record.format("phread_qual")
            average_qual_scores.append((read,mean(quals))
    return average_qual_scores
            

def samples_ex():
    logger.info("Checking to see if samples exist")
    fastq = glob.glob('./*.fastq')
    if not fastq:

        fastq_tar_gz = glob.glob('./*.fastq.gz')
        if not fastq_tar_gz:
            logging.critical(
                "Please check to make sure fastq or fastq.gz files exist")
        else:
            logger.info("Unzipping files")
            os.system("gunzip *.fastq.gz")
    else:
        logger.info(str(len(fastq)) + " Fastq files detected " +
                    str((len(fastq)/2)) + " Samples expected ")
        if (len(fastq) % 2 != 0):
            logger.error("if using paired reads check for missing reads")


# Check to see if trainset16_022016.pds.fasta, and trainset16_022016.pds.tax exists or downloads them


def check_rdp():
    if (os.path.exists("trainset16_022016.pds.fasta") and os.path.exists("trainset16_022016.pds.tax")):
        logger.info("rdp exists")
    else:
        logger.info("Will collect RDP files")
        os.system(
            "wget https://www.mothur.org/w/images/c/c3/Trainset16_022016.pds.tgz")
        os.system("tar -xzf Trainset16_022016.pds.tgz")
        os.chdir("trainset16_022016.pds")
        os.system("cp trainset16_022016.pds.fasta ..")
        os.system("cp trainset16_022016.pds.tax ..")
        os.chdir("..")
        os.system("rm -rf Trainset16_022016.pds.tgz")
        logger.info("rdp files collected")
# Checks to see if silva.bacteria.fasta exists or downloads it


def check_silva():
    if os.path.exists("silva.bacteria.fasta"):
        logger.info("silva exists")
    else:
        logger.info("Will collect silva files")
        os.system("wget https://www.mothur.org/w/images/9/98/Silva.bacteria.zip")
        os.system("unzip Silva.bacteria.zip")
        os.chdir("silva.bacteria")
        os.system("cp silva.bacteria.fasta ..")
        os.chdir("..")
        os.system("rm -rf silva.bacteria")
        os.system("rm -rf Silva.bacteria.zip")
        os.system("rm -rf __MACOSX")
        logger.info("files collected")


# TODO Check that the mothur run was sucessful
def mothurCompleted():
    # TODO Check if the three files exist in the end
    print("ok")

    # TODO move them into a dir called output for r analysis


def main(arg):
    logger.info("running setup checks")
    setup()
    mothur_batch.set_up_logger(arg.quiet)
    readTableNumpy.set_up_logger(arg.quiet)
    logger.info("calling mothur batch")
    mothur_batch.mothur_batch(project_name=arg.job_name, standard=not arg.custom, max_len=arg.max_len,
                              pre_clust_val=arg.pre_clust, design="design.txt", sub_samp_size=arg.sub_sample)


if __name__ == "__main__":
    # Build Argument Parser in order to facilitate ease of use for user
    parser = argparse.ArgumentParser(
        description="Perform Analysis of 16s Microbial Data using Mothur")
    parser.add_argument('-n', action='store', required=True,
                        help="name for analysis, will be filename for resultant files", dest='job_name')
    parser.add_argument('-c', action='store_true', default=False, required=False,
                        help="flag to use custom parameters, ignore for kelly lab parameters", dest='custom')
    parser.add_argument('-l', action='store', default=300, type=int, required=False,
                        help='determines the longest value permitted for sequences', dest='max_len')
    parser.add_argument('-p', action='store', default=2, type=int, required=False,
                        help='pre cluster value, higher is more stringent', dest='pre_clust')
    parser.add_argument('-s', action='store', type=int, required=False,
                        help='sub sample value, only applicable for custom runs', dest='sub_sample')
    parser.add_argument('-q', '--quiet', action='store_true', default=False,
                        help="Reduces the amount of text printed to terminal, check logfiles more often", dest='quiet')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s 1.0')

    args = parser.parse_args()
    # print(args)
    main(args)
