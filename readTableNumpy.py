# readTableNumpy.py a script to make a human readable table of decreasing sequeces over the course of mothur analysis
# Ben Lorentz, created on 12.4.2019
import numpy
import csv
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
# Logging handler which catches EVERYTHING
file_logger = logging.FileHandler('mum.log')
file_logger.setLevel(logging.DEBUG)
# Logging handler which logs less
console_logger = logging.StreamHandler()


def set_up_logger(quiet):
    if quiet:
        console_logger.setLevel(logging.WARNING)
    else:
        console_logger.setLevel(logging.INFO)


# Formats the logs so they are pretty
logFormatter = '%(asctime)s- %(name)s - %(lineno)s - %(levelname)s - %(message)s'
formatter = logging.Formatter(logFormatter)
file_logger.setFormatter(formatter)
console_logger.setFormatter(formatter)

# adds handlers to logger
logger.addHandler(file_logger)
logger.addHandler(console_logger)


def create_seq_table(project_name):
    # creates reads_table and populates header
    logger.info("created the empty table and added headers")
    reads_table = []
    reads_table.append(["Step", "Ambigs", "No. Seq", "% Total Seq"])

    # prepares the sum of seqs
    logger.info("finds count of the first set of seqs")
    max_seqs = 0
    assembly_3 = numpy.loadtxt(
        project_name+".trim.contigs.summary", usecols=(4, 6), skiprows=1)
    assembly_3_vals = numpy.sum(assembly_3, axis=0)

    # populates the max_seqs variable for calculating the percent remaining seqs
    max_seqs = assembly_3_vals[1]
    # appends the data from assembly step to table
    logger.info("adds info about step 3 to table")
    reads_table.append(["Assembly_3", numpy.around(assembly_3_vals[0], decimals=0),
                        numpy.around(assembly_3_vals[1], decimals=0), numpy.around((assembly_3_vals[1]/max_seqs)*100, decimals=2)])

    # calculates the sum of seqs at ambiguous step 4 and then appends
    ambig_5 = numpy.loadtxt(
        project_name+".trim.contigs.good.summary", usecols=(4, 6), skiprows=1)
    ambig_5_vals = numpy.sum(ambig_5, axis=0)
    logger.info("adds the info from step 5 to the table")
    reads_table.append(["Ambigous_5", numpy.around(ambig_5_vals[0], decimals=0),
                        numpy.around(ambig_5_vals[1], decimals=0), numpy.around((ambig_5_vals[1]/max_seqs)*100, decimals=2)])

    # calculates the sum of seqs at Unique step 8 and then appends
    unique_8 = numpy.loadtxt(
        project_name+".trim.contigs.good.unique.summary", usecols=(4, 6), skiprows=1)
    unique_8_vals = numpy.sum(unique_8, axis=0)
    logger.info("adds info from step 8 to table")
    reads_table.append(["Unique_8", numpy.around(unique_8_vals[0], decimals=0),
                        numpy.around(unique_8_vals[1], decimals=0), numpy.around((unique_8_vals[1]/max_seqs)*100, decimals=2)])

    # calculates the sum of seqs at aligned step 13 and then appends
    aligned_13 = numpy.loadtxt(
        project_name+".trim.contigs.good.unique.summary", usecols=(4, 6), skiprows=1)
    aligned_13_vals = numpy.sum(aligned_13, axis=0)
    logger.info("adds info from step 13 to table")
    reads_table.append(["Aligned_13", numpy.around(aligned_13_vals[0], decimals=0),
                        numpy.around(aligned_13_vals[1], decimals=0), numpy.around((aligned_13_vals[1]/max_seqs)*100, decimals=2)])

    # calculates the sum of seqs at chimera step 20 and then appends
    chimera_20 = numpy.loadtxt(
        project_name+".trim.contigs.good.unique.good.filter.unique.precluster.pick.summary", usecols=(4, 6), skiprows=1)
    chimera_20_vals = numpy.sum(chimera_20, axis=0)
    logger.info("adds info from step 20 to table")
    reads_table.append(["Chimera_20", numpy.around(chimera_20_vals[0], decimals=0),
                        numpy.around(chimera_20_vals[1], decimals=0), numpy.around((chimera_20_vals[1]/max_seqs)*100, decimals=2)])

    return reads_table


def print_seq_table(seq_table):
    for line in seq_table:
        print(*line,)


def write_seq_table_to_file(seq_table, project_name):
    with open("seqTable"+project_name+".txt", "a+") as seq_table_out:
        for line in seq_table:
            print(*line, file=seq_table_out)

    logger.info("Table written to: seqTable"+project_name+".txt")

# Test case
#temp = create_seq_table("test")
# print_seq_table(temp)
# write_seq_table_to_file(temp,"test")
