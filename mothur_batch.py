# mothur_batch.py a script for running the microbial community analysis tool Mothur with python following the Kelly Lab Protocol uses the mothur_py package to be called through python
# Ben Lorentz blorentz@luc.edu created on 12.4.2019

import readTableNumpy
import os
import csv
import logging


def set_up_logger(quiet):
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    # Logging handler which catches EVERYTHING
    file_logger = logging.FileHandler('mum.log')
    file_logger.setLevel(logging.DEBUG)
    # Logging handler which logs less
    console_logger = logging.StreamHandler()
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


try:
    from mothur_py import Mothur
except ImportError:
    logger.critical("install mothur-py using pip")


def mothur_batch(project_name, standard, max_len, pre_clust_val, design, sub_samp_size):
    # This one uses the trainset16_022016.pds.fasta
    # Uses the silva.bacteria.fasta

    from mothur_py import Mothur
    m = Mothur()
    # figures out what files are in the folder
    m.make.file(inputdir=".", type="fastq", prefix=project_name)
    logger.info("Made File: (1/38)")
    # figures out how many processors exist and uses 40% of them
    os.system("nproc > processors.txt")
    logger.debug("nproc > processors.txt")
    with open("processors.txt") as temp:
        proc = temp.read()
    # print(float(proc))
    # puts forward and backwards reads together into one
    nproc = int(float(proc)*.4)
    logger.debug(str(nproc) + " processors")
    # print(nproc)
    os.system("rm -rf processors.txt")
    logger.debug("rm -rf processors.txt")
    m.make.contigs(file=project_name+".files", processors=nproc)
    logger.info("Made Contigs: (2/38)")
    # read len should be around 292 or so
    m.summary.seqs(fasta=project_name+".trim.contigs.fasta")
    logger.info("Summarized Seqs: (3/38)")
    # checks if a Kelly lab run or if has different parameters
    if standard == True:
        logger.info("standard Run")
        max_length = 300
    else:
        logger.info("custom run")
        max_length = max_len
    # removes ambiguous seqs or seqs that are too long
    m.screen.seqs(fasta=project_name+".trim.contigs.fasta",
                  group=project_name+".contigs.groups", maxambig=0, maxlength=max_length)
    logger.info("Screened Seqs: (4/38)")
    # Seqs should be shorter than 300 and have no ambigs
    m.summary.seqs(fasta=project_name+".trim.contigs.good.fasta")
    logger.info("Summarized Seqs: (5/38)")
    # removes duplicates
    m.unique.seqs(fasta=project_name+".trim.contigs.good.fasta")
    logger.info("Removed Duplicates: (6/38)")
    # rows are names of unique seqs columns are names of groups num of times that unique seq shows up
    m.count.seqs(name=project_name+".trim.contigs.good.names",
                 group=project_name+".contigs.good.groups")
    logger.info("Counted Seqs: (7/38)")
    # number of unique seqs and total number of seqs
    m.summary.seqs(count=project_name+".trim.contigs.good.count_table")
    logger.info("Summarized Seqs: (8/38)")
    # aligns the silva file to oligos
    m.pcr.seqs(fasta="silva.bacteria.fasta", oligos="oligos.txt")
    logger.info("Trimmed Seqs Based on Oligos: (9/38)")
    # gives the file a new name
    m.rename.file(input="silva.bacteria.pcr.fasta", new="silva.v4.fasta")
    logger.info("Renamed Trimmed Seqs: (10/38)")
    # 253 bases and gives start and end values
    m.summary.seqs(fasta="silva.v4.fasta")
    logger.info("Summarized Seqs: (11/38)")
    # align my seqs to the database sequences
    m.align.seqs(fasta=project_name+".trim.contigs.good.unique.fasta",
                 reference="silva.v4.fasta")
    logger.info("Aligns Samples to Generated File: (12/38)")
    # See what the alignment looks like should be about 253
    m.summary.seqs(fasta=project_name+".trim.contigs.good.unique.align",
                   count=project_name+".trim.contigs.good.count_table")
    logger.info("Summarized Seqs: (13/38)")
    # Pulled from from project_name.trim.contigs.good.unique.summary
    # start=13862
    # end=23444
    start = 0
    end = 0
    count = 0
    with open(project_name+".trim.contigs.good.unique.summary") as summary:

        contents = csv.reader(summary, delimiter='\t')

        for line in contents:
            try:

                new_1 = int(float(line[1]))
                new_2 = int(float(line[2]))
                start = start+new_1
                end = end+new_2
                count = count + 1

            except ValueError:
                pass
    f_start = int(start/count)
    f_end = int(end/count)
    logger.info("start: "+str(f_start)+" end: "+str(f_end))

    # cuts seqs based on the numbers above
    m.screen.seqs(fasta=project_name+".trim.contigs.good.unique.align", count=project_name+".trim.contigs.good.count_table",
                  summary=project_name+".trim.contigs.good.unique.summary", start=f_start, end=f_end, maxhomop=8)
    logger.info("Removed Seqs Outside of Alignment Range: (14/38)")
    #summary.seqs(fasta=robson.trim.contigs.good.unique.good.align, count=robson.trim.contigs.good.count_table)
    # remove overhangs on both ends
    m.filter.seqs(fasta=project_name +
                  ".trim.contigs.good.unique.good.align", vertical=True, trump=".")
    logger.info("Removed Overhangs: (15/38)")
    #summary.seqs(fasta=robson.trim.contigs.good.unique.good.filter.fasta, count=robson.trim.contigs.good.count_table)
    # check for duplicate sequences again
    m.unique.seqs(fasta=project_name+".trim.contigs.good.unique.good.filter.fasta",
                  count=project_name+".trim.contigs.good.count_table")
    logger.info("Removed Duplicate Seqs: (16/38)")
    # sets up the clustering
    m.pre.cluster(fasta=project_name+".trim.contigs.good.unique.good.filter.unique.fasta",
                  count=project_name+".trim.contigs.good.unique.good.filter.count_table", diffs=pre_clust_val)
    logger.info("Pre Clustered: (17/38)")
    # finds the chimeras
    m.chimera.vsearch(fasta=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.fasta",
                      count=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.count_table", dereplicate=True)
    logger.info("Detected Chimeras: (18/38)")
    # Gets rid of the chimeras found
    m.remove.seqs(fasta=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.fasta",
                  accnos=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos")
    logger.info("Removed Detected Chimeras: (19/38)")
    # checks out the file after seqs removed
    m.summary.seqs(fasta=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta",
                   count=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table")
    logger.info("Summarized Seqs: (20/38)")
    # assigns taxa to seqs
    m.classify.seqs(fasta=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta", count=project_name +
                    ".trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table", reference="trainset16_022016.pds.fasta", taxonomy="trainset16_022016.pds.tax", cutoff=80)
    logger.info("Classified Seqs: (21/38)")
    # removed taxa we don't care about
    m.remove.lineage(fasta=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta", count=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table",
                     taxonomy=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy", taxon="Chloroplast - Mitochondria - unknown - Archaea - Eukaryota")
    logger.info("Removed Taxa: (22/38)")
    m.count.groups(count=project_name +
                   ".trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table")
    logger.info("Counted Groups: (23/38)")

    #sub_samp_size = 15338
    if standard == True:
        logger.info("standard run")
        sub_sample_choices = []
        with open(project_name+".trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count.summary") as count_summary:

            contents = csv.reader(count_summary, delimiter='\t')

            for line in contents:
                try:
                    new_1 = int(float(line[1]))
                    sub_sample_choices.append(new_1)
                except ValueError:
                    pass
        sub_samp_size = sub_sample_choices[0]
        for val in sub_sample_choices:
            if val < sub_samp_size:
                sub_samp_size = val
    else:
        logger.info("custom run")
        #sub_samp_size = input("Please enter a value to Subsample at > ")
        sub_samp_size = sub_samp_size
    logger.info("Subsampling size: " + str(sub_samp_size))

    # Normalizes to a picked value
    m.sub.sample(fasta=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta", count=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table",
                 taxonomy=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy", size=sub_samp_size, persample=True)
    logger.info("Sub Sampled Based on Lowest Sample: (24/38)")
    # Check that it worked
    m.count.groups(count=project_name +
                   ".trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.subsample.count_table")
    logger.info("Counted Groups: (25/38)")
    # ake a look at the sequences
    m.summary.seqs(fasta=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.subsample.fasta",
                   count=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.subsample.count_table")
    logger.info("Summarized Seqs: (26/38)")
    # end of preprocessing
    logger.info("End of Preprocessing")

    # Creates a sequence Table and prints it to the terminal as well as to file with the filename seqTablePROJECTNAME.txt
    sequence_table = readTableNumpy.create_seq_table(project_name)
    readTableNumpy.print_seq_table(sequence_table)
    readTableNumpy.write_seq_table_to_file(sequence_table, project_name)

    # seqTable(project_name)

    # Rename files to more managable names
    m.rename.file(input=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.subsample.fasta",
                  new=project_name+".final.subs.fasta")
    logger.info("Renamed Subsampled Seqs: (27/38)")
    m.rename.file(input=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.subsample.count_table",
                  new=project_name+".final.subs.count_table")
    logger.info("Renamed Count Table: (28/38)")
    m.rename.file(input=project_name+".trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.subsample.taxonomy",
                  new=project_name+".final.subs.taxonomy")
    logger.info("Renamed Taxonomy File: (29/38)")
    # create a taxonomy summary
    m.summary.tax(taxonomy=project_name+".final.subs.taxonomy",
                  count=project_name+".final.subs.count_table")
    logger.info("Summarized Taxonomy: (30/38)")
    # bin sequences to otu
    m.cluster.split(fasta=project_name+".final.subs.fasta", count=project_name+".final.subs.count_table",
                    taxonomy=project_name+".final.subs.taxonomy", splitmethod="classify", taxlevel=4, cutoff=0.03)
    logger.info("Binned Sequences to OTU: (31/38)")
    # determine taxonomy for each OTU
    m.classify.otu(list=project_name+".final.subs.opti_mcc.list", count=project_name +
                   ".final.subs.count_table", taxonomy=project_name+".final.subs.taxonomy", label=0.03)
    logger.info("Determied Tax for OTU: (32/38)")
    # Make a shred file
    m.make.shared(list=project_name+".final.subs.opti_mcc.list",
                  count=project_name+".final.subs.count_table", label=0.03)
    logger.info("Made Shared File: (33/38)")
    # Asses Alpha diversity
    m.summary.single(shared=project_name+".final.subs.opti_mcc.shared",
                     calc="nseqs - coverage - sobs - shannon")
    logger.info("Assesed Alpha Diversity: (34/38)")
    # calculate distance between samples
    m.dist.shared(shared=project_name+".final.subs.opti_mcc.shared",
                  calc="braycurtis - thetayc")
    logger.info("Calculated Distance Between Samples: (35/38)")
    # create ordination plot
    m.pcoa(phylip=project_name+".final.subs.opti_mcc.braycurtis.0.03.lt.dist")
    logger.info("Made Ordination Plot: (36/38)")
    # Check if groups are significantly different
    m.amova(phylip=project_name +
            ".final.subs.opti_mcc.braycurtis.0.03.lt.dist", design=design)
    logger.info("Determined if Groups are Significantly Different: (37/38)")
    # OTUs sig different between groups
    m.metastats(shared=project_name +
                ".final.subs.opti_mcc.shared", design=design)
    logger.info("Determined if OTU are Different Between Groups: (38/38)")
