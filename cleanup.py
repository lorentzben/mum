# Cleanup.py A small script to remove files used for testing but not dev
# Ben Lorentz 11.20.2019
import os
files = ["oligos.txt", "Silva.nr_v132.tgz", "trainset16_022016.rdp", "*.fastq",
         "design.txt", "Kelly*", "Test*", "current_files.summary", "silva*", "trainset*"]

for item in files:
    # print(item)
    command = "rm -rf "+item
    os.system(command)

print("All Clean")
