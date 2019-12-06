import csv
project_name = "test"

start = 0
end = 0
count = 0
with open(project_name+".trim.contigs.good.unique.summary") as summary:
    start = 0
    end = 0
    count = 0
    contents = csv.reader(summary, delimiter='\t')

    for line in contents:
        try:
            new1 = int(float(line[1]))
            new2 = int(float(line[2]))
            start = start+new1
            end = end+new2
            count = count + 1
        except ValueError:
            pass
        # print(line[1])
        # print(line[2])

    #print(str(int(start/count)) + " " + str(int(end/count)))
subSampleChoices = []
with open(project_name+".trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count.summary") as countsummary:

    contents = csv.reader(countsummary, delimiter='\t')

    for line in contents:
        try:
            new1 = int(float(line[1]))
            subSampleChoices.append(new1)
        except ValueError:
            pass
    for val in subSampleChoices:
        pass
        # print(val)

# project_name+".trim.contigs.summary" (assembly 3), project_name+".trim.contigs.good.summary" (Ambig 5), project_name+".trim.contigs.good.unique.summary" (Unique 8), project_name+".trim.contigs.good.unique.summary" (Aligned 13),
# project_name+".trim.contigs.good.unique.good.filter.unique.precluster.pick.summary" (Chimera 20)
readsTable = []
readsTable.append(["Step", "No. Seq", "% Total Seq"])
assmb3 = 0
with open(project_name+".trim.contigs.summary") as assembly3:
    contents = csv.reader(assembly3, delimiter='\t')

    for line in contents:
        try:
            assmb3 = assmb3 + int(float(line[6]*100))

        except ValueError:
            pass
readsTable.append(["Assembly 3", assmbl3, float(assmb3/assmb3)*100])
am5 = 0
with open(project_name+".trim.contigs.good.summary") as ambig5:
    contents = csv.reader(ambig5, delimiter='\t')

    for line in contents:
        try:
            am5 = am5 + int(float(line[6]*100))

        except ValueError:
            pass
readsTable.append(["Ambigous 5", am5, float(am5/assmb3)*100])

unq8 = 0
with open(project_name+".trim.contigs.good.unique.summary") as unique8:
    contents = csv.reader(unique8, delimiter='\t')

    for line in contents:
        try:
            unq8 = unq8 + int(float(line[6]))

        except ValueError:
            pass
readsTable.append(["Unique 8", unq8, float(unq8/assmb3)*100])

alig13 = 0
with open(project_name+".trim.contigs.good.unique.summary") as aligned13:
    contents = csv.reader(aligned13, delimiter='\t')

    for line in contents:
        try:
            alig13 = alig13 + int(float(line[6]))

        except ValueError:
            pass
readsTable.append(["Aligned 13", alig13, float(alig13/assmb3)*100])


with open(project_name+".trim.contigs.good.unique.good.filter.unique.precluster.pick.summary") as chimera20:
    contents = csv.reader(chimera20, delimiter='\t')

    for line in contents:
        try:
            chim20 = chim20 + int(float(line[6]))

        except ValueError:
            pass
readsTable.append(["Chimera 20", chim20, float(chim20/assmb3)*100])
