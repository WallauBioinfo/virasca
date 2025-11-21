import csv

a = input()
b = input()

if a == "":
    a = "./out.emapper.annotations"
if b == "":
    b = "./out.emapper.decorated.gff"

ann = []
f = open(a)
reader = f.readlines()
for row in reader:
    if row[0] == "#":
        continue
    row = row.split('\t')
    if row[7] != '-':
        temp = row[0].split('_')
        ann.append([temp[len(temp) - 1], row[7]])

vis = []
for i in range(len(ann)):
    vis.append(ann[i][0])

f = open(b, encoding = "utf-8")
s = f.readlines()
s = s[4:]
x = []
final = []
for i in range(len(s)):
    s[i] = s[i].strip()
    temp = s[i].split('\t')
    x.append([temp[3], temp[4], temp[6], temp[8]])
for i in range(len(x)):
    description = ""
    if "em_ID" in x[i][3]:
        x[i][3] = x[i][3].split(';')
        temp = x[i][3][13]
        temp = temp.split('_')
        val = temp[len(temp) - 1]
        if val in vis:
            description = ann[vis.index(val)][1]
    final.append([description, "CDS", x[i][0], x[i][1], x[i][2]])
final_name = []
final_no_name = []
for i in range(len(final)):
    if final[i][0] == "":
        final_no_name.append(final[i])
    else:
        final_name.append(final[i])

# Generate the annotation file with annotated genes.
with open('final.csv', 'w', newline='') as csvfile:
    fieldnames = ["name", "type", "start", "stop", "strand"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for i in range(len(final_name)):
        writer.writerow({"name": final_name[i][0], "type": final_name[i][1], "start": final_name[i][2], "stop": final_name[i][3], "strand": final_name[i][4]})

# Generate the annotation file with unannotated genes.
with open('final2.csv', 'w', newline='') as csvfile:
    fieldnames = ["name", "type", "start", "stop", "strand"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for i in range(len(final_no_name)):
        writer.writerow({"name": final_no_name[i][0], "type": final_no_name[i][1], "start": final_no_name[i][2], "stop": final_no_name[i][3], "strand": final_no_name[i][4]})