with open("AF.imputed.thin.bim", "r") as file:
    lines = file.readlines()

newlines = []
for line in lines:
    columns = line.split()
    chromosome_id = columns[0]
    chromosome = int(chromosome_id[7:9]) - 23
    newlines.append(line.replace(chromosome_id, str(chromosome)))

with open("AF.imputed.thin.bim", "w") as file:
    file.writelines(newlines)