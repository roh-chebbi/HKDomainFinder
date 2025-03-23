import sys

file_name = sys.argv[1]

sequences=[]

f = open(file_name, "r")
for line in f:
    if not line.startswith(">"):
        sequences.append(line.rstrip())

matches=[]
for i in range(len(sequences[1])):
    if sequences[0][i] == sequences[1][i]:
        matches.append("|")
    else:
        matches.append(" ")

print(sequences[0])
print("".join(matches))
print(sequences[1])