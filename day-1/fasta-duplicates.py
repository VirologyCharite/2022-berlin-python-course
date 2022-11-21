import sys

ids = set()

for line in sys.stdin:
    if line.startswith(">"):
        if line in ids:
            print("Found duplicate:", line, end="")
        else:
            ids.add(line)
    else:
        print("SEQ", line)
