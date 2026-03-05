from pathlib import Path
import sys

aln_in = Path(sys.argv[1])
aln_out = Path(sys.argv[2])

seqs = []
names = []
cur_name = None
cur_seq = []

for line in aln_in.read_text().splitlines():
    if line.startswith(">"):
        if cur_name is not None:
            names.append(cur_name)
            seqs.append("".join(cur_seq))
        cur_name = line[1:].strip()
        cur_seq = []
    else:
        cur_seq.append(line.strip())

if cur_name is not None:
    names.append(cur_name)
    seqs.append("".join(cur_seq))

if not seqs:
    raise SystemExit("No sequences found in alignment.")

L = len(seqs[0])
if any(len(s) != L for s in seqs):
    raise SystemExit("Alignment sequences are not the same length. MSA failed?")

def is_gap(c): return c in "-."

# trim only ends: first/last column where ALL sequences are non-gap
left = next((i for i in range(L) if all(not is_gap(s[i]) for s in seqs)), None)
right = next((i for i in range(L-1, -1, -1) if all(not is_gap(s[i]) for s in seqs)), None)

if left is None or right is None or right <= left:
    raise SystemExit("Could not find a shared non-gapped core window.")

with aln_out.open("w") as f:
    for n, s in zip(names, seqs):
        f.write(f">{n}\n{s[left:right+1]}\n")