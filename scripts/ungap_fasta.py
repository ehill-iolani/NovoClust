from pathlib import Path
import sys

inp = Path(sys.argv[1])
out = Path(sys.argv[2])

name = None
seq = []

def flush(nm, sq, handle):
    if nm is None:
        return
    s = "".join(sq).replace("-", "").replace(".", "")
    handle.write(f">{nm}\n{s}\n")

with inp.open() as f, out.open("w") as out_handle:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            flush(name, seq, out_handle)
            name = line[1:].strip()
            seq = []
        else:
            seq.append(line)
    flush(name, seq, out_handle)