from Bio import SeqIO

primers = list(SeqIO.parse("data/SARS-COV-2_spikePrimers.fasta", "fasta"))

primers_rc = []
for s in primers:
    primers_rc.append("".join(s.reverse_complement()))

ref, = SeqIO.parse("data/MN908947.3.fa", "fasta")
ref_str = "".join(ref)

assert all([s in ref_str for s in primers_rc])

primer_pos = [ref_str.find(s) for s in primers_rc]

chrom = "MN908947.3"
with open("data/SARS-COV-2_spikePrimers.bed", "w") as f:
    for i in range(len(primers)):
        name = primers
        rc = primers_rc[i]
        pos = primer_pos[i]
        print(chrom, pos, pos + len(rc),
              primers[i].name, 60, "-", file=f, sep="\t")
