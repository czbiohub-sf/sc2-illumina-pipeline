from Bio import SeqIO

primers = list(SeqIO.parse("data/SARS-COV-2_spikePrimers.fasta", "fasta"))

ref, = SeqIO.parse("data/MN908947.3.fa", "fasta")
ref_str = "".join(ref)

primer_pos = []
for s in primers:
    rc = "".join(s.reverse_complement())
    assert rc in ref_str
    primer_pos.append(ref_str.find(rc))

sorted_pos_seq = sorted(zip(primer_pos, primers), key=lambda x: x[0])

chrom = "MN908947.3"
with open("data/SARS-COV-2_spikePrimers.bed", "w") as f:
    for pos, seq in sorted_pos_seq:
        print(chrom, pos, pos + len(seq),
              seq.name, 60, "+", file=f, sep="\t")
