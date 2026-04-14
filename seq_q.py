from Bio import SeqIO
from statistics import mean
import sys

def sliding_window_lowq(quals, window=int(sys.argv[3]), step=int(sys.argv[4]), cutoff=int(sys.argv[2])):
    results = []
    for i in range(0, len(quals) - window + 1, step):
        chunk = quals[i:i+window]
        mean_q = sum(chunk) / window
        if mean_q < cutoff:
            results.append((i+1, i+window, mean_q))
    return results

def main(fastq_file):
    for record in SeqIO.parse(fastq_file, "fastq"):
        quals = record.letter_annotations["phred_quality"]
        read_len = len(quals)
        midpoint = read_len / 2
        mean_q= mean(quals)
        windows = sliding_window_lowq(quals)

        last_before = None
        first_after = None

        for start, end, m in windows:
            center = (start + end) / 2

            if center <= midpoint:
                last_before = end  # pega a última base antes da metade
            elif center > midpoint and first_after is None:
                first_after = start  # primeira base depois da metade
        print(f"id\tread_len\tmean_q\tlast_before\tfirst_after")
        print(f"{record.id}\t{read_len}\t{mean_q}\t{last_before}\t{first_after}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Uso: python script.py arquivo.fastq Q_cutoff window step")
        sys.exit(1)

    main(sys.argv[1])
