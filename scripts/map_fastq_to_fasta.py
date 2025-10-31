import os
import subprocess
import csv
import argparse
from multiprocessing import cpu_count
import re
from collections import defaultdict


def extract_species(header):
    pattern = r">\S+\s+([A-Z][a-z]+\s[a-z]+)"
    match = re.search(pattern, header)
    if match:
        return match.group(1)
    return None


def parse_cigar(cigar):
    if cigar == "*" or not cigar:
        return {op: 0 for op in "MIDNSHP=X"}
    counts = {op: 0 for op in "MIDNSHP=X"}
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        counts[op] += int(length)
    return counts


def get_optional_tag(tags, key):
    prefix = key + ":"
    for t in tags:
        if t.startswith(prefix):
            parts = t.split(":", 2)
            if len(parts) == 3:
                _k, _type, val = parts
                if _type == "i":
                    try:
                        return int(val)
                    except ValueError:
                        return None
                return val
    return None


def is_unmapped(flag):
    return (flag & 0x4) != 0


def is_secondary_or_supp(flag):
    return (flag & 0x100) != 0 or (flag & 0x800) != 0


def map_reads(
    target_folder,
    query_folder,
    threads,
    software,
    output_dir,
    min_identity,
    max_size_diff,
    min_mapq=20
):
    os.makedirs(output_dir, exist_ok=True)

    if software not in ["minimap2", "pblat"]:
        raise ValueError("Unsupported software. Choose 'minimap2' or 'pblat'.")

    summary_file = os.path.join(output_dir, "mapping_summary.csv")

    with open(summary_file, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Query_File",
                             "Target_File",
                             "Unique_Mappers",
                             "Total_Mappers"])

        for target_file in os.listdir(target_folder):
            if not target_file.endswith(".fasta"):
                continue
            target_path = os.path.join(target_folder, target_file)

            for query_file in os.listdir(query_folder):
                if not query_file.endswith(".fastq"):
                    continue
                query_path = os.path.join(query_folder, query_file)

                output_file = os.path.join(
                    output_dir,
                    f"{os.path.splitext(query_file)[0]}_to_{os.path.splitext(target_file)[0]}.txt"
                )

                if software == "minimap2":
                    cmd = [
                        "minimap2",
                        "-ax", "map-ont",
                        "-t", str(threads),
                        target_path,
                        query_path
                    ]
                else:
                    cmd = [
                        "pblat",
                        target_path,
                        query_path,
                        output_file
                    ]

                try:
                    if software == "minimap2":
                        with open(output_file, "w") as output:
                            subprocess.run(
                                cmd,
                                stdout=output,
                                stderr=subprocess.PIPE,
                                check=True
                            )
                    else:
                        subprocess.run(
                            cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            check=True
                        )
                except subprocess.CalledProcessError as e:
                    print(f"Error running {software} on {query_file} vs {target_file}: {e}")
                    csv_writer.writerow([query_file, target_file, 0, 0])
                    continue

                query_to_targets = defaultdict(set)

                if software == "minimap2":
                    with open(output_file) as f:
                        for line in f:
                            if not line or line[0] == "@":
                                continue
                            cols = line.rstrip("\n").split("\t")
                            if len(cols) < 11:
                                continue

                            qname = cols[0]
                            try:
                                flag = int(cols[1])
                            except ValueError:
                                continue
                            rname = cols[2]
                            try:
                                mapq = int(cols[4])
                            except ValueError:
                                mapq = 0
                            cigar = cols[5]
                            seq = cols[9]
                            tags = cols[11:] if len(cols) > 11 else []

                            if is_unmapped(flag):
                                continue
                            if is_secondary_or_supp(flag):
                                continue
                            if rname == "*":
                                continue
                            if mapq < min_mapq:
                                continue

                            c = parse_cigar(cigar)
                            read_aln_len = c["M"] + c["="] + c["X"] + c["I"]
                            target_aln_len = c["M"] + c["="] + c["X"] + c["D"]

                            if read_aln_len == 0 and seq and seq != "*":
                                read_aln_len = len(seq)

                            nm = get_optional_tag(tags, "NM")
                            identity_ok = True
                            if nm is not None and read_aln_len > 0:
                                identity = max(0.0, (1.0 - (nm / read_aln_len)) * 100.0)
                                identity_ok = identity >= min_identity

                            size_ok = True
                            if read_aln_len > 0 and target_aln_len > 0:
                                size_diff = abs(read_aln_len - target_aln_len) / max(read_aln_len, target_aln_len) * 100.0
                                size_ok = size_diff <= max_size_diff

                            if identity_ok and size_ok:
                                species = extract_species(rname)
                                query_to_targets[qname].add(sometimes_alias_species(species) if False else species)

                else:
                    with open(output_file) as f:
                        lines = f.readlines()[5:]

                    for line in lines:
                        columns = line.split()
                        if len(columns) < 17:
                            continue

                        matches = int(columns[0])
                        mismatches = int(columns[1])
                        qname = columns[9]
                        qsize = int(columns[10])
                        tname_raw = columns[13]
                        tsize = int(columns[14])

                        identity = matches / (matches + mismatches) * 100
                        size_diff = abs(qsize - tsize) / max(qsize, tsize) * 100

                        if identity >= min_identity and size_diff <= max_size_diff:
                            species = extract_species(tname_raw)
                            query_to_targets[qname].add(species)

                unique_mappers = sum(1 for targets in query_to_targets.values() if len(targets) == 1)
                total_mappers = len(query_to_targets)

                csv_writer.writerow([
                    query_file,
                    target_file,
                    unique_mappers,
                    total_mappers
                ])

    return summary_file


def main():
    parser = argparse.ArgumentParser(description="Map query reads to target sequences and evaluate unique mappers.")
    parser.add_argument("--target_folder", required=True, help="Folder containing target fasta sequences.")
    parser.add_argument("--query_folder", required=True, help="Folder containing fastq query files.")
    parser.add_argument("--threads", type=int, default=cpu_count(), help="Number of threads to use (default: all available).")
    parser.add_argument("--software", choices=["minimap2", "pblat"], required=True, help="Software to use for mapping.")
    parser.add_argument("--output_dir", required=True, help="Directory to store output files.")
    parser.add_argument("--min_identity", type=float, default=95.0, help="Minimum % identity to keep an alignment (pblat only).")
    parser.add_argument("--max_size_diff", type=float, default=5.0, help="Maximum % length difference between query and target (pblat only).")

    args = parser.parse_args()

    summary_file = map_reads(
        target_folder=args.target_folder,
        query_folder=args.query_folder,
        threads=args.threads,
        software=args.software,
        output_dir=args.output_dir,
        min_identity=args.min_identity,
        max_size_diff=args.max_size_diff
    )

    print(f"Mapping summary written to: {summary_file}")

if __name__ == "__main__":
    main()

