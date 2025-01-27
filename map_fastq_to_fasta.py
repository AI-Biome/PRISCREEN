import os
import subprocess
import csv
import argparse
from multiprocessing import cpu_count

def map_reads(target_folder, query_folder, threads, software, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    if software not in ["minimap2", "pblat"]:
        raise ValueError("Unsupported software. Choose 'minimap2' or 'pblat'.")

    summary_file = os.path.join(output_dir, "mapping_summary.csv")

    with open(summary_file, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Query_File", "Target_File", "Unique_Mappers", "Total_Mappers"])

        for target_file in os.listdir(target_folder):
            if not target_file.endswith(".fasta"):
                continue
            target_path = os.path.join(target_folder, target_file)

            for query_file in os.listdir(query_folder):
                if not query_file.endswith(".fastq"):
                    continue
                query_path = os.path.join(query_folder, query_file)
                output_file = os.path.join(output_dir, f"{os.path.splitext(query_file)[0]}_to_{os.path.splitext(target_file)[0]}.txt")

                if software == "minimap2":
                    cmd = [
                        "minimap2",
                        "-ax", "map-ont",
                        f"--threads={threads}",
                        target_path,
                        query_path
                    ]
                elif software == "pblat":
                    cmd = [
                        "pblat",
                        target_path,
                        query_path,
                        output_file
                    ]

                try:
                    with open(output_file, "w") as output:
                        subprocess.run(cmd, stdout=output, stderr=subprocess.PIPE, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error running {software} on {query_file} vs {target_file}: {e}")
                    continue

                # todo: evaluate output

    return summary_file

def main():
    parser = argparse.ArgumentParser(description="Map query reads to target sequences and evaluate unique mappers.")
    parser.add_argument("--target_folder", required=True, help="Folder containing target fasta sequences.")
    parser.add_argument("--query_folder", required=True, help="Folder containing fastq query files.")
    parser.add_argument("--threads", type=int, default=cpu_count(), help="Number of threads to use (default: all available).")
    parser.add_argument("--software", choices=["minimap2", "pblat"], required=True, help="Software to use for mapping.")
    parser.add_argument("--output_dir", required=True, help="Directory to store output files.")

    args = parser.parse_args()

    summary_file = map_reads(
        target_folder=args.target_folder,
        query_folder=args.query_folder,
        threads=args.threads,
        software=args.software,
        output_dir=args.output_dir
    )

    print(f"Mapping summary written to: {summary_file}")

if __name__ == "__main__":
    main()

