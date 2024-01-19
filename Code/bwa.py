#works, parallel processing, 3 subprocesses at a time
import concurrent.futures
import subprocess
import glob
import os

# Map reads to reference SRR3401416_1.fastq
input_dir = "./"
input_pattern = "SRR*_1.fastq.gz"
reference_fa = "/users/rugarem/volatile/chapter3/data/reference/hg38"
threads = "64"

fastq_files = glob.glob(f"{input_dir}/{input_pattern}")

def map_reads(fastq):
    fastq_prefix = fastq[:-11]  # Remove "_1.fastq.gz" from the file name
    output_sam = f"{fastq_prefix}.sam"

    bwamem = ["bwa-mem2", "mem", "-t", threads, reference_fa, fastq, f"{fastq_prefix}_2.fastq.gz"]

    print(f"Mapping reads from {fastq} to reference...", flush=True)
    with open(output_sam, "w") as outfile:
        subprocess.run(bwamem, check=True, stdout=outfile)
    print(f"Mapping of reads from {fastq} completed.", flush=True)

with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:
    executor.map(map_reads, fastq_files)

print("All reads mapped to the reference.", flush=True)


#convert sam to bam
sam_files = glob.glob("*.sam")

def process_sam(sam_file):
    bam_file = os.path.splitext(sam_file)[0] + ".bam"
    sorted_bam_file = os.path.splitext(sam_file)[0] + ".sort.bam"

    # Convert SAM to BAM
    samtools_convert_cmd = ["samtools", "view", "-@", threads, "-b", "-o", bam_file, sam_file]
    print(f"Converting {sam_file} to BAM...", flush=True)
    subprocess.run(samtools_convert_cmd, check=True)
    print(f"Conversion of {sam_file} to BAM completed.", flush=True)

    # Sort the BAM file
    samtools_sort_cmd = ["samtools", "sort", "-o", sorted_bam_file, bam_file]
    print(f"Sorting BAM file {bam_file}...", flush=True)
    subprocess.run(samtools_sort_cmd, check=True)
    print(f"Sorting of {bam_file} completed.", flush=True)

    # Delete the unsorted BAM file
    os.remove(bam_file)
    print(f"Deleted unsorted BAM file: {bam_file}", flush=True)

    # Delete the SAM file
    os.remove(sam_file)
    print(f"Deleted SAM file: {sam_file}", flush=True)  

with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:
    executor.map(process_sam, sam_files)

print("All SAM files processed.", flush=True)
