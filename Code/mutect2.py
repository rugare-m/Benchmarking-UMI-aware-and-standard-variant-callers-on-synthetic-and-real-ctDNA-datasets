#works
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

def process_file(alignments_file):
    sample_name = os.path.basename(alignments_file).split(".")[0]
    reference_file = "/users/rugarem/volatile/chapter3/data/reference/hg38.fa.gz"
    output_file = f"mutect2.{sample_name}.vcf"
    filtered = f"filtered_mutect2.{sample_name}.vcf"
    log_file = f"mutect2.{sample_name}.log"  # Log file name for caller output

    mutect2_command = ["gatk", "Mutect2", "-R", reference_file, "-I", alignments_file, "-O", output_file]
    compress = ["bgzip", output_file]
    mutect2_filter = ["gatk", "FilterMutectCalls", "-R", reference_file, "-V", output_file, "-O", filtered]
    
    with open(log_file, "w") as log:
        process = subprocess.run(mutect2_command, stdout=subprocess.PIPE, stderr=log)
        process2 = subprocess.run(mutect2_filter, stdout=subprocess.PIPE, stderr=log)

    subprocess.run(compress, stdout=subprocess.PIPE)
bam_files = [f for f in os.listdir() if f.endswith(".bqsr.bam")]

with ProcessPoolExecutor(max_workers=3) as executor:
    futures = [executor.submit(process_file, bam_file) for bam_file in bam_files]
    for future in as_completed(futures):
        pass


if not os.path.exists("logs"):
    os.makedirs("logs")

for log_file in os.listdir():
    if log_file.endswith(".log"):
        os.rename(log_file, os.path.join("logs", log_file))
