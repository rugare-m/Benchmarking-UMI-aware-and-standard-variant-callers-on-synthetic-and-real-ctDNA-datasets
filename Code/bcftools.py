#works
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

def process_file(alignments_file):
    sample_name = os.path.basename(alignments_file).split(".")[0]
    reference_file = "/users/rugarem/volatile/chapter3/data/reference/hg38.fa.gz"
    output_file = f"bcftools.{sample_name}.vcf.gz"
    log_file = f"bcftools.{sample_name}.log"  # Log file name for caller output

    command1 = ["bcftools", "mpileup", "-Ob", "-f", reference_file, alignments_file]
    command2 = ["bcftools", "call", "-mv", "-Oz", "-o", output_file]

    with open(log_file, "w") as log:
        process1 = subprocess.Popen(command1, stdout=subprocess.PIPE, stderr=log)
        process2 = subprocess.Popen(command2, stdin=process1.stdout, stderr=log)
        process1.stdout.close()
        process2.communicate()

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
