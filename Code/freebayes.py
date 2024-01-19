#works
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

def process_file(alignments_file):
    sample_name = os.path.basename(alignments_file).split(".")[0]
    reference_file = "/users/rugarem/volatile/chapter3/data/reference/uhg38.fa"
    #reference_file = "/users/rugarem/volatile/chapter3/data/reference/hg38.fa.gz"
    output_file = f"fbayes.{sample_name}.vcf"
    log_file = f"fbayes.{sample_name}.log"  # Log file name for caller output

    fbayes_command = ["freebayes", "-f", reference_file, alignments_file]

    with open(log_file, "w") as log:
        process = subprocess.Popen(fbayes_command, stdout=open(output_file, "w"), stderr=log)
        process.communicate()
    
    process2 = subprocess.Popen(["bgzip", output_file], stdout=subprocess.PIPE)

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
