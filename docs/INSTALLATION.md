# GermCall: installation guide

## 1. conda (recommended)
```bash
# install conda once (skip if you have conda/mamba)
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh          # restart the terminal afterwards

cd GermCall
conda config --set channel_priority strict        # avoids mixed-channel library errors
conda env create -f environment.yml               # 10-30 minutes
conda activate germcall
```
Check:
```bash
fastqc --version; multiqc --version; trim_galore --version | head -3; cutadapt --version
bwa 2>&1 | head -3; samtools --version | head -1; bcftools --version | head -1; gatk --version; java -version
```
GATK4 needs Java 17, which the environment provides.

## 2. Without conda (Ubuntu / Debian)
```bash
sudo apt update
sudo apt install fastqc trim-galore cutadapt bwa samtools bcftools tabix openjdk-17-jre-headless python3-yaml
pip install multiqc
```
GATK is not in the Ubuntu repositories: download the latest GATK4 release (`gatk-4.x.zip`) from https://github.com/broadinstitute/gatk/releases, unzip it, and either add its folder to your `PATH` or set `tools.gatk` to the full path of the `gatk` launcher in the config. (GATK's launcher also needs Python 3.)

## 3. What each stage needs
| Stage | Tools |
|---|---|
| prepare | bwa, samtools, GATK, bcftools, htslib (bgzip, tabix), Trim Galore (+ cutadapt) |
| qc | FastQC |
| trim_align | Trim Galore, (FastQC), BWA, samtools, GATK (MarkDuplicates) |
| gvcf / joint / filter | GATK, bcftools, htslib |
| report | MultiQC |

Python 3 with PyYAML is always needed (the scripts read the YAML config).

## 4. Test the installation
```bash
conda activate germcall
bash examples/run_test.sh        # simulated data, a few minutes; creates ./germcall_test/ (safe to delete)
```
It runs every stage on four simulated samples and ends with `All requested stages finished.` If a tool is missing, the error names it. The scripts were developed and used on Linux (Ubuntu); macOS and Windows (WSL) are not tested.

## 5. Common problems
| Problem | Fix |
|---|---|
| `ModuleNotFoundError: dnaio._core` when Trim Galore runs | a broken system `cutadapt`: use the conda environment (or `pip install cutadapt` in a virtual environment) |
| `PyYAML is required` | `pip install pyyaml` or `conda install pyyaml` |
| Library errors such as `libicui18n.so.58` or `GLIBCXX` | packages from different conda channels got mixed: create the environment with `conda config --set channel_priority strict` (this is a global conda setting) and without the `defaults` channel, as `environment.yml` does |
| conda solving takes very long or fails | run `conda config --set channel_priority strict`, install `mamba`, and use `mamba env create -f environment.yml` |
| GATK: `A USER ERROR has occurred: ... reference ... dictionary` | the `.dict` file is missing: run `--steps prepare` (the reference folder must be writable) |
