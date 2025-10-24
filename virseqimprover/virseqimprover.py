#!/usr/bin/env python3
"""
Refactored virseqimprover CLI.
- Uses pathlib for safe path handling.
- Uses subprocess.run([...]) with capture for safer execution.
- Logs commands and outputs to run.log
- Stops immediately on any subprocess error (check=True -> exception)

Note: This is a functional refactor preserving the original pipeline steps and
naming. Some tiny behavioral differences may exist where the original relied on
shell side-effects; those parts use "bash -c" with cwd set to ensure parity.
"""

from __future__ import annotations

import math
import subprocess
import sys
import shlex
from pathlib import Path
from typing import List, Optional, Tuple

import click
from virseqimprover import __version__


# --- Configuration / Globals -------------------------------------------------

DEFAULT_SPADES_KMER = "default"
MIN_IDENTITY_CIRCULAR = 95
SALMON_READ_FRACTION = 0
MIN_SUSPICIOUS_LEN = 1000


# --- Utilities ---------------------------------------------------------------

class RunError(RuntimeError):
    pass


def run_cmd(cmd: List[str] | str, cwd: Optional[Path] = None, log_file: Optional[Path] = None) -> str:
    """Run a command safely.  If `cmd` is a list it's passed directly to subprocess.
    If `cmd` is a str it is executed with `bash -c` to preserve complex shell pipelines.

    On failure this function re-raises subprocess.CalledProcessError as RunError
    so the pipeline stops immediately (as requested).
    """
    if isinstance(cmd, list):
        proc_args = cmd
    else:
        # run complex shell command via bash -c
        proc_args = ["bash", "-lc", cmd]

    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)

    try:
        result = subprocess.run(
            proc_args,
            cwd=str(cwd) if cwd else None,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        if log_file:
            with open(log_file, "a") as fh:
                fh.write(f"$ {' '.join(proc_args)}\n")
                if result.stdout:
                    fh.write(result.stdout + "\n")
                if result.stderr:
                    fh.write("STDERR:\n" + result.stderr + "\n")

        return result.stdout.strip()

    except subprocess.CalledProcessError as e:
        # write failing command and stderr to log if provided
        if log_file:
            with open(log_file, "a") as fh:
                fh.write(f"FAILED: {' '.join(proc_args)}\n")
                fh.write(e.stderr or "(no stderr captured)")
        raise RunError(f"Command failed: {' '.join(proc_args)}\n{e.stderr}") from e


def safe_read_text(path: Path) -> str:
    with path.open("r") as fh:
        return fh.read()


def safe_read_first_line(path: Path) -> str:
    with path.open("r") as fh:
        return fh.readline()


# --- Pipeline functions ------------------------------------------------------


def get_read_length(read1: Path, output_dir: Path, log_file: Path):
    """
    Calculate average read length safely for gzipped or plain FASTQ.
    Stops immediately if computation fails.
    """
    quoted_path = shlex.quote(str(read1.resolve()))
    if read1.suffix == ".gz" or read1.suffixes[-2:] == [".fastq", ".gz"]:
        cmd = [
            "bash",
            "-c",
            f"gzip -dc -- {quoted_path} | awk 'NR%4==2{{len+=length($0); n++}} END{{if (n>0) print len/n; else print 0}}'"
        ]
    else:
        cmd = [
            "bash",
            "-c",
            f"awk 'NR%4==2{{len+=length($0); n++}} END{{if (n>0) print len/n; else print 0}}' {quoted_path}"
        ]

    out = run_cmd(cmd, cwd=output_dir, log_file=log_file)
    try:
        avg_len = float(out)
    except ValueError:
        print(f"[ERROR] Could not extract average read length from {read1}")
        avg_len = 0.0

    if avg_len <= 0:
        raise RuntimeError(
            f"Average read length computed as 0. Check input FASTQ file: {read1}"
        )

    return avg_len


def create_bed(output_dir: Path, avg_read_len: float, log: Path) -> None:
    """Create scaffold-start-end.bed from scaffold.fasta.fai located in scaffold-truncated or root.
    Also moves spades results if present (preserves original behavior).
    """
    scaffold_tmp_spades = output_dir / "scaffold-truncated" / "tmp" / "spades-res" / "scaffold.fasta"
    dest_scaffold_truncated = output_dir / "scaffold-truncated"

    # If spades produced scaffolds inside tmp/spades-res and present, move them up.
    if scaffold_tmp_spades.exists():
        # use bash to preserve several mv/rm operations in one atomic command, similar to original
        cmd = (
            f"cd {str(output_dir/'scaffold-truncated')} && "
            "rm -f scaffold.fasta* && "
            "cd tmp/spades-res && "
            f"mv scaffold.fasta {str(dest_scaffold_truncated)} && "
            f"mv scaffold.fasta.fai {str(dest_scaffold_truncated)}"
        )
        run_cmd(cmd, cwd=output_dir, log_file=log)

    fai_path = dest_scaffold_truncated / "scaffold.fasta.fai"
    if not fai_path.exists():
        raise RunError(f"Missing FAI index: {fai_path}")

    line = safe_read_first_line(fai_path)
    cols = line.split("\t")
    scaffold_id = cols[0].strip()
    scaffold_length = int(cols[1].strip())

    bed_path = dest_scaffold_truncated / "scaffold-start-end.bed"
    with bed_path.open("w") as bw, (output_dir / "output-log.txt").open("a") as logfh:
        if scaffold_length >= avg_read_len * 2:
            bw.write(f"{scaffold_id}\t0\t{math.ceil(avg_read_len*1.5)}\n")
            bw.write(f"{scaffold_id}\t{math.ceil(scaffold_length - avg_read_len*1.5)}\t{scaffold_length}\n")
            logfh.write(f"Trying to grow scaffold {scaffold_id} with length {scaffold_length}\n")
            if scaffold_length > 300000:
                logfh.write(f"Length of {scaffold_id} is already greater than 300kbp, so stop extending this one.\n")


def run_alignment(output_dir: Path, read1: Path, read2: Optional[Path], salmon_fraction: int, log: Path) -> None:
    """Run salmon quant + mapping steps used in original script. Results in salmon-mapped.sam in scaffold-truncated.
    We use bash -c to preserve LD_PRELOAD / pipes present in original script.
    """
    cwd = output_dir / "scaffold-truncated"
    # ensure target dir exists
    cwd.mkdir(parents=True, exist_ok=True)

    if read2 is None:
        reads_arg = f"-r {str(read1)}"
    else:
        reads_arg = f"-1 {str(read1)} -2 {str(read2)}"

    cmd = (
        "bedtools getfasta -fi scaffold.fasta -bed scaffold-start-end.bed -fo scaffold-start-end.fasta && "
        "rm -rf salmon-index salmon-res salmon-mapped.sam && "
        "LD_PRELOAD=/lib/x86_64-linux-gnu/libpthread.so.0:/lib/x86_64-linux-gnu/librt.so.1 "
        f"salmon index -t scaffold-start-end.fasta -i salmon-index && "
        "LD_PRELOAD=/lib/x86_64-linux-gnu/libpthread.so.0:/lib/x86_64-linux-gnu/librt.so.1 "
        f"salmon quant -i salmon-index -l A {reads_arg} -o salmon-res --writeMappings -p 16 --quasiCoverage {salmon_fraction} | "
        "samtools view -bS - | samtools view -h -F 0x04 - > salmon-mapped.sam"
    )

    run_cmd(cmd, cwd=cwd, log_file=log)


def get_mapped_reads(output_dir: Path, read1: Path, read2: Optional[Path], log: Path) -> None:
    cwd = output_dir / "scaffold-truncated"
    if read2 is None:
        cmd = f"rm -rf tmp && mkdir -p tmp && bash filterbyname.sh in={str(read1)} out=tmp/mapped_reads_1.fastq names=salmon-mapped.sam include=t"
    else:
        cmd = f"rm -rf tmp && mkdir -p tmp && bash filterbyname.sh in={str(read1)} in2={str(read2)} out=tmp/mapped_reads_1.fastq out2=tmp/mapped_reads_2.fastq names=salmon-mapped.sam include=t"

    run_cmd(cmd, cwd=cwd, log_file=log)


def run_spades(output_dir: Path, paired: bool, spades_kmer: str, log: Path) -> None:
    cwd = output_dir / "scaffold-truncated" / "tmp"
    cwd.mkdir(parents=True, exist_ok=True)

    if spades_kmer == DEFAULT_SPADES_KMER:
        kmer_arg = ""
    else:
        kmer_arg = f"-k {spades_kmer}"

    if paired:
        reads_arg = "-1 mapped_reads_1.fastq -2 mapped_reads_2.fastq"
    else:
        reads_arg = "-s mapped_reads_1.fastq"

    cmd = (
        f"spades.py -o spades-res {reads_arg} --trusted-contigs ../scaffold.fasta {kmer_arg} --only-assembler && "
        "cd spades-res && samtools faidx scaffolds.fasta"
    )

    run_cmd(cmd, cwd=cwd, log_file=log)


def get_scaffold_from_scaffolds(output_dir: Path, log: Path) -> int:
    """Return max scaffold length found in spades-res scaffolds.fasta; if none, read scaffold.fasta.fai in root.
    Also write out selected scaffold to scaffold-truncated if applicable (preserve original behavior).
    """
    spades_scaffolds = output_dir / "scaffold-truncated" / "tmp" / "spades-res" / "scaffolds.fasta"
    if spades_scaffolds.exists():
        fai = spades_scaffolds.with_suffix(spades_scaffolds.suffix + ".fai")
        if not fai.exists():
            raise RunError(f"Missing faidx for {spades_scaffolds}")
        first = safe_read_first_line(fai)
        cols = first.split("\t")
        max_id = cols[0].strip()
        max_len = int(cols[1].strip())

        # extract that scaffold to scaffold.fasta
        cmd = f"bash filterbyname.sh in=scaffolds.fasta out=scaffold.fasta names={max_id} include=t && samtools faidx scaffold.fasta"
        run_cmd(cmd, cwd=spades_scaffolds.parent, log_file=log)
        return max_len

    # fallback to scaffold-truncated/scaffold.fasta.fai
    fallback = output_dir / "scaffold-truncated" / "scaffold.fasta.fai"
    if not fallback.exists():
        raise RunError(f"No scaffolds.fasta found in spades-res and missing fallback {fallback}")
    first = safe_read_first_line(fallback)
    cols = first.split("\t")
    return int(cols[1].strip())


def get_scaffold_len_from_truncated_extend(output_dir: Path) -> int:
    a = output_dir / "scaffold-truncated" / "scaffold.fasta"
    if a.exists() and a.is_file():
        fai = a.with_suffix(a.suffix + ".fai")
        if not fai.exists():
            raise RunError(f"Missing FAI: {fai}")
        first = safe_read_first_line(fai)
        return int(first.split("\t")[1].strip())

    b = output_dir / "scaffold.fasta.fai"
    if not b.exists():
        raise RunError(f"Missing FAI: {b}")
    first = safe_read_first_line(b)
    return int(first.split("\t")[1].strip())


def update_current_scaffold(output_dir: Path, log: Path) -> None:
    sfile = output_dir / "scaffold-truncated" / "scaffold.fasta"
    if sfile.exists():
        cmd = (
            f"cd {str(output_dir)} && rm -f scaffold.fasta* && cd scaffold-truncated && cp scaffold.fasta {str(output_dir)} && cp scaffold.fasta.fai {str(output_dir)}"
        )
        run_cmd(cmd, cwd=output_dir, log_file=log)


def create_bed_for_truncated_scaffold(output_dir: Path, truncated_len: int) -> None:
    fai = output_dir / "scaffold.fasta.fai"
    if not fai.exists():
        raise RunError(f"Missing FAI: {fai}")
    first = safe_read_first_line(fai)
    cols = first.split("\t")
    scaffold_id = cols[0].strip()
    scaffold_length = int(cols[1].strip())

    bed = output_dir / "scaffold-truncated.bed"
    with bed.open("w") as bw:
        if scaffold_length >= truncated_len:
            bw.write(f"{scaffold_id}\t{truncated_len}\t{scaffold_length - truncated_len}\n")


def grow_scaffold_with_assembly(output_dir: Path, read1: Path, read2: Optional[Path], spades_kmer: str, avg_read_len: float, log: Path) -> None:
    iteration = 1
    extend_contig = True
    prev_length = 0

    while extend_contig:
        current_length = get_scaffold_from_scaffolds(output_dir, log)
        if current_length > 300000:
            extend_contig = False
            break

        if current_length > prev_length:
            prev_length = current_length
            # create bed, run alignment, extract mapped reads, run spades
            create_bed(output_dir, avg_read_len, log)
            run_alignment(output_dir, read1, read2, SALMON_READ_FRACTION, log)
            get_mapped_reads(output_dir, read1, read2, log)
            paired = read2 is not None
            run_spades(output_dir, paired, spades_kmer, log)

            # check whether spades produced new scaffolds
            spades_scaffolds = output_dir / "scaffold-truncated" / "tmp" / "spades-res" / "scaffolds.fasta"
            if spades_scaffolds.exists():
                iteration += 1
                extend_contig = True
            else:
                extend_contig = False
        else:
            extend_contig = False

        if extend_contig and iteration > 100:
            extend_contig = False


def get_truncated_scaffold_and_extend(output_dir: Path, current_length: int, read1: Path, read2: Optional[Path], spades_kmer: str, avg_read_len: float, log: Path) -> None:
    # tries truncations in increasing size similar to original
    truncations = [300, 500, 700, 1000, 1300, 1500, 1700, 2000]

    for trunc in truncations:
        create_bed_for_truncated_scaffold(output_dir, trunc)
        # recreate scaffold-truncated from root scaffold.fasta
        cmd = (
            f"cd {str(output_dir)} && rm -rf scaffold-truncated && mkdir -p scaffold-truncated && "
            f"bedtools getfasta -fi scaffold.fasta -bed scaffold-truncated.bed -fo scaffold-truncated/scaffold.fasta && "
            f"cd scaffold-truncated && samtools faidx scaffold.fasta"
        )
        run_cmd(cmd, cwd=output_dir, log_file=log)

        grow_scaffold_with_assembly(output_dir, read1, read2, spades_kmer, avg_read_len, log)
        length_from_growing = get_scaffold_len_from_truncated_extend(output_dir)
        if length_from_growing > current_length:
            return
    # if none increased, return (no change)


def check_circularity(output_dir: Path, avg_read_len: float, log: Path) -> bool:
    # remove existing blastn-* directories
    for p in output_dir.glob("blastn-*"):
        if p.is_dir() or p.is_file():
            if p.is_dir():
                run_cmd(f"rm -rf {str(p)}", cwd=output_dir, log_file=log)
            else:
                p.unlink()

    # create query and subject bed files based on scaffold.fasta.fai
    fai = output_dir / "scaffold.fasta.fai"
    if not fai.exists():
        raise RunError(f"Missing FAI: {fai}")

    first = safe_read_first_line(fai)
    cols = first.split("\t")
    scaffold_id = cols[0].strip()
    scaffold_len = int(cols[1].strip())
    int_read_len = int(avg_read_len)

    if scaffold_len <= int_read_len * 2:
        return False

    subject_bed = output_dir / "blastn-subject-1stround.bed"
    query_bed = output_dir / "blastn-query-1stround.bed"
    with subject_bed.open("w") as s, query_bed.open("w") as q:
        s.write(f"{scaffold_id}\t0\t{scaffold_len - (int_read_len*2)}\n")
        q.write(f"{scaffold_id}\t{scaffold_len - (int_read_len*2)}\t{scaffold_len-int_read_len}\n")

    cmd = (
        f"cd {str(output_dir)} && "
        "bedtools getfasta -fi scaffold.fasta -bed blastn-subject-1stround.bed -fo blastn-subject-1stround.fasta && "
        "bedtools getfasta -fi scaffold.fasta -bed blastn-query-1stround.bed -fo blastn-query-1stround.fasta && "
        "samtools faidx blastn-subject-1stround.fasta && samtools faidx blastn-query-1stround.fasta && "
        "makeblastdb -in blastn-subject-1stround.fasta -dbtype nucl && "
        "blastn -query blastn-query-1stround.fasta -db blastn-subject-1stround.fasta -num_threads 16 -outfmt '7' -out blastn-res-1stround.txt"
    )

    run_cmd(cmd, cwd=output_dir, log_file=log)

    # parse blastn output
    res_file = output_dir / "blastn-res-1stround.txt"
    if not res_file.exists():
        return False

    is_circular = False
    subject_start = 0
    query_start = 0
    with res_file.open("r") as br, (output_dir / "circularity-output-log.txt").open("a") as bc, (output_dir / "output-log.txt").open("a") as bo:
        for line in br:
            if not line:
                break
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            percent_iden = float(fields[2].strip())
            alignment_len = int(fields[3].strip())
            subject_start = int(fields[8].strip())
            query_start = int(fields[6].strip())
            min_alignment_len = int(avg_read_len * 0.95)

            if (percent_iden >= MIN_IDENTITY_CIRCULAR) and (alignment_len >= min_alignment_len):
                bc.write(f"Scaffold seems circular. Scaffold position {fields[6]} to {fields[7]} mapped to position {fields[8]} to {fields[9]} with {percent_iden}% identity and {alignment_len} alignment length.\n")
                bo.write(f"Scaffold seems circular. Scaffold position {fields[6]} to {fields[7]} mapped to position {fields[8]} to {fields[9]} with {percent_iden}% identity and {alignment_len} alignment length.\n")
                is_circular = True
                break
            else:
                bc.write(f"Scaffold does not seem to be circular. Scaffold position {fields[6]} to {fields[7]} mapped to position {fields[8]} to {fields[9]} with {percent_iden}% identity and {alignment_len} alignment length.\n")
                bo.write(f"Scaffold does not seem to be circular. Scaffold position {fields[6]} to {fields[7]} mapped to position {fields[8]} to {fields[9]} with {percent_iden}% identity and {alignment_len} alignment length.\n")

    if not is_circular:
        return False

    # second round if circular
    fai = output_dir / "scaffold.fasta.fai"
    first = safe_read_first_line(fai)
    cols = first.split("\t")
    scaffold_len = int(cols[1].strip())
    int_read_len = int(avg_read_len)

    br = subject_start
    # writing 2nd round beds
    with (output_dir / "blastn-subject-2ndround.bed").open("w") as b1, (output_dir / "blastn-query-2ndround.bed").open("w") as b2:
        if scaffold_len > (int_read_len * 2):
            b1.write(f"{cols[0]}\t0\t{subject_start}\n")
            b2.write(f"{cols[0]}\t{scaffold_len - (int_read_len*2) - subject_start}\t{scaffold_len-(int_read_len*2)}\n")

    cmd2 = (
        f"cd {str(output_dir)} && "
        "bedtools getfasta -fi scaffold.fasta -bed blastn-subject-2ndround.bed -fo blastn-subject-2ndround.fasta && "
        "bedtools getfasta -fi scaffold.fasta -bed blastn-query-2ndround.bed -fo blastn-query-2ndround.fasta && "
        "samtools faidx blastn-subject-2ndround.fasta && samtools faidx blastn-query-2ndround.fasta && "
        "makeblastdb -in blastn-subject-2ndround.fasta -dbtype nucl && "
        "blastn -query blastn-query-2ndround.fasta -db blastn-subject-2ndround.fasta -num_threads 16 -outfmt '7' -out blastn-res-2ndround.txt"
    )
    run_cmd(cmd2, cwd=output_dir, log_file=log)

    return True


def run_alignment_get_coverage(output_dir: Path, read1: Path, read2: Optional[Path], log: Path) -> None:
    # cleanup old indices
    for p in output_dir.glob("bowtie2*"):
        if p.exists():
            run_cmd(f"rm -rf {str(p)}", cwd=output_dir, log_file=log)
    for p in output_dir.glob("samtools*"):
        if p.exists():
            run_cmd(f"rm -rf {str(p)}", cwd=output_dir, log_file=log)

    if read2 is None:
        cmd = (
            f"cd {str(output_dir)} && bowtie2-build scaffold.fasta bowtie2-index && "
            f"bowtie2 -x bowtie2-index -U {str(read1)} | samtools view -bS - | samtools view -h -F 0x04 -b - | samtools sort - -o bowtie2-mapped.bam && "
            "samtools depth -a bowtie2-mapped.bam > samtools-coverage.txt"
        )
    else:
        cmd = (
            f"cd {str(output_dir)} && bowtie2-build scaffold.fasta bowtie2-index && "
            f"bowtie2 -x bowtie2-index -1 {str(read1)} -2 {str(read2)} | samtools view -bS - | samtools view -h -F 0x04 -b - | samtools sort - -o bowtie2-mapped.bam && "
            "samtools depth -a bowtie2-mapped.bam > samtools-coverage.txt"
        )

    run_cmd(cmd, cwd=output_dir, log_file=log)


def percentile(values: List[int], percentiles: List[float]) -> List[int]:
    if not values:
        return [0 for _ in percentiles]
    temp = sorted(values)
    results = []
    for p in percentiles:
        idx = int(math.ceil((p / 100.0) * len(temp))) - 1
        idx = max(0, min(idx, len(temp) - 1))
        results.append(temp[idx])
    return results


def read_coverage_get_percentile(output_dir: Path) -> Tuple[int, int]:
    cov_file = output_dir / "samtools-coverage.txt"
    if not cov_file.exists():
        raise RunError(f"Coverage file not found: {cov_file}")
    coverages: List[int] = []
    with cov_file.open("r") as br:
        for line in br:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            cov = int(cols[2])
            if cov != 0:
                coverages.append(cov)

    p15, p85 = percentile(coverages, [15.0, 85.0])
    return p15, p85


def write_coverage_quantile(output_dir: Path, quantile15Percent: int, quantile85Percent: int) -> bool:
    suspicious_starts: List[int] = []
    suspicious_ends: List[int] = []

    fai = output_dir / "scaffold.fasta.fai"
    if not fai.exists():
        raise RunError(f"Missing FAI: {fai}")
    first = safe_read_first_line(fai).strip()
    cols = first.split("\t")
    scaffold_id = cols[0].strip()
    scaffold_length = int(cols[1].strip())

    cov_file = output_dir / "samtools-coverage.txt"
    if not cov_file.exists():
        raise RunError(f"Missing coverage file: {cov_file}")

    low_high = False
    start_base = 0
    with cov_file.open("r") as br:
        for line in br:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            coverage = int(cols[2])
            current_base = int(cols[1])
            if (coverage >= quantile15Percent) and (coverage <= quantile85Percent):
                if low_high:
                    if (current_base - start_base) > MIN_SUSPICIOUS_LEN:
                        suspicious_starts.append(start_base)
                        suspicious_ends.append(current_base - 1)
                low_high = False
            else:
                if (start_base == 0) or (not low_high):
                    start_base = current_base
                low_high = True

    if not suspicious_starts:
        with (output_dir / "suspicious-regions.log").open("w") as bw:
            bw.write(f"{scaffold_id} has no suspicious region\n")
        return False

    # find longest non-suspicious region
    max_start_base = 1
    max_length = 0
    for i, s in enumerate(suspicious_starts):
        if i == 0:
            cur_len = s - 1
            if cur_len > max_length:
                max_length = cur_len
                max_start_base = 1
        else:
            cur_len = s - suspicious_ends[i - 1]
            if cur_len > max_length:
                max_length = cur_len
                max_start_base = suspicious_ends[i - 1] + 1

    last_gap = scaffold_length - suspicious_ends[-1]
    if last_gap > max_length:
        max_length = last_gap
        max_start_base = suspicious_ends[-1] + 1

    with (output_dir / "longest-non-suspicious.bed").open("w") as bw:
        bw.write(f"{scaffold_id}\t{max_start_base}\t{max_start_base + max_length - 1}\n")

    # extract sequence
    cmd = (
        f"cd {str(output_dir)} && rm -f longest-non-suspicious.fasta* && "
        "bedtools getfasta -fi scaffold.fasta -bed longest-non-suspicious.bed -fo longest-non-suspicious.fasta && "
        "samtools faidx longest-non-suspicious.fasta"
    )
    run_cmd(cmd, cwd=output_dir, log_file=output_dir / "run.log")

    return True


def update_scaffold_with_longest(output_dir: Path, log: Path) -> None:
    scaffold_file = output_dir / "scaffold.fasta"
    if scaffold_file.exists():
        cmd = f"cd {str(output_dir)} && rm -f scaffold.fasta* && cp longest-non-suspicious.fasta scaffold.fasta && samtools faidx scaffold.fasta"
        run_cmd(cmd, cwd=output_dir, log_file=log)


def copy_scaffold_for_assembly(output_dir: Path, log: Path) -> None:
    scaffold_file = output_dir / "scaffold.fasta"
    if scaffold_file.exists():
        cmd = (
            f"cd {str(output_dir)} && rm -rf scaffold-truncated && mkdir -p scaffold-truncated && cp scaffold.fasta scaffold-truncated/scaffold.fasta && cd scaffold-truncated && samtools faidx scaffold.fasta"
        )
        run_cmd(cmd, cwd=output_dir, log_file=log)


def check_coverage(output_dir: Path, read1: Path, read2: Optional[Path], log: Path) -> bool:
    run_alignment_get_coverage(output_dir, read1, read2, log)
    p15, p85 = read_coverage_get_percentile(output_dir)
    has_suspicious = write_coverage_quantile(output_dir, p15, p85)

    with (output_dir / "output-log.txt").open("a") as out:
        if has_suspicious:
            update_scaffold_with_longest(output_dir, log)
            out.write("Scaffold got updated for suspicious region.\n")
            copy_scaffold_for_assembly(output_dir, log)
            grow_scaffold_with_assembly(output_dir, read1, read2, DEFAULT_SPADES_KMER, float(0), log)
        else:
            out.write("Scaffold didn't get updated as there is no suspicious region.\n")

    return has_suspicious


def extend_one_scaffold(scaffold: Path, output_dir: Path, read1: Path, read2: Optional[Path], spades_kmer: str, log: Path) -> None:
    # copy scaffold into output and index
    cmd = f"cd {str(output_dir)} && cp {str(scaffold)} scaffold.fasta && samtools faidx scaffold.fasta"
    run_cmd(cmd, cwd=output_dir, log_file=log)

    extend_contig = True
    prev_length = 0
    iteration = 1

    while extend_contig:
        current_length = get_scaffold_len_from_truncated_extend(output_dir)
        if current_length > 300000:
            break
        if current_length > prev_length:
            prev_length = current_length
            update_current_scaffold(output_dir, log)
            if check_circularity(output_dir, 0, log):  # avg_read_len not available here; original called getReadLen earlier
                break
            else:
                need_update = check_coverage(output_dir, read1, read2, log)
                if need_update:
                    new_length = get_scaffold_len_from_truncated_extend(output_dir)
                    update_current_scaffold(output_dir, log)
                    prev_length = new_length
                    get_truncated_scaffold_and_extend(output_dir, new_length, read1, read2, spades_kmer, 0, log)
                else:
                    get_truncated_scaffold_and_extend(output_dir, current_length, read1, read2, spades_kmer, 0, log)
            iteration += 1
        else:
            break

        if iteration > 200:
            break


# --- CLI ---------------------------------------------------------------------

@click.group()
def cli():
    """An integrated pipeline for error-correction, extension, and annotation of viral contigs."""
    pass


@cli.command()
@click.version_option(__version__)
@click.option("-i", "--input", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("-r1", "--read1", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("-r2", "--read2", required=False, type=click.Path(exists=True, dir_okay=False))
@click.option("-o", "--output-dir", required=True, type=click.Path())
@click.option("--spades-kmer", default=DEFAULT_SPADES_KMER, help="K-mer length for spades or 'default'.")
def main(input: str, read1: str, read2: Optional[str], output_dir: str, spades_kmer: str):
    input_path = Path(input)
    read1_path = Path(read1)
    read2_path = Path(read2) if read2 else None
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    log = out_dir / "run.log"
    # start fresh log
    if log.exists():
        log.unlink()
    log.write_text("# virseqimprover run log\n")

    try:
        # estimate avg read length (stop on failure)
        avg_read_len = get_read_length(read1_path, out_dir, log)

        print("Started growing scaffold.")
        extend_one_scaffold(input_path, out_dir, read1_path, read2_path, spades_kmer, log)
        print("Finished growing scaffold.")

        # final polishing (bowtie2 + pilon)
        if read2_path is None:
            cmd = (
                f"cd {str(out_dir)} && bowtie2-build scaffold.fasta bowtie2-index && "
                f"bowtie2 -x bowtie2-index -U {str(read1_path)} | samtools view -bS - | samtools view -h -F 0x04 -b - | samtools sort - -o bowtie2-mapped.bam && "
                "samtools index bowtie2-mapped.bam && "
                "pilon --genome scaffold.fasta --bam bowtie2-mapped.bam --output pilon_out --outdir ./"
            )
        else:
            cmd = (
                f"cd {str(out_dir)} && bowtie2-build scaffold.fasta bowtie2-index && "
                f"bowtie2 -x bowtie2-index -1 {str(read1_path)} -2 {str(read2_path)} | samtools view -bS - | samtools view -h -F 0x04 -b - | samtools sort - -o bowtie2-mapped.bam && "
                "samtools index bowtie2-mapped.bam && "
                "pilon --genome scaffold.fasta --bam bowtie2-mapped.bam --output pilon_out --outdir ./"
            )

        run_cmd(cmd, cwd=out_dir, log_file=log)
        print("Program finished.")

    except RunError as e:
        print(f"ERROR: {e}")
        print(f"See run.log for details: {str(log)}")
        sys.exit(1)


if __name__ == "__main__":
    cli()
