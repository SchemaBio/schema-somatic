"""Genotype curated chemotherapy SNVs in one indexed BAM (1-based catalog coordinates)."""
import argparse
import csv
from pathlib import Path
from tempfile import TemporaryDirectory

import pysam
import pysam.bcftools as bcftools


def inspect_site(bam, fasta, site, build, work, min_depth):
    chrom = site["chrom"]
    pos = int(site[build])
    ref, alt = site["ref"], site["alt"]
    evidence = dict(status="no_coverage", depth=0, ref_reads=0, alt_reads=0, gt="./.")
    chrom = chrom if chrom in fasta.references else "chr" + chrom
    if chrom not in fasta.references or fasta.fetch(chrom, pos - 1, pos).upper() != ref:
        raise ValueError(f"{build} reference mismatch at {site['rsid']} ({chrom}:{pos})")
    if chrom not in bam.references:
        return dict(evidence, status="missing_bam_contig")

    pile, calls = str(work / "pile.vcf"), str(work / "calls.vcf")
    bcftools.mpileup("-f", fasta.filename.decode(), "-r", f"{chrom}:{pos}-{pos}",
                     "-q", "20", "-Q", "20", "-d", "100000", "-I",
                     "-a", "FORMAT/DP,FORMAT/AD", "-Ov", "-o", pile,
                     bam.filename.decode(), catch_stdout=False)
    with pysam.VariantFile(pile) as records:
        record = next(iter(records), None)
        if record is None:
            return evidence
        sample = next(iter(record.samples.values()))
        counts = dict(zip(record.alleles, sample["AD"]))
        evidence.update(depth=sample["DP"], ref_reads=counts[ref], alt_reads=counts.get(alt, 0))

    if evidence["depth"] < min_depth:
        return dict(evidence, status="low_depth")
    bcftools.call("-m", "-A", "-f", "GQ", "-Ov", "-o", calls, pile, catch_stdout=False)
    with pysam.VariantFile(calls) as records:
        record = next(iter(records))
        sample = next(iter(record.samples.values()))
        alleles = sample.alleles
        if not alleles or any(a not in (ref, alt) for a in alleles):
            return dict(evidence, status="other_or_uncertain_genotype")
        copies = alleles.count(alt)
        if copies and (evidence["alt_reads"] < 3 or (sample.get("GQ") or 0) < 20):
            return dict(evidence, status="uncertain_genotype")
        if copies == 0 and evidence["ref_reads"] < 0.95 * evidence["depth"]:
            return dict(evidence, status="mixed_evidence")
        return dict(evidence, status=("reference", "one_copy", "two_copies")[copies],
                    gt="/".join("1" if a == alt else "0" for a in alleles))


def interpret(evidence, site, role):
    status = evidence["status"]
    if role == "tumor":
        if status in ("one_copy", "two_copies") or (status == "mixed_evidence" and evidence["alt_reads"] >= 3):
            return "肿瘤检出候选等位基因；胚系来源待确认", "待胚系确认", "待胚系确认"
    if status in ("one_copy", "two_copies"):
        toxicity = site["toxicity_one_copy" if status == "one_copy" else "toxicity_two_copies"]
        return "位点级关联；须结合适用人群解读", site["efficacy_if_carrier"], toxicity
    if status == "reference":
        return "本位点未检出指定等位基因", "未判定", "未判定"
    return "证据不足", "未判定", "未判定"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("bam", "fasta", "catalog", "output", "sample"):
        parser.add_argument("--" + key, required=True)
    parser.add_argument("--assembly", choices=("GRCh37", "GRCh38"), required=True)
    parser.add_argument("--role", choices=("tumor", "normal"), required=True)
    parser.add_argument("--min-depth", type=int, default=20)
    args = parser.parse_args()
    with pysam.AlignmentFile(args.bam) as bam, pysam.FastaFile(args.fasta) as fasta, \
            open(args.catalog, encoding="utf-8") as catalog, TemporaryDirectory() as work, \
            open(args.output, "w", newline="", encoding="utf-8") as out:
        sites = csv.DictReader(catalog, delimiter="\t")
        fields = ["sample", "source_sample", "assembly", "position", "status", "gt", "depth",
                  "ref_reads", "alt_reads", "interpretation", "efficacy", "toxicity"]
        report = csv.DictWriter(out, fields + sites.fieldnames, delimiter="\t", lineterminator="\n")
        report.writeheader()
        for site in sites:
            evidence = inspect_site(bam, fasta, site, args.assembly, Path(work), args.min_depth)
            note, efficacy, toxicity = interpret(evidence, site, args.role)
            if args.role == "tumor":
                evidence["gt"] = "./."
                if evidence["status"] in ("one_copy", "two_copies"):
                    evidence["status"] = "tumor_alt_detected"
            report.writerow(dict(site, **evidence, sample=args.sample, source_sample=args.role,
                                 assembly=args.assembly, position=site[args.assembly],
                                 interpretation=note, efficacy=efficacy, toxicity=toxicity))


if __name__ == "__main__":
    main()
