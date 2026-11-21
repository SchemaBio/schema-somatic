version 1.2

task FixBed {
    input {
        File bed
        Directory ref_dir
        String fasta
    }

    command <<<
        set -euo pipefail
        awk 'BEGIN {FS=OFS="\t"}
            NR==FNR {contigs[$1]=1; next}
            /^#|^track |^browser |^[[:space:]]*$/ {next}
            {
                chrom=$1
                if (!(chrom in contigs)) {
                    sub(/^chr/, "", chrom)
                    if (chrom=="M") chrom="MT"
                    if (!(chrom in contigs)) chrom=(chrom=="MT" ? "chrM" : "chr" chrom)
                }
                if (!(chrom in contigs) || NF<3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3<$2) {
                    print "Invalid BED interval: " $0 > "/dev/stderr"; exit 1
                }
                $1=chrom; if ($3==$2) $3++; print
            }' "~{ref_dir}/~{fasta}.fai" "~{bed}" > fixed.bed
        test -s fixed.bed
    >>>

    output {
        File fixed_bed = "fixed.bed"
    }

    requirements {
        container: "schemabio/mapping:v1.1.0"
    }
}

task QCReport {
    input {
        String prefix
        File fastp_stats
        File xamdst_json
        File gatk_metric
        File hs_metric
        File sry_result
        Int sry_cutoff
        File stats
    }

    command <<<
        python3 - "~{prefix}" "~{fastp_stats}" "~{xamdst_json}" "~{gatk_metric}" "~{hs_metric}" "~{sry_result}" ~{sry_cutoff} "~{stats}" <<'PY'
        import csv
        import json
        import sys
        from pathlib import Path

        sample, fastp, coverage, alignment, capture, sry, cutoff, duplicates = sys.argv[1:]

        def metrics(path):
            lines = iter(Path(path).read_text().splitlines())
            next(line for line in lines if line.startswith("## METRICS CLASS"))
            rows = []
            for line in lines:
                if not line.strip():
                    break
                rows.append(line)
            return list(csv.DictReader(rows, delimiter="\t"))

        raw_coverage = json.loads(Path(coverage).read_text())
        target = raw_coverage["target"]
        insert_size = raw_coverage.get("insert_size", {})
        deduplicated = target.get("coverage_rmdup", {})
        xamdst = {
            "mapped_reads": raw_coverage["total"].get("mapped_reads"),
            "mapped_reads_fraction": raw_coverage["total"].get("mapped_reads_fraction"),
            "insert_size_average": insert_size.get("average"),
            "insert_size_median": insert_size.get("median"),
            "region_length": target.get("region_length"),
            "target_data_fraction_all": target.get("target_data_fraction_all"),
            "average_depth": target.get("average_depth"),
            "average_depth_rmdup": target.get("average_depth_rmdup"),
            "coverage_gt_0x": deduplicated.get("gt_0x"),
            "coverage_gte_30x": deduplicated.get("gte_30x"),
            "coverage_gte_100x": deduplicated.get("gte_100x"),
            "coverage_gte_1000x": deduplicated.get("custom", {}).get("gte_1000x"),
            "coverage_gt_0_2_avg": target.get("coverage", {}).get("gt_0_2_avg"),
        }
        count = int(Path(sry).read_text())
        report = {
            "sample": sample,
            "fastp": json.loads(Path(fastp).read_text())["summary"],
            "coverage": raw_coverage,
            "xamdst": xamdst,
            "alignment": metrics(alignment),
            "capture": metrics(capture),
            "sry_reads": count,
            "sry_cutoff": int(cutoff),
            "sry_above_cutoff": count >= int(cutoff),
            "duplicate_metrics": Path(duplicates).read_text().strip(),
        }
        Path(f"{sample}.qc.json").write_text(json.dumps(report, indent=2) + "\n")
        PY
    >>>

    output {
        File qc_report = "~{prefix}.qc.json"
    }

    requirements {
        container: "python:3.12-slim"
    }
}
