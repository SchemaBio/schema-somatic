#!/usr/bin/env python3
"""
参考基因组数据库下载脚本
用于下载GRCh37和GRCh38参考基因组及相关数据库文件
"""

import os
import argparse
import urllib.request
import tarfile
import gzip
import shutil
from pathlib import Path


# 参考基因组下载链接
REFGENOME_SOURCES = {
    'GRCh38': {
        'fasta': 'https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz',
        'gff': 'https://ftp.ensembl.org/pub/release-98/gff3/homo_sapiens/Homo_sapiens.GRCh38.98.gff3.gz',
        'vep': 'https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_merged_vep_115_GRCh38.tar.gz',
        'evoscore2': 'https://huggingface.co/datasets/pzweuj/EVOScore-2/blob/main/hg38_VESM_3B_scores.vcf.gz',
        'alphamissense': '',
        'gnomad': '',
        'clnvar': '',
        'exon_bed': '',
        'civic': '',
        'pangolin': '',
        'cytoband': '',
        'mskcc_hotspot': ''
    },
    'GRCh37': {
        'fasta': 'https://ftp.ensembl.org/pub/grch37/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz',
        'gff': 'https://ftp.ensembl.org/pub/release-75/gff3/homo_sapiens/Homo_sapiens.GRCh37.75.gff3.gz',
        'vep': 'https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache/homo_sapiens_merged_vep_115_GRCh37.tar.gz',
        'evoscore2': 'https://huggingface.co/datasets/pzweuj/EVOScore-2/blob/main/hg19_VESM_3B_scores.vcf.gz',
        'alphamissense': '',
        'gnomad': '',
        'clnvar': '',
        'exon_bed': '',
        'civic': '',
        'pangolin': '',
        'cytoband': '',
        'mskcc_hotspot': ''
    }
}


def download_file(url, output_path, desc='Downloading'):
    """下载文件"""
    print(f"{desc}: {url}")
    print(f"Saving to: {output_path}")

    try:
        urllib.request.urlretrieve(url, output_path)
        print(f"Download complete: {output_path}")
        return True
    except Exception as e:
        print(f"Error downloading: {e}")
        return False


def decompress_gz(input_file, output_file):
    """解压gz文件"""
    print(f"Decompressing: {input_file}")
    with gzip.open(input_file, 'rb') as f_in:
        with open(output_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f"Decompressed to: {output_file}")


def download_refgenome(genome_version, output_dir):
    """下载参考基因组"""
    output_path = Path(output_dir) / genome_version
    output_path.mkdir(parents=True, exist_ok=True)

    sources = REFGENOME_SOURCES.get(genome_version)
    if not sources:
        print(f"Error: Unknown genome version: {genome_version}")
        return False

    # 下载fasta
    fasta_gz = output_path / f"{genome_version}.fa.gz"
    fasta = output_path / f"{genome_version}.fa"

    if not fasta.exists():
        if not download_file(sources['fasta'], str(fasta_gz), f"Downloading {genome_version} FASTA"):
            return False
        decompress_gz(str(fasta_gz), str(fasta))
        os.remove(fasta_gz)
    else:
        print(f"FASTA already exists: {fasta}")

    return True


def download_annotation(genome_version, output_dir):
    """下载注释文件"""
    output_path = Path(output_dir) / genome_version

    sources = REFGENOME_SOURCES.get(genome_version)
    if not sources:
        return False

    # 下载GFF3
    gff_gz = output_path / f"{genome_version}.gff3.gz"
    gff = output_path / f"{genome_version}.gff3"

    if not gff.exists():
        if not download_file(sources['gff'], str(gff_gz), f"Downloading {genome_version} GFF3"):
            return False
        decompress_gz(str(gff_gz), str(gff))
        os.remove(gff_gz)
    else:
        print(f"GFF3 already exists: {gff}")

    return True


def download_vep_cache(genome_version, output_dir):
    """下载VEP缓存"""
    output_path = Path(output_dir) / genome_version / 'vep_cache'
    output_path.mkdir(parents=True, exist_ok=True)

    vep_cache_url = REFGENOME_SOURCES.get(genome_version, {}).get('vep')

    if not vep_cache_url:
        print(f"No VEP cache available for {genome_version}")
        return True

    cache_tar = output_path / 'vep_cache.tar.gz'

    if not (output_path / 'homo_sapiens').exists():
        if not download_file(vep_cache_url, str(cache_tar), f"Downloading {genome_version} VEP cache"):
            return False

        print("Extracting VEP cache...")
        with tarfile.open(cache_tar, 'r:gz') as tar:
            tar.extractall(output_path)
        os.remove(cache_tar)
    else:
        print("VEP cache already exists")

    return True


def main():
    parser = argparse.ArgumentParser(description='Download reference genome database')
    parser.add_argument('--genome', '-g', choices=['GRCh37', 'GRCh38', 'both'], default='both',
                        help='Genome version to download')
    parser.add_argument('--output', '-o', default='./reference',
                        help='Output directory')
    parser.add_argument('--skip-vep', action='store_true',
                        help='Skip VEP cache download')

    args = parser.parse_args()

    genome_versions = ['GRCh37', 'GRCh38'] if args.genome == 'both' else [args.genome]

    for genome in genome_versions:
        print(f"\n{'='*50}")
        print(f"Processing {genome}")
        print(f"{'='*50}\n")

        # 下载参考基因组
        if not download_refgenome(genome, args.output):
            print(f"Failed to download {genome} reference genome")
            continue

        # 下载注释文件
        download_annotation(genome, args.output)

        # 下载VEP缓存
        if not args.skip_vep:
            download_vep_cache(genome, args.output)

    print(f"\n{'='*50}")
    print("Reference genome database download complete!")
    print(f"Output directory: {args.output}")
    print(f"{'='*50}")


if __name__ == '__main__':
    main()