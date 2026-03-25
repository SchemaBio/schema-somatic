#!/usr/bin/env python3
"""
EGFR MNP 合并脚本
用途：合并 EGFR 19del 和 20ins 区域的相邻变异（MNP）
特性：
  - 只处理 EGFR 19del 和 20ins 区域
  - 跳过 T790M 和 C797S 位点的合并
  - 支持 GRCh37 和 GRCh38 参考基因组
  - 调用 /usr/local/bin/merge_mnp.py 和 bcftools 进行处理
"""

import argparse
import gzip
import os
import subprocess
import sys
import tempfile
from pathlib import Path


# EGFR 关键位点坐标
# exon 坐标来源：
#   hg19: exon19=chr7:55242415-55242513, exon20=chr7:55248986-55249171
#   hg38: exon19=chr7:55174722-55174820, exon20=chr7:55181293-55181478
# T790M (c.2369C>T) 和 C797S (c.2389T>C) 都在 exon 20，不应参与 MNP 合并
EGFR_COORDS = {
    'GRCh38': {
        'chr': 'chr7',
        # exon 19 范围
        'exon19': (55174722, 55174820),
        # 19del 典型缺失区域（exon19 内）
        'exon19_del': (55174771, 55174791),
        # exon 20 范围
        'exon20': (55181293, 55181478),
        # 20ins 典型插入区域（exon20 内）
        'exon20_ins': (55181320, 55181380),
        # T790M: c.2369C>T, p.Thr790Met
        'T790M': 55181378,
        # C797S: c.2389T>C, p.Cys797Ser
        'C797S': 55181398,
    },
    'GRCh37': {
        'chr': 'chr7',
        # exon 19 范围
        'exon19': (55242415, 55242513),
        # 19del 典型缺失区域（exon19 内）
        'exon19_del': (55242464, 55242484),
        # exon 20 范围
        'exon20': (55248986, 55249171),
        # 20ins 典型插入区域（exon20 内）
        'exon20_ins': (55249071, 55249131),
        # T790M: c.2369C>T, p.Thr790Met
        'T790M': 55249071,
        # C797S: c.2389T>C, p.Cys797Ser
        'C797S': 55249091,
    }
}

# 需要跳过合并的位点
SKIP_POSITIONS = ['T790M', 'C797S']

# 外部工具路径
MERGE_MNP_SCRIPT = '/usr/local/bin/merge_mnp.py'


def run_cmd(cmd, description=''):
    """执行 shell 命令"""
    if description:
        print(f"[INFO] {description}")
    print(f"[CMD] {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] Command failed: {result.stderr}")
        sys.exit(1)
    return result.stdout


def parse_vcf(vcf_path):
    """解析 VCF 文件，返回 header 和变异列表"""
    header_lines = []
    variants = []

    open_func = gzip.open if vcf_path.endswith('.gz') else open
    mode = 'rt' if vcf_path.endswith('.gz') else 'r'

    with open_func(vcf_path, mode) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('#'):
                header_lines.append(line)
            else:
                fields = line.split('\t')
                if len(fields) >= 5:
                    variants.append(fields)

    return header_lines, variants


def write_vcf(output_path, header_lines, variants, compress=False):
    """写入 VCF 文件"""
    if compress or output_path.endswith('.gz'):
        import gzip
        with gzip.open(output_path, 'wt') as f:
            for line in header_lines:
                f.write(line + '\n')
            for var in variants:
                f.write('\t'.join(var) + '\n')
    else:
        with open(output_path, 'w') as f:
            for line in header_lines:
                f.write(line + '\n')
            for var in variants:
                f.write('\t'.join(var) + '\n')


def should_skip_position(pos, coords, tolerance=0):
    """检查该位置是否应该跳过（T790M 或 C797S）"""
    for skip_name in SKIP_POSITIONS:
        skip_pos = coords[skip_name]
        if abs(pos - skip_pos) <= tolerance:
            return True, skip_name
    return False, None


def filter_skip_variants(vcf_path, output_path, coords, skip_positions):
    """
    从 VCF 中移除指定位置的变异，输出到新文件
    返回被移除的变异列表
    """
    header_lines, variants = parse_vcf(vcf_path)

    filtered_variants = []
    skipped_variants = []

    for var in variants:
        pos = int(var[1])
        should_skip, skip_name = should_skip_position(pos, coords)
        if should_skip:
            print(f"[INFO] Skipping {skip_name} at position {pos}")
            skipped_variants.append((skip_name, var))
        else:
            filtered_variants.append(var)

    write_vcf(output_path, header_lines, filtered_variants)
    return skipped_variants, header_lines


def merge_region(input_vcf, output_vcf, reference, region_str, chr_name):
    """调用 merge_mnp.py 合并指定区域的变异"""
    # 先用 bcftools 提取区域
    region_vcf = output_vcf.replace('.vcf', '.region.vcf')

    cmd = ['bcftools', 'view', '-r', f"{chr_name}:{region_str}", input_vcf, '-o', region_vcf]
    run_cmd(cmd, f"Extracting region {chr_name}:{region_str}")

    # 调用 merge_mnp.py
    cmd = ['python', MERGE_MNP_SCRIPT, region_vcf, reference, '--out_file', output_vcf]
    run_cmd(cmd, f"Merging variants in region {chr_name}:{region_str}")

    # 清理临时文件
    if os.path.exists(region_vcf):
        os.remove(region_vcf)


def merge_region_with_skip(input_vcf, output_vcf, reference, coords, region_name):
    """
    合并指定区域的变异，但跳过 T790M 和 C797S
    """
    chr_name = coords['chr']
    region = coords[region_name]
    region_str = f"{region[0]}-{region[1]}"

    print(f"[INFO] Processing {region_name}: {chr_name}:{region_str}")

    # 用 bcftools 提取区域
    region_vcf = output_vcf.replace('.vcf', f'.{region_name}.vcf')
    cmd = ['bcftools', 'view', '-r', f"{chr_name}:{region_str}", input_vcf, '-o', region_vcf]
    run_cmd(cmd, f"Extracting {region_name}")

    # 检查区域内的变异，分离需要跳过的
    header_lines, variants = parse_vcf(region_vcf)

    # 分离：跳过的 vs 需要合并的
    to_merge = []
    skipped = []

    for var in variants:
        pos = int(var[1])
        should_skip, skip_name = should_skip_position(pos, coords)
        if should_skip:
            print(f"[INFO] Skipping {skip_name} at position {pos}")
            skipped.append(var)
        else:
            to_merge.append(var)

    if to_merge:
        # 写入需要合并的变异
        merge_input = output_vcf.replace('.vcf', f'.{region_name}.to_merge.vcf')
        write_vcf(merge_input, header_lines, to_merge)

        # 调用 merge_mnp.py
        merge_output = output_vcf.replace('.vcf', f'.{region_name}.merged.vcf')
        cmd = ['python', MERGE_MNP_SCRIPT, merge_input, reference, '--out_file', merge_output]
        run_cmd(cmd, f"Merging variants in {region_name}")

        # 读取合并后的结果
        _, merged_variants = parse_vcf(merge_output)

        # 合并结果和跳过的变异
        final_variants = merged_variants + skipped
        # 按位置排序
        final_variants.sort(key=lambda x: int(x[1]))

        write_vcf(output_vcf, header_lines, final_variants)

        # 清理临时文件
        if os.path.exists(merge_input):
            os.remove(merge_input)
        if os.path.exists(merge_output):
            os.remove(merge_output)
    else:
        # 没有需要合并的变异，直接输出跳过的
        write_vcf(output_vcf, header_lines, skipped)

    # 清理
    if os.path.exists(region_vcf):
        os.remove(region_vcf)

    return output_vcf


def process_vcf(input_vcf, output_vcf, reference, genome_build):
    """主处理函数"""
    coords = EGFR_COORDS.get(genome_build)
    if not coords:
        print(f"[ERROR] Unknown genome build: {genome_build}")
        sys.exit(1)

    chr_name = coords['chr']

    # 创建临时目录
    with tempfile.TemporaryDirectory() as tmpdir:
        # 1. 处理 19del 区域
        exon19_vcf = os.path.join(tmpdir, 'exon19.merged.vcf')
        merge_region_with_skip(input_vcf, exon19_vcf, reference, coords, 'exon19_del')

        # 2. 处理 20ins 区域
        exon20_vcf = os.path.join(tmpdir, 'exon20.merged.vcf')
        merge_region_with_skip(input_vcf, exon20_vcf, reference, coords, 'exon20_ins')

        # 3. 提取其他区域（排除 19del 和 20ins）
        del_region = coords['exon19_del']
        ins_region = coords['exon20_ins']

        # 创建排除区域的 BED 文件
        exclude_bed = os.path.join(tmpdir, 'exclude.bed')
        with open(exclude_bed, 'w') as f:
            f.write(f"{chr_name}\t{del_region[0] - 1}\t{del_region[1]}\n")
            f.write(f"{chr_name}\t{ins_region[0] - 1}\t{ins_region[1]}\n")

        other_vcf = os.path.join(tmpdir, 'other.vcf.gz')
        cmd = ['bcftools', 'view', '-T', f'^{exclude_bed}', input_vcf, '-O', 'z', '-o', other_vcf]
        run_cmd(cmd, "Extracting other regions")

        # 索引
        cmd = ['tabix', '-p', 'vcf', other_vcf]
        run_cmd(cmd, "Indexing other regions VCF")

        # 4. 合并所有区域
        # 准备合并的文件列表
        files_to_concat = []
        if os.path.exists(exon19_vcf) and os.path.getsize(exon19_vcf) > 0:
            files_to_concat.append(exon19_vcf)
        if os.path.exists(exon20_vcf) and os.path.getsize(exon20_vcf) > 0:
            files_to_concat.append(exon20_vcf)
        files_to_concat.append(other_vcf)

        if len(files_to_concat) > 1:
            cmd = ['bcftools', 'concat', '-a'] + files_to_concat + ['-O', 'z', '-o', output_vcf]
            run_cmd(cmd, "Concatenating all regions")
        else:
            # 只有一个文件，直接复制
            import shutil
            shutil.copy(files_to_concat[0], output_vcf)

        # 5. 创建索引
        cmd = ['tabix', '-p', 'vcf', output_vcf]
        run_cmd(cmd, "Creating final index")

    print(f"[INFO] Done. Output: {output_vcf}")


def main():
    parser = argparse.ArgumentParser(
        description='Merge EGFR MNP variants (19del and 20ins regions), skipping T790M and C797S'
    )
    parser.add_argument('input_vcf', help='Input VCF file (can be gzipped)')
    parser.add_argument('reference', help='Reference FASTA file')
    parser.add_argument('--out_file', '-o', required=True, help='Output VCF file (.vcf.gz)')
    parser.add_argument('--genome', '-g', choices=['GRCh37', 'GRCh38'], default='GRCh38',
                        help='Genome build (default: GRCh38)')

    args = parser.parse_args()

    # 检查输入文件
    if not os.path.exists(args.input_vcf):
        print(f"[ERROR] Input VCF not found: {args.input_vcf}")
        sys.exit(1)

    if not os.path.exists(args.reference):
        print(f"[ERROR] Reference FASTA not found: {args.reference}")
        sys.exit(1)

    # 确保输出文件以 .gz 结尾
    if not args.out_file.endswith('.gz'):
        args.out_file += '.gz'

    process_vcf(
        input_vcf=args.input_vcf,
        output_vcf=args.out_file,
        reference=args.reference,
        genome_build=args.genome
    )


if __name__ == '__main__':
    main()