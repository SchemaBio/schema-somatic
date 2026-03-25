# Schema-Somatic Pipeline

基于 Nextflow 的肿瘤体细胞突变检测流程。

！流程建立中，尚不可用！

## 分析模块

| 步骤 | 模块 | 功能 | 输出目录 |
|------|------|------|----------|
| 1 | FASTP | 质量控制、去接头 | 01.QC |
| 2 | BWAMEM | 序列比对 | 02.Alignment |
| 3 | INDEXALIGNMENT | 比对文件索引 | 02.Alignment |
| 4.1 | MUTECT2 → WHATSHAP → EGFRHAP → VEP | 体细胞突变检测、单倍型分型、EGFR区域合并、注释 | 03.Mutations, 04.Annotation |
| 4.2 | HOTSPOT → VEP_HOTSPOT | 热点区域变异检测、注释 | 03.Mutations, 04.Annotation |
| 4.3 | DEEPVARIANT → VEP_GERMLINE | 种系变异检测、注释 | 03.Mutations, 04.Annotation |
| 5 | MANTA | 结构变异检测 | 03.Mutations/SV |
| 6 | GRIDSS | 高精度结构变异 | 03.Mutations/SV |
| 7 | CNVKIT | 拷贝数变异分析 | 03.Mutations/CNV |
| 8 | MSISENSORPRO | 微卫星不稳定性分析 | 03.Mutations/MSI |

## 输出目录结构

```
${params.output}/
├── 01.QC/                                      # 质控结果
│   ├── {sample_id}.clean_1.fq.gz               # 过滤后FASTQ
│   ├── {sample_id}.clean_2.fq.gz
│   ├── {sample_id}.fastp.stat.json             # FASTP统计
│   └── {sample_id}.fastp.stat.html
│
├── 02.Alignment/                               # 比对结果
│   ├── {sample_id}.cram                        # CRAM文件
│   └── {sample_id}.cram.crai                   # 索引
│
├── 03.Mutations/                               # 变异检测结果
│   ├── SNV_InDel/                              # SNP/Indel
│   │   ├── {sample_id}.mutect2.filtered.normalized.vcf.gz
│   │   ├── {sample_id}.phased.vcf.gz           # 单倍型分型
│   │   ├── {sample_id}.egfr.merged.vcf.gz      # EGFR区域合并
│   │   └── {sample_id}.hotspot.vcf.gz          # 热点区域
│   │
│   ├── Germline/                               # 种系变异
│   │   └── {sample_id}.deepvariant.vcf.gz
│   │
│   ├── SV/                                     # 结构变异
│   │   ├── {sample_id}.manta.vcf.gz
│   │   └── {sample_id}.gridss.vcf.gz
│   │
│   ├── CNV/                                    # 拷贝数变异
│   │   └── {sample_id}.call.cns
│   │
│   └── MSI/                                    # 微卫星不稳定性
│       └── {sample_id}_msisensor/
│
├── 04.Annotation/                              # 注释结果（最终报告）
│   ├── {sample_id}.egfr.vcf                    # EGFR注释结果
│   ├── {sample_id}.egfr.html
│   ├── {sample_id}.hotspot.vcf                 # 热点注释结果
│   ├── {sample_id}.hotspot.html
│   ├── {sample_id}.germline.vcf                # 种系变异注释结果
│   └── {sample_id}.germline.html
│
└── reports/                                    # 流程报告
    ├── pipeline_report.html
    ├── timeline.html
    └── pipeline_dag.html
```

## 运行参数

### 基本参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--input` | FASTQ输入目录 | null |
| `--output` | 输出目录 | results |
| `--reference` | 参考基因组FASTA | null |
| `--genome_version` | 基因组版本 | GRCh38 |
| `--analysis_mode` | 分析模式 (single/paired) | single |

### 变异检测参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--run_mutect2` | 运行MUTECT2体细胞突变检测 | false |
| `--run_hotspot` | 运行热点区域检测 | false |
| `--hotspot` | 热点区域VCF文件 | null |
| `--run_deepvariant` | 运行DeepVariant种系变异检测 | false |
| `--run_manta` | 运行Manta结构变异检测 | false |
| `--run_gridss` | 运行GRIDSS结构变异检测 | false |
| `--run_cnvkit` | 运行CNVkit分析 | false |
| `--run_msisensorpro` | 运行MSIsensor-pro分析 | false |

### 注释参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--vep_database` | VEP数据库路径 | null |

### CNV/MSI参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--cnv_baseline` | CNV基线文件 | null |
| `--msi_baseline` | MSI基线文件 | null |
| `--user_bed` | 捕获区域BED文件 | null |

## 使用示例

```bash
# 单样本全流程分析
nextflow run main.nf --analysis_mode single \
    --input /path/to/fastq \
    --output results \
    --reference /path/to/reference.fa \
    --run_mutect2 true \
    --run_hotspot true \
    --hotspot /path/to/hotspot.vcf \
    --run_deepvariant true \
    --vep_database /path/to/vep_cache \
    --cnv_baseline /path/to/cnv.reference.cnn \
    --msi_baseline /path/to/msi_baseline
```

## 依赖环境

- Nextflow >= 21.0
- Docker / Singularity / Conda
- 参考基因组 (GRCh37/GRCh38)
- VEP数据库