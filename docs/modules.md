# Schema-Somatic 可用模块列表

## 1. 质控模块

### FASTP
- **路径**: `modules/fastp/main.nf`
- **用途**: 原始数据质控和接头去除
- **输入**:
  - `tuple val(sample_id), path(read1), path(read2)` - 样本ID和双端FASTQ
- **输出**:
  - `clean_reads` - 过滤后的FASTQ
  - `json_report` - JSON质控报告
  - `html_report` - HTML质控报告

## 2. UMI处理模块

### EXTRACTUMI (FGBIO)
- **路径**: `modules/fgbio/main.nf`
- **用途**: 提取分子标签
- **输入**:
  - `tuple val(sample_id), path(read1), path(read2)` - FASTQ文件
  - `val read_structure` - UMI结构（如"5M5S+T +T"）
  - `path fasta` - 参考基因组
  - `path fasta_index` - 参考基因组索引
- **输出**:
  - `unmapped_bam` - 包含UMI的未比对BAM

## 3. 比对模块

### BWAMEM / BWAMEM2
- **路径**: `modules/bwamem/main.nf`
- **用途**: 序列比对
- **输入**:
  - `tuple val(sample_id), path(reads)` - FASTQ文件
  - `path fasta` - 参考基因组
  - `path fasta_index` - 参考基因组索引
  - `val output_format` - 输出格式('bam'或'cram')
  - `val rgid` - Read Group ID
- **输出**:
  - `cram/bam` - 比对文件
  - `crai/bai` - 索引文件

## 4. Samtools工具

### INDEXALIGNMENT
- **路径**: `modules/samtools/main.nf`
- **用途**: 创建BAM/CRAM索引
- **输入**:
  - `tuple val(sample_id), path(alignment_file)` - 比对文件
  - `val publish_dir` - 发布目录
- **输出**:
  - `bai/crai` - 索引文件

### BAM2FASTQ
- **路径**: `modules/samtools/main.nf`
- **用途**: BAM转FASTQ
- **输入**:
  - `tuple val(sample_id), path(bam)` - BAM文件
  - `val publish_dir` - 发布目录
- **输出**:
  - `reads` - FASTQ文件

## 5. GATK变异检测模块

### MARKDUPLICATES
- **路径**: `modules/gatk/main.nf`
- **用途**: 标记PCR重复
- **输入**:
  - `tuple val(sample_id), path(cram), path(crai)` - CRAM文件
  - `path fasta` - 参考基因组
- **输出**:
  - `cram` - 去重后的CRAM
  - `metrics` - 去重指标

### MARKDUPLICATESUMI
- **路径**: `modules/gatk/main.nf`
- **用途**: 基于UMI标记重复
- **输入**: 同MARKDUPLICATES
- **输出**:
  - `cram` - 去重后的CRAM
  - `metrics` - 去重指标
  - `umi_metrics` - UMI指标

### MUTECT2
- **路径**: `modules/gatk/main.nf`
- **用途**: 体细胞突变检测
- **输入**:
  - `tuple val(sample_id), path(cram), path(crai)` - 肿瘤样本
  - `path fasta` - 参考基因组
  - `path fasta_index` - 参考基因组索引
  - `val normal_cram` - 正常样本（可选）
  - `val normal_crai` - 正常样本索引（可选）
  - `val gnomad` - germline resource（可选）
  - `val pon` - panel of normals（可选）
  - `val bed` - 捕获区域（可选）
- **输出**:
  - `vcf` - 变异VCF
  - `tbi` - VCF索引
  - `left_vcf` - 标准化后的VCF
  - `left_tbi` - 标准化VCF索引

## 6. 单倍型分析模块

### WHATSHAP
- **路径**: `modules/whatshap/main.nf`
- **用途**: 单倍型分型
- **输入**:
  - `tuple val(sample_id), path(vcf), path(vcf_index), path(bam), path(bai)` - VCF和BAM
  - `path fasta` - 参考基因组
  - `path fasta_index` - 参考基因组索引
- **输出**:
  - `phased_vcf` - 分型后的VCF

### EGFRHAP
- **路径**: `modules/whatshap/main.nf`
- **用途**: EGFR区域MNP合并
- **输入**:
  - `tuple val(sample_id), path(vcf), path(vcf_index)` - 分型后的VCF
  - `path fasta` - 参考基因组
  - `path fasta_index` - 参考基因组索引
  - `val genome_build` - 基因组版本('GRCh37'或'GRCh38')
- **输出**:
  - `merged_vcf` - 合并后的VCF

## 7. 种系变异检测模块

### DEEPVARIANT
- **路径**: `modules/deepvariant/main.nf`
- **用途**: 深度学习变异检测
- **输入**:
  - `tuple val(sample_id), path(cram), path(crai)` - CRAM文件
  - `path fasta` - 参考基因组
  - `path fasta_index` - 参考基因组索引
  - `val model_type` - 模型类型('WGS'/'WES'/'PACBIO')
- **输出**:
  - `vcf` - 变异VCF
  - `gvcf` - gVCF文件

## 8. 变异注释模块

### VEP
- **路径**: `modules/vep/main.nf`
- **用途**: 变异功能注释
- **输入**:
  - `tuple val(sample_id), path(vcf), path(vcf_index)` - VCF文件
  - `path fasta` - 参考基因组
  - `path fasta_index` - 参考基因组索引
  - `val vep_database` - VEP数据库路径
  - `val genome_build` - 基因组版本('GRCh37'或'GRCh38')
- **输出**:
  - `annotated_vcf` - 注释后的VCF
  - `html_report` - HTML报告

## 9. 结构变异检测模块

### MANTA
- **路径**: `modules/manta/main.nf`
- **用途**: 结构变异检测
- **输入**:
  - `tuple val(sample_id), path(alignment), path(index)` - BAM/CRAM
  - `path fasta` - 参考基因组
  - `path fasta_index` - 参考基因组索引
- **输出**:
  - `vcf` - SV VCF
  - `candidate_vcf` - 候选SV

### GRIDSS
- **路径**: `modules/gridss/main.nf`
- **用途**: 高精度断点检测
- **输入**:
  - `tuple val(sample_id), path(alignment), path(index)` - BAM/CRAM
  - `path fasta` - 参考基因组
  - `path fasta_index` - 参考基因组索引
- **输出**:
  - `vcf` - SV VCF
  - `assembly_bam` - Assembly BAM

## 10. 拷贝数变异检测模块

### CNVKIT
- **路径**: `modules/cnvkit/main.nf`
- **用途**: CNV检测
- **输入**:
  - `tuple val(sample_id), path(alignment), path(index)` - BAM/CRAM
  - `path fasta` - 参考基因组
  - `path fasta_index` - 参考基因组索引
  - `val genome_version` - 基因组版本('GRCh38'或'GRCh37')
  - `path user_bed` - 用户自定义BED（可选，与预置target取交集）
  - `path baseline` - 基线参考
- **预置文件**:
  - `assets/cnvkit/GRCh38_target.bed`
  - `assets/cnvkit/GRCh38_antitarget.bed`
  - `assets/cnvkit/GRCh37_target.bed`
  - `assets/cnvkit/GRCh37_antitarget.bed`
- **输出**:
  - `call_result` - CNV结果
  - `segment_result` - 分段结果

### CNVKITBASELINE
- **路径**: `modules/cnvkit/main.nf`
- **用途**: 多样本创建CNV基线
- **输入**:
  - `path normal_bams` - 多个正常样本BAM
  - `path normal_bais` - 索引文件
  - `path fasta` - 参考基因组
  - `path fasta_index` - 参考基因组索引
  - `val genome_version` - 基因组版本
  - `path user_bed` - 用户自定义BED（可选）
  - `val prefix` - 输出前缀
- **输出**:
  - `baseline` - 基线参考文件

### CNVKITNORMAL
- **路径**: `modules/cnvkit/main.nf`
- **用途**: 单样本创建CNV基线
- **输入**: 同CNVKITBASELINE（单个样本）
- **输出**:
  - `baseline` - 基线参考文件

### CNVKITSEX
- **路径**: `modules/cnvkit/main.nf`
- **用途**: 性别判断
- **输入**:
  - `path cnr` - CNVkit比率文件
  - `val sample_id` - 样本ID
- **输出**:
  - `cnvkit_sex` - 性别判断结果

## 11. MSI检测模块

### MSISENSORPRO_SCAN
- **路径**: `modules/msisensorpro/main.nf`
- **用途**: 扫描参考基因组获取MSI位点
- **输入**:
  - `path fasta` - 参考基因组
  - `path fasta_index` - 参考基因组索引
  - `val scan_name` - 扫描名称
- **输出**:
  - `loci` - MSI位点文件

### MSISENSORPRO
- **路径**: `modules/msisensorpro/main.nf`
- **用途**: MSI分析
- **输入**:
  - `tuple val(sample_id), path(alignment), path(index)` - BAM/CRAM
  - `path fasta` - 参考基因组
  - `path fasta_index` - 参考基因组索引
  - `path baseline` - MSI基线
- **输出**:
  - `msi_result` - MSI分析结果

### MSISENSORPRO_BASELINE
- **路径**: `modules/msisensorpro/main.nf`
- **用途**: 创建MSI基线
- **输入**:
  - `path normal_bams` - 多个正常样本BAM
  - `path normal_bais` - 索引文件
  - `path fasta` - 参考基因组
  - `path fasta_index` - 参考基因组索引
  - `path loci` - MSI位点文件
  - `val prefix` - 输出前缀
- **输出**:
  - `baseline` - MSI基线文件

## 12. HLA分型模块

### OPTITYPE
- **路径**: `modules/optitype/main.nf`
- **用途**: HLA分型
- **输入**:
  - `tuple val(sample_id), path(read1), path(read2)` - FASTQ文件
- **输出**:
  - `hla_result` - HLA分型结果

## 13. LOH分析模块

### LOHHLA
- **路径**: `modules/lohhla/main.nf`
- **用途**: 检测HLA Loss of Heterozygosity (LOH)
- **输入**:
  - `tuple val(sample_id), path(tumor_bam), path(tumor_bai)` - 肿瘤样本
  - `tuple val(sample_id), path(normal_bam), path(normal_bai)` - 正常样本
  - `path solutions` - 拷贝数结果文件
  - `path hla_results` - HLA分型结果
  - `val genomeAssembly` - 基因组版本('hg19'或'grch38')
- **输出**:
  - `lohhla_result` - LOHHLA原始结果
  - `lohhla_filtered` - 过滤后的LOH结果

## 14. VCF过滤模块

### BCFTOOLS
- **路径**: `modules/bcftools/main.nf`
- **用途**: VCF过滤
- **输入**:
  - `tuple val(sample_id), path(vcf)` - VCF文件
  - `val publish_dir` - 发布目录
  - `val minAD` - 最小等位基因深度
  - `val minAF` - 最小等位基因频率
- **输出**:
  - `filtered_vcf` - 过滤后的VCF

## 工作流

### SINGLE - 单样本分析
- **路径**: `workflows/single/main.nf`
- **用途**: 单样本全流程分析

### PAIRED - 配对样本分析
- **路径**: `workflows/paired/main.nf`
- **用途**: 肿瘤-正常配对样本分析

### CNV_BASELINE - CNV基线建立
- **路径**: `workflows/cnv_baseline/main.nf`
- **用途**: 建立CNV分析基线

### MSI_BASELINE - MSI基线建立
- **路径**: `workflows/msi_baseline/main.nf`
- **用途**: 建立MSI分析基线