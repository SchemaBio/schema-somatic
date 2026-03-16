// LOHHLA HLA Loss of Heterozygosity 分析模块
// 用途：检测肿瘤样本的HLA LOH事件
// 用法：
//   - 输入：肿瘤和正常样本BAM/CRAM文件及索引、HLA分型结果、拷贝数结果
//   - 输出：HLA LOH预测结果（原始结果和过滤结果）

process LOHHLA {
    tag "LOHHLA on $sample_id"
    label 'lohhla'
    publishDir "${params.output}/03.Mutations/LOHHLA", mode: 'copy'

    input:
        tuple val(sample_id), path(tumor_bam), path(tumor_bai)
        tuple val(sample_id), path(normal_bam), path(normal_bai)
        path solutions  // 拷贝数结果文件
        path hla_results  // HLA分型结果
        val genomeAssembly  // 基因组版本：hg19 或 grch38

    output:
        path("${sample_id}.lohhla.txt"), emit: lohhla_result
        path("${sample_id}.loh.txt"), emit: lohhla_filtered

    script:
    def genome = genomeAssembly ?: 'hg19'
    """
    date=\$(date +%Y%m%d)
    runningDir=\$PWD
    mkdir -p ${sample_id}_tmp/bam

    # 创建符号链接
    ln -s ${tumor_bam} ${runningDir}/${sample_id}_tmp/bam/${sample_id}.bam
    ln -s ${tumor_bai} ${runningDir}/${sample_id}_tmp/bam/${sample_id}.bam.bai
    ln -s ${normal_bam} ${runningDir}/${sample_id}_tmp/bam/${sample_id}_normal.bam
    ln -s ${normal_bai} ${runningDir}/${sample_id}_tmp/bam/${sample_id}_normal.bam.bai

    # 运行LOHHLA分析
    LOHHLAscript.R \\
        --patientId ${sample_id} \\
        --outputDir ${runningDir}/${sample_id}_tmp \\
        --normalBAMfile ${runningDir}/${sample_id}_tmp/bam/${sample_id}_normal.bam \\
        --tumorBAMfile ${runningDir}/${sample_id}_tmp/bam/${sample_id}.bam \\
        --BAMDir ${runningDir}/${sample_id}_tmp/bam \\
        --hlaPath ${hla_results} \\
        --CopyNumLoc ${solutions} \\
        --genomeAssembly ${genome} \\
        --minCoverageFilter 10 \\
        --mappingStep TRUE \\
        --fishingStep TRUE \\
        --cleanUp FALSE

    # 复制原始结果
    cp ${runningDir}/${sample_id}_tmp/${sample_id}.10.DNA.HLAlossPrediction_CI.\${date}.tsv ${sample_id}.lohhla.txt

    # 过滤结果 (使用R)
    Rscript -e '
    library(tools)
    input_file <- "${sample_id}.lohhla.txt"
    output_file <- "${sample_id}.loh.txt"

    loh_results <- file(output_file, "w")
    writeLines(c("HLAType", "HLACopyNumWithBAFBin", "pval", "LOHstat"), loh_results)

    lines <- readLines(input_file)
    for (line in lines) {
        if (!startsWith(line, "message")) {
            fields <- strsplit(line, "\t")[[1]]
            if (length(fields) >= 34) {
                HLAtype <- substr(fields[2], 1, 5)

                if (startsWith(fields[1], "homozygous_alleles")) {
                    HLAtype2copyNumWithBAFBin <- "-"
                    pVal <- "-"
                    LOHstat <- "FALSE"
                } else {
                    HLAtype2copyNumWithBAFBin <- fields[29]
                    pVal <- fields[34]

                    if (pVal == "NA") {
                        LOHstat <- "-"
                    } else {
                        tryCatch({
                            pval_num <- as.numeric(pVal)
                            copy_num <- as.numeric(HLAtype2copyNumWithBAFBin)
                            if (pval_num < 0.01 && copy_num < log2(0.5)) {
                                LOHstat <- "TRUE"
                            } else {
                                LOHstat <- "FALSE"
                            }
                        }, error = function(e) {
                            LOHstat <<- "-"
                        })
                    }
                }
                writeLines(c(HLAtype, HLAtype2copyNumWithBAFBin, pVal, LOHstat), loh_results)
            }
        }
    }
    close(loh_results)
    cat("LOH分析完成\n")
    '
    """
}