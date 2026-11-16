version 1.2

task FixBed {
    input {
        File bed
    }

    command <<<
        python3 /opt/schema-somatic/scripts/process_bed.py -i ~{bed} -o fixed.bed
    >>>

    output {
        File fixed_bed = "fixed.bed"
    }
}
