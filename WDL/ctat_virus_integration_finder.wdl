version 1.0

workflow ctat_virus_integration_finder {
    input {
        String sample_name
        File genome_lib_tar
        File viral_db_fasta
        File? rnaseq_aligned_bam
        File? fastq_pair_tar_gz
        File? left_fq
        File? right_fq
        String docker
        String memory =  "50G"
        String disk = "500 SSD"
        Int cpu = 10
        Int preemptible = 2
        Boolean remove_duplicates = false
    }

    if (defined(rnaseq_aligned_bam)) {

        call CTAT_BAM_TO_FASTQ {
            input:
                input_bam=select_first([rnaseq_aligned_bam]),
                sample_name=sample_name,
                docker=docker,
                preemptible=preemptible
        }
    }

    if (defined(fastq_pair_tar_gz)) {

        call CTAT_UNTAR_FASTQS {
            input:
                fastq_pair_tar_gz=select_first([fastq_pair_tar_gz]),
                docker=docker,
                preemptible=preemptible
        }
    }

    File left_fq_use = select_first([left_fq, CTAT_UNTAR_FASTQS.left_fq, CTAT_BAM_TO_FASTQ.left_fq])
    File right_fq_use = select_first([right_fq, CTAT_UNTAR_FASTQS.right_fq, CTAT_BAM_TO_FASTQ.right_fq])

    call CTAT_VIF {
        input:
            sample_name=sample_name,
            genome_lib_tar=genome_lib_tar,
            viral_db_fasta=viral_db_fasta,
            left_fq=left_fq_use,
            right_fq=right_fq_use,
            remove_duplicates=remove_duplicates,
            docker=docker,
            cpu=cpu,
            memory=memory,
            preemptible=preemptible,
            disk=disk
    }
    output {
        File virus_insertion_candidates_tsv=CTAT_VIF.virus_insertion_candidates_tsv
        File virus_read_counts_summary=CTAT_VIF.virus_read_counts_summary
        File vif_igv_report_html=CTAT_VIF.vif_igv_report_html
        File vif_aligned_reads_bam=CTAT_VIF.vif_aligned_reads_bam
        File vif_chimeric_targets_bed=CTAT_VIF.vif_chimeric_targets_bed
        File vif_chimeric_targets_fa=CTAT_VIF.vif_chimeric_targets_fa
        File virus_aligned_reads_bam=CTAT_VIF.virus_aligned_reads_bam
        File virus_igv_report_html=CTAT_VIF.virus_igv_report_html
        File report_html=CTAT_VIF.report_html
    }

}

task CTAT_BAM_TO_FASTQ {
    input {
        File input_bam
        String sample_name
        String docker
        Int preemptible
    }

    command <<<

        set -e

        # initial potential cleanup of read names in the bam file
        /usr/local/bin/sam_readname_cleaner.py ~{input_bam} ~{input_bam}.cleaned.bam


        # revert aligned bam
        java -Xmx1000m -jar /usr/local/src/picard.jar \
        RevertSam \
        INPUT=~{input_bam}.cleaned.bam \
        OUTPUT_BY_READGROUP=false \
        VALIDATION_STRINGENCY=SILENT \
        SORT_ORDER=queryname \
        OUTPUT=~{sample_name}.reverted.bam


        # bam to fastq
        java -jar /usr/local/src/picard.jar \
        SamToFastq I=~{sample_name}.reverted.bam \
        F=~{sample_name}_1.fastq F2=~{sample_name}_2.fastq \
        INTERLEAVE=false NON_PF=true \
        CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2

    >>>

    output {
        File left_fq="~{sample_name}_1.fastq"
        File right_fq="~{sample_name}_2.fastq"
    }

    runtime {
        docker: docker
        disks: "local-disk 500 SSD"
        memory: "20G"
        cpu: "16"
        preemptible: preemptible

    }

}

task CTAT_UNTAR_FASTQS {
    input {
        File fastq_pair_tar_gz
        String docker
        Int preemptible
    }

    command <<<

        set -e

        # untar the fq pair
        tar xvf ~{fastq_pair_tar_gz}
    >>>

    output {
        File left_fq = select_first(glob("*_1.fastq*"))
        File right_fq = select_first(glob("*_2.fastq*"))
    }

    runtime {
        docker: docker
        disks: "local-disk 500 SSD"
        memory: "2G"
        cpu: 1
        preemptible: preemptible

    }
}

task CTAT_VIF {
    input {
        File genome_lib_tar
        String sample_name
        File left_fq
        File right_fq
        File viral_db_fasta
        String docker
        Int cpu
        String memory
        Int preemptible
        String disk
        Boolean remove_duplicates

    }
    String output_prefix = if(remove_duplicates) then "VIF/" + sample_name + ".DupsRm"  else "VIF/" + sample_name

    command <<<

        set -e

        # untar the genome lib
        tar xvf ~{genome_lib_tar}

        ctat-VIF.py \
        --left_fq ~{left_fq} \
        --right_fq ~{right_fq} \
        --viral_db_fasta ~{viral_db_fasta} \
        --CPU $(nproc) \
        --genome_lib_dir `pwd`/ctat_genome_lib_build_dir \
        -O VIF --out_prefix ~{sample_name} \
        ~{true='--remove_duplicates' false='' remove_duplicates}
    >>>

    output {
        File virus_insertion_candidates_tsv="~{output_prefix}.insertion_site_candidates.tsv"
        File virus_read_counts_summary="~{output_prefix}.virus_read_counts_summary.tsv"
        File vif_aligned_reads_bam="~{output_prefix}.reads.bam"
        File vif_chimeric_targets_bed="~{output_prefix}.bed"
        File vif_chimeric_targets_fa="~{output_prefix}.fa"
        File virus_aligned_reads_bam="~{output_prefix}.virus.reads.bam"
        File virus_igv_report_html="~{output_prefix}.virus.igvjs.html"
        File vif_igv_report_html="~{output_prefix}.igvjs.html"
        File report_html="~{output_prefix}.html"
    }

    runtime {
        docker: docker
        disks: "local-disk ~{disk}"
        memory: memory
        cpu: cpu
        preemptible: preemptible
    }

}

