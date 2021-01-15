version 1.0

workflow ctat_virus_LRsplice {
    input {
        String sample_name
        File genome_lib_tar
        File? rnaseq_aligned_bam
        File? fastq_pair_tar_gz
        File? left_fq
        File? right_fq
        String output_prefix
        File virus_insertions_tsv
        Int flank
      
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

    call virus_LRsplice {
        input:
            sample_name=sample_name,
            output_prefix=output_prefix,
            virus_insertions_tsv=virus_insertions_tsv,
            flank=flank,
            genome_lib_tar=genome_lib_tar,
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
      

      File LRsplice_gtf = virus_LRsplice.LRsplice_gtf
      File LRsplice_fasta = virus_LRsplice.LRsplice_fasta
      File LRsplice_virus_only_gtf = virus_LRsplice.LRsplice_virus_only_gtf
      File LRsplice_host_only_bam = virus_LRsplice.LRsplice_host_only_bam
      File LRsplice_host_only_bai = virus_LRsplice.LRsplice_host_only_bai
      File LRsplice_virus_only_bam = virus_LRsplice.LRsplice_virus_only_bam
      File LRsplice_virus_only_bai = virus_LRsplice.LRsplice_virus_only_bai
      File LRsplice_host_virus_fusions_bam = virus_LRsplice.LRsplice_host_virus_fusions_bam
      File LRsplice_host_virus_fusions_bai = virus_LRsplice.LRsplice_host_virus_fusions_bai

            
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

task virus_LRsplice {
    input {
        File genome_lib_tar
        String sample_name
        String output_prefix
        File virus_insertions_tsv
        File left_fq
        File right_fq
        String docker
        Int cpu
        String memory
        Int preemptible
        String disk
        Boolean remove_duplicates
        Int flank

    }
    
    command <<<

        set -e

        # untar the genome lib
        tar xvf ~{genome_lib_tar}

        /usr/local/bin/util/ctat-VIF.longrange_virus_splice.py \
            --left_fq ~{left_fq} \
            --right_fq ~{right_fq} \
            --virus_insertions_tsv ~{virus_insertions_tsv} \
            --viral_db_fasta `pwd`/ctat_genome_lib/viruses/HPVs_db.fasta \
            --viral_db_gtf `pwd`/ctat_genome_lib/viruses/HPVs_db.annots.gtf \
            --flank ~{flank} \
            --CPU ~{cpu} \
            --genome_lib_dir `pwd`/ctat_genome_lib \
            -O ~{sample_name} --out_prefix ~{output_prefix}.LRsplice \
            ~{true='--remove_duplicates' false='' remove_duplicates}
            
    >>>

    output {

      File LRsplice_gtf = "~{sample_name}/~{output_prefix}.LRsplice.gtf"
      File LRsplice_fasta = "~{sample_name}/~{output_prefix}.LRsplice.fasta"
      File LRsplice_virus_only_gtf = "~{sample_name}/~{output_prefix}.LRsplice.virus-only.gtf"
      File LRsplice_host_only_bam = "~{sample_name}/~{output_prefix}.LRsplice.host-only.bam"
      File LRsplice_host_only_bai = "~{sample_name}/~{output_prefix}.LRsplice.host-only.bam.bai"
      File LRsplice_virus_only_bam = "~{sample_name}/~{output_prefix}.LRsplice.virus-only.bam"
      File LRsplice_virus_only_bai = "~{sample_name}/~{output_prefix}.LRsplice.virus-only.bam.bai"
      File LRsplice_host_virus_fusions_bam = "~{sample_name}/~{output_prefix}.LRsplice.host-virus.fusions.bam"
      File LRsplice_host_virus_fusions_bai = "~{sample_name}/~{output_prefix}.LRsplice.host-virus.fusions.bam.bai"

            
    }

    runtime {
        docker: docker
        disks: "local-disk ~{disk}"
        memory: memory
        cpu: cpu
        preemptible: preemptible
    }

}

