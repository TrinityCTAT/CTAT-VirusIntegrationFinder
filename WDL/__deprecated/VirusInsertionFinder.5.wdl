task CTAT_BAM_TO_FASTQ {

    File input_bam
    String sample_name
    String Docker
  
    
    command {

    set -e

    # initial potential cleanup of read names in the bam file
    /usr/local/bin/sam_readname_cleaner.py ${input_bam} ${input_bam}.cleaned.bam


    # revert aligned bam
    java -Xmx1000m -jar /usr/local/src/picard.jar \
        RevertSam \
        INPUT=${input_bam}.cleaned.bam \
        OUTPUT_BY_READGROUP=false \
        VALIDATION_STRINGENCY=SILENT \
        SORT_ORDER=queryname \
        OUTPUT=${sample_name}.reverted.bam 


    # bam to fastq
    java -jar /usr/local/src/picard.jar \
        SamToFastq I=${sample_name}.reverted.bam \
        F=${sample_name}_1.fastq F2=${sample_name}_2.fastq \
        INTERLEAVE=false NON_PF=true \
        CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2

   }
    
    output {
      File left_fq="${sample_name}_1.fastq"
      File right_fq="${sample_name}_2.fastq"
    }
    
    runtime {
            docker: Docker
            disks: "local-disk 500 SSD"
            memory: "20G"
            cpu: "16"
            preemptible: 0
            maxRetries: 3
    }
    
}



task CTAT_UNTAR_FASTQS {

  File fastq_pair_tar_gz
  String Docker

  command {

     set -e
    
     # untar the fq pair
     tar xvf ${fastq_pair_tar_gz}
  }

  output {
    File left_fq = select_first(glob("*_1.fastq*"))
    File right_fq = select_first(glob("*_2.fastq*"))
  }

  runtime {
            docker: Docker
            disks: "local-disk 500 SSD"
            memory: "10G"
            cpu: "4"
            preemptible: 0
            maxRetries: 3
    }
}



task CTAT_VIF {

    File genome_lib_tar
    String sample_name
    File left_fq
    File right_fq
    File viral_db_fasta
    String Docker
    
    command {

    set -e

    # untar the genome lib
    tar xvf ${genome_lib_tar}
      
    ctat-VIF.py \
         --left_fq ${left_fq} \
         --right_fq ${right_fq} \
         --viral_db_fasta ${viral_db_fasta} \
         --CPU 10 \
         --genome_lib_dir `pwd`/ctat_genome_lib_build_dir \
         -O VIF --out_prefix ${sample_name}
      

    }
    
    output {
      File  prelim_insertion_candidates_tsv="VIF/prelim.vif.rmdups-False.abridged.tsv"
      Array[File] virus_insertion_candidates_tsv=glob("VIF/${sample_name}.insertion_site_candidates.tsv")
      Array[File] prelim_virus_insertion_plot=glob("VIF/prelim.vif.rmdups-False.abridged.tsv.genome_plot.png")
      Array[File] virus_insertion_plot=glob("VIF/${sample_name}.insertion_site_candidates.genome_plot.png")
      Array[File] virus_read_counts_summary=glob("VIF/${sample_name}.virus_read_counts_summary.tsv")
      Array[File] virus_genome_coverage_plots=glob("VIF/${sample_name}.virus_coverage_plots.pdf")
      Array[File] vif_igv_report_html=glob("VIF/${sample_name}.igvjs.html")
      Array[File] vif_aligned_reads_bam=glob("VIF/${sample_name}.reads.bam")
      Array[File] vif_chimeric_targets_bed=glob("VIF/${sample_name}.bed")
      Array[File] vif_chimeric_targets_fa=glob("VIF/${sample_name}.fa")
    }
    

    runtime {
            docker: Docker
            disks: "local-disk 500 SSD"
            memory: "50G"
            cpu: "10"
            preemptible: 0
    }
    
    
}



workflow ctat_VIF_wf {

    String sample_name
    File genome_lib_tar
    File viral_db_fasta
    File? rnaseq_aligned_bam
    File? fastq_pair_tar_gz
    File? left_fq
    File? right_fq
	String? Docker

    
  
    if (defined(rnaseq_aligned_bam)) {
      
      call CTAT_BAM_TO_FASTQ {
            input:
              input_bam=rnaseq_aligned_bam,
              sample_name=sample_name,
              Docker=Docker
        }
    }

    if (defined(fastq_pair_tar_gz)) {

      call CTAT_UNTAR_FASTQS {
        input:
          fastq_pair_tar_gz=fastq_pair_tar_gz,
          Docker=Docker
      }
    }

    
    File? left_fq_use = select_first([left_fq, CTAT_UNTAR_FASTQS.left_fq, CTAT_BAM_TO_FASTQ.left_fq])
    File? right_fq_use = select_first([right_fq, CTAT_UNTAR_FASTQS.right_fq, CTAT_BAM_TO_FASTQ.right_fq])

    
    call CTAT_VIF {
      input:
        sample_name=sample_name,
        genome_lib_tar=genome_lib_tar,
        viral_db_fasta=viral_db_fasta,
        left_fq=left_fq_use,
        right_fq=right_fq_use,
        Docker=Docker
    }
    
    output {
      File  wf_prelim_insertion_candidates_tsv=CTAT_VIF.prelim_insertion_candidates_tsv
      Array[File] wf_virus_insertion_candidates_tsv=CTAT_VIF.virus_insertion_candidates_tsv
      Array[File] wf_prelim_virus_insertion_plot=CTAT_VIF.prelim_virus_insertion_plot
      Array[File] wf_virus_insertion_plot=CTAT_VIF.virus_insertion_plot
      Array[File] wf_virus_read_counts_summary=CTAT_VIF.virus_read_counts_summary
      Array[File] wf_virus_genome_coverage_plots=CTAT_VIF.virus_genome_coverage_plots
      Array[File] wf_vif_igv_report_html=CTAT_VIF.vif_igv_report_html
      Array[File] wf_vif_aligned_reads_bam=CTAT_VIF.vif_aligned_reads_bam
      Array[File] wf_vif_chimeric_targets_bed=CTAT_VIF.vif_chimeric_targets_bed
      Array[File] wf_vif_chimeric_targets_fa=CTAT_VIF.vif_chimeric_targets_fa
    }

}

