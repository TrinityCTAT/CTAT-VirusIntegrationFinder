version 1.0

import "https://raw.githubusercontent.com/broadinstitute/CTAT-VirusIntegrationFinder/Terra-1.0.1c/WDL/ctat_VIF.wdl" as ctat_VIF_wf


struct CTAT_VIF_config {

  File ref_genome_fasta
  File ref_genome_gtf
  File viral_fasta
  File star_index_human_only
  File star_index_human_plus_virus
  File NULL_file

}


workflow ctat_VIF_Terra {

  input {
    String sample_id
    File? left
    File? right
    File? drs_path_fastqs
    Boolean clean_reads = true
    String docker = docker
    CTAT_VIF_config pipe_inputs_config
    Int preemptible
    
    }



   if (defined(drs_path_fastqs)) {
     call unpack_drs {
       input:
         sample_id = sample_id,
         drs_path_fastqs = select_first([drs_path_fastqs]),
         docker = docker,
         preemptible = preemptible
     }
  }
  
  call ctat_VIF_wf.ctat_vif as vif {
    input:     
      sample_id = sample_id,
      left = select_first([unpack_drs.left_fq, left]),
      right = select_first([unpack_drs.right_fq, right, pipe_inputs_config.NULL_file]),
      clean_reads = clean_reads,
      docker = docker,
      preemptible = preemptible,
    
      ref_genome_fasta = pipe_inputs_config.ref_genome_fasta,
      ref_genome_gtf = pipe_inputs_config.ref_genome_gtf,
      viral_fasta = pipe_inputs_config.viral_fasta,
      star_index_human_only = pipe_inputs_config.star_index_human_only,
      star_index_human_plus_virus = pipe_inputs_config.star_index_human_plus_virus
    

   }


   output {

     # STAR_init_hgOnly        
     File? star_init_hgOnly_bam = vif.star_init_hgOnly_bam
     File? star_init_hgOnly_bam_index = vif.star_init_hgOnly_bam_index
     File? star_init_hgOnly_log_final = vif.star_init_hgOnly_log_final
     File? star_init_hgOnly_SJ = vif.star_init_hgOnly_SJ
     File? star_init_hgOnly_ummapped_left_fq = vif.star_init_hgOnly_ummapped_left_fq
     File? star_init_hgOnly_unmapped_right_fq = vif.star_init_hgOnly_unmapped_right_fq

     
     # STAR_init_hgPlusVirus
     File? star_init_hgPlusVirus_bam = vif.star_init_hgPlusVirus_bam
     File? star_init_hgPlusVirus_bam_index = vif.star_init_hgPlusVirus_bam_index
     File? star_init_hgPlusVirus_log_final = vif.star_init_hgPlusVirus_log_final
     File? star_init_hgPlusVirus_SJ = vif.star_init_hgPlusVirus_SJ
     File? star_init_hgPlusVirus_chimeric_junction = vif.star_init_hgPlusVirus_chimeric_junction
     
     File? insertion_site_candidates_full = vif.insertion_site_candidates_full
     File? insertion_site_candidates_filtered = vif.insertion_site_candidates_filtered
     File? insertion_site_candidates_full_abridged = vif.insertion_site_candidates_full_abridged
     File? insertion_site_candidates_filtered_abridged = vif.insertion_site_candidates_filtered_abridged
     File? insertion_site_candidates_genome_chimeric_evidence_reads_bam = vif.insertion_site_candidates_genome_chimeric_evidence_reads_bam
     File? insertion_site_candidates_genome_chimeric_evidence_reads_bai = vif.insertion_site_candidates_genome_chimeric_evidence_reads_bai

     File? genome_abundance_plot = vif.genome_abundance_plot
     File? virus_coverage_read_counts_summary = vif.virus_coverage_read_counts_summary
     File? virus_coverage_read_counts_image = vif.virus_coverage_read_counts_image
     File? virus_coverage_read_counts_log_image = vif.virus_coverage_read_counts_log_image
     Array[File]? virus_coverage_virus_images = vif.virus_coverage_virus_images
     File? igv_virus_report_html = vif.igv_virus_report_html

     File? fasta_extract = vif.fasta_extract
     File? gtf_extract = vif.gtf_extract

     File? star_validate_inserts_bam = vif.star_validate_inserts_bam
     File? star_validate_inserts_bam_index = vif.star_validate_inserts_bam_index

     File? evidence_counts =  vif.evidence_counts
     File? evidence_bam = vif.evidence_bam
     File? evidence_bai = vif.evidence_bai

     File? refined_counts = vif.refined_counts
     File? genome_abundance_refined_plot = vif.genome_abundance_refined_plot
     File? igv_report_html = vif.igv_report_html


   }


}


task unpack_drs {
  input {
    String sample_id
    File drs_path_fastqs
    String docker
    Int preemptible
  }

  command <<<

    set -ex

    python <<CODE
    import sys, os, re
    import glob
    import subprocess

    sample_id = "~{sample_id}"
    fastqs_tar = "~{drs_path_fastqs}"

    os.makedirs("fastq")
    
    subprocess.check_call(f"tar -xvf {fastqs_tar} -C fastq", shell=True)

    fq_files = sorted(glob.glob("fastq/*"))

    def is_gzipped(filename):
        if re.search(".gz$", filename):
            return True
        else:
            return False
    
    if len(fq_files) == 1:
        if is_gzipped(fq_files[0]):
            os.rename(fq_files[0], sample_id + "_1.fq.gz")
        else:
            os.rename(fq_files[0], sample_id + "_1.fq")
            subprocess.check_call("gzip " + sample_id + "_1.fq", shell=True)

    

    elif len(fq_files) == 2:
        if is_gzipped(fq_files[0]):
            os.rename(fq_files[0], sample_id + "_1.fq.gz")
            os.rename(fq_files[1], sample_id + "_2.fq.gz")
        else:
            os.rename(fq_files[0], sample_id + "_1.fq")
            os.rename(fq_files[1], sample_id + "_2.fq")
            subprocess.check_call(f"gzip {sample_id}_*.fq", shell=True)

    elif len(fq_files) == 4:
        method = "cat"
        if is_gzipped(fq_files[0]):
            method = "zcat"

        subprocess.check_call(f"{method} {fq_files[0]} {fq_files[1]} | gzip -c > {sample_id}_1.fq.gz", shell=True)
        subprocess.check_call(f"{method} {fq_files[2]} {fq_files[3]} | gzip -c > {sample_id}_2.fq.gz", shell=True)

    CODE

    >>>

    output {
      File left_fq = "~{sample_id}_1.fq.gz"
      File? right_fq = "~{sample_id}_2.fq.gz"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(50 + size(drs_path_fastqs, "GB")*10 + 1) + " HDD"
        docker: docker
        cpu: 1
        memory: "4GB"
        maxRetries: 3
    }
}

