version 1.0

import "https://raw.githubusercontent.com/broadinstitute/CTAT-VirusIntegrationFinder/Terra-1.0.2/WDL/ctat_VIF.wdl" as ctat_VIF_wf


struct CTAT_VIF_config {

  File ref_genome_fasta
  File ref_genome_gtf
  File viral_fasta
  File star_index_human_only
  File star_index_human_plus_virus

}


workflow ctat_VIF_Terra {

  input {
    String sample_id
    File left
    File? right
    String docker = docker
    CTAT_VIF_config pipe_inputs_config

    }
  
  call ctat_VIF_wf.ctat_vif as vif {
    input:     
      sample_id = sample_id,
      left = left,
      right = right,
      docker = docker,
    
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

