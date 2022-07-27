version 1.0

import "https://raw.githubusercontent.com/broadinstitute/CTAT-VirusIntegrationFinder/Terra-1.5.0-predev/WDL/Terra/ctat_VIF.Terra.wdl" as ctat_VIF_Terra
#import "ctat_VIF.Terra.wdl" as ctat_VIF_Terra


workflow ctat_VIF_Terra_hg38 {

  input {
    String sample_id
    File? left
    File? right
    File? drs_path_fastqs
    Boolean clean_reads = true
    Int max_hits = 50
    String docker = "trinityctat/ctat_vif:latest"
    Int preemptible = 1


    Int star_cpu = 12
    Float star_init_memory = 50
    Float star_validate_memory= 75

    
    ####################
    ## Kickstart options:
      
    # run stage 2 only
    File? hg_unmapped_left_fq
    File? hg_unmapped_right_fq

    # run state 3 only (needs stage 2 inputs above too!)
    File? human_virus_chimJ
    File? human_virus_bam
    File? human_virus_bai

    
    CTAT_VIF_config pipe_inputs_config = {
      "ref_genome_fasta" : "gs://ctat_genome_libs/GRCh38_gencode_v22/03-01-2021/ref_genome.fa",
      "ref_genome_gtf" : "gs://ctat_genome_libs/GRCh38_gencode_v22/03-01-2021/ref_annot.gtf",
      "star_index_human_only" : "gs://ctat_genome_libs/GRCh38_gencode_v22/03-01-2021/ref_genome.fa.star.idx.tar",
      "star_index_human_plus_virus" : "gs://ctat_genome_libs/GRCh38_gencode_v22/06-08-2022/hg_plus_viraldb.fasta.star.idx.tar",
	  "viral_fasta" : "gs://ctat_genome_libs/GRCh38_gencode_v22/06-08-2022/virus_db.fasta",
	  "NULL_file" : "gs://ctat_genome_libs/null" 
    }

  }
  
  call ctat_VIF_Terra.ctat_VIF_Terra as vif_terra {
    input:     
      sample_id = sample_id,
      left = left,
      right = right,
      clean_reads = clean_reads,
      max_hits = max_hits,
      drs_path_fastqs = drs_path_fastqs,
      docker = docker,
      pipe_inputs_config = pipe_inputs_config,
      preemptible = preemptible,

      star_cpu = star_cpu,
      star_init_memory = star_init_memory,
      star_validate_memory = star_validate_memory,
      
      hg_unmapped_left_fq = hg_unmapped_left_fq,
      hg_unmapped_right_fq = hg_unmapped_right_fq,

      human_virus_chimJ = human_virus_chimJ,
      human_virus_bam = human_virus_bam,
      human_virus_bai = human_virus_bai
    
   }


output {

     # STAR_init_hgOnly        
     File? star_init_hgOnly_bam = vif_terra.star_init_hgOnly_bam
     File? star_init_hgOnly_bam_index = vif_terra.star_init_hgOnly_bam_index
     File? star_init_hgOnly_log_final = vif_terra.star_init_hgOnly_log_final
     File? star_init_hgOnly_SJ = vif_terra.star_init_hgOnly_SJ
     File? star_init_hgOnly_ummapped_left_fq = vif_terra.star_init_hgOnly_ummapped_left_fq
     File? star_init_hgOnly_unmapped_right_fq = vif_terra.star_init_hgOnly_unmapped_right_fq

     
     # STAR_init_hgPlusVirus
     File? star_init_hgPlusVirus_bam = vif_terra.star_init_hgPlusVirus_bam
     File? star_init_hgPlusVirus_bam_index = vif_terra.star_init_hgPlusVirus_bam_index
     File? star_init_hgPlusVirus_log_final = vif_terra.star_init_hgPlusVirus_log_final
     File? star_init_hgPlusVirus_SJ = vif_terra.star_init_hgPlusVirus_SJ
     File? star_init_hgPlusVirus_chimeric_junction = vif_terra.star_init_hgPlusVirus_chimeric_junction
     
     File? insertion_site_candidates_full = vif_terra.insertion_site_candidates_full
     File? insertion_site_candidates_filtered = vif_terra.insertion_site_candidates_filtered
     File? insertion_site_candidates_full_abridged = vif_terra.insertion_site_candidates_full_abridged
     File? insertion_site_candidates_filtered_abridged = vif_terra.insertion_site_candidates_filtered_abridged
     File? insertion_site_candidates_genome_chimeric_evidence_reads_bam = vif_terra.insertion_site_candidates_genome_chimeric_evidence_reads_bam
     File? insertion_site_candidates_genome_chimeric_evidence_reads_bai = vif_terra.insertion_site_candidates_genome_chimeric_evidence_reads_bai
     File? insertion_site_candidates_human_virus_chimJ = vif_terra.insertion_site_candidates_human_virus_chimJ
     
     File? genome_abundance_plot = vif_terra.genome_abundance_plot
     File? virus_coverage_read_counts_summary = vif_terra.virus_coverage_read_counts_summary
     File? virus_coverage_read_counts_image = vif_terra.virus_coverage_read_counts_image
     File? virus_coverage_read_counts_log_image = vif_terra.virus_coverage_read_counts_log_image
     Array[File]? virus_coverage_virus_images = vif_terra.virus_coverage_virus_images
     File? igv_virus_report_html = vif_terra.igv_virus_report_html
     File? virus_alignments_bam = vif_terra.virus_alignments_bam
     File? virus_alignments_bai = vif_terra.virus_alignments_bai

     File? fasta_extract = vif_terra.fasta_extract
     File? gtf_extract = vif_terra.gtf_extract

     File? star_validate_inserts_bam = vif_terra.star_validate_inserts_bam
     File? star_validate_inserts_bam_index = vif_terra.star_validate_inserts_bam_index

     File? evidence_counts =  vif_terra.evidence_counts
     File? evidence_bam = vif_terra.evidence_bam
     File? evidence_bai = vif_terra.evidence_bai

     File? prelim_refined_counts =  vif_terra.prelim_refined_counts
     File? refined_counts = vif_terra.refined_counts
     File? refined_distilled = vif_terra.refined_distilled
     File? genome_abundance_refined_plot = vif_terra.genome_abundance_refined_plot
     File? igv_report_html = vif_terra.igv_report_html


   }

   
}

