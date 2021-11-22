version 1.0


workflow ctat_vif {
    input {
        String sample_id
      
        File left
        File? right

        File ref_genome_fasta
        File ref_genome_gtf

        File viral_fasta
        File? viral_gtf

        Boolean star_init_only = false

        Boolean remove_duplicates = true
        Boolean generate_reports = true

        Int min_reads = 10

        # star indices needed (local: point to directory name, on cloud: give tar file)
        File? star_index_human_only
        String? star_index_human_only_dirpath
        
        File? star_index_human_plus_virus
        String? star_index_human_plus_virus_dirpath 

      
        String igv_virus_reports_memory = "14GB"
        String igv_reports_memory = "1GB"

        # where to find utilities in the docker image.
        String util_dir = "/usr/local/src/CTAT-VirusIntegrationFinder/util"
        String picard = "/usr/local/src/picard.jar"

        Float star_extra_disk_space = 30
        Float star_fastq_disk_space_multiplier = 10
        Int sjdb_overhang = 150

        Boolean autodetect_cpu = true # auto-detect number of cpus for STAR as # of requested CPUs might not equal actual CPUs, depending on memory
        Boolean star_use_ssd = false

        Int star_cpu = 12

        # star init settings
        Float star_init_memory = 50
        String star_init_two_pass_mode = "Basic"

        # star validate settings
        Float star_validate_memory = 75
        String star_validate_two_pass_mode = "Basic" # or None


        # general runtime settings
        Int preemptible = 2
        String docker = "trinityctat/ctat_vif:latest"

        # run stage 2 only
        File? insertion_site_candidates
    }

    parameter_meta {
        left:{help:"One of the two paired RNAseq samples"}
        right:{help:"One of the two paired RNAseq samples"}
        min_reads:{help:"Filter insertion sites candidates that do not have at least 'min_reads'"}
        remove_duplicates:{help:"Remove duplicate alignments"}
        ref_genome_fasta:{help:"Host fasta"}
        ref_genome_gtf:{help:"Host annotations GTF"}
        viral_fasta:{help:"Viral fasta"}
        ref_genome_incl_viral:{help:"The star index already includes the viral genomes"}
        insertion_site_candidates:{help:"Previously generated candidates"}
        star_reference:{help:"STAR index archive containing both host and viral genomes"}
        star_cpu:{help:"STAR aligner number of CPUs"}
        star_memory:{help:"STAR aligner memory"}
        util_dir:{help:"Path to util directory (for non-Docker use)"}
        star_init_only:{help:"Only perform initial STAR chimeric junction analysis"}
        docker:{help:"Docker image"}
    }

    output {

        # STAR_init_hgOnly        
        File? star_init_hgOnly_bam = select_first([STAR_init_hgOnly.bam])
        File? star_init_hgOnly_bam_index = select_first([STAR_init_hgOnly.bai])
        File? star_init_hgOnly_log_final = STAR_init_hgOnly.output_log_final
        File? star_init_hgOnly_SJ = STAR_init_hgOnly.output_SJ
        File? star_init_hgOnly_chimeric_junction = STAR_init_hgOnly.chimeric_junction
        File? star_init_hgOnly_ummapped_left_fq = STAR_init_hgOnly.Unmapped_left_fq
        File? star_init_hgOnly_unmapped_right_fq = STAR_init_hgOnly.Unmapped_right_fq
        
        # STAR_init_hgPlusVirus
        File? star_init_hgPlusVirus_bam = select_first([STAR_init_hgPlusVirus.bam])
        File? star_init_hgPlusVirus_bam_index = select_first([STAR_init_hgPlusVirus.bai])
        File? star_init_hgPlusVirus_log_final = STAR_init_hgPlusVirus.output_log_final
        File? star_init_hgPlusVirus_SJ = STAR_init_hgPlusVirus.output_SJ
        File? star_init_hgPlusVirus_chimeric_junction = STAR_init_hgPlusVirus.chimeric_junction
      
        File? insertion_site_candidates_full = InsertionSiteCandidates.full
        File? insertion_site_candidates_full_filtered = InsertionSiteCandidates.full_filtered
        File? insertion_site_candidates_abridged = InsertionSiteCandidates.abridged
        File? insertion_site_candidates_abridged_filtered = InsertionSiteCandidates.abridged_filtered
        File? insertion_site_candidates_abridged_detailed = InsertionSiteCandidates.abridged_detailed
        File? insertion_site_candidates_full_read_stats = InsertionSiteCandidates.full_read_stats
        File? insertion_site_candidates_genome_chimeric_evidence_reads_bam = InsertionSiteCandidates.genome_chimeric_evidence_reads_bam
        File? insertion_site_candidates_genome_chimeric_evidence_reads_bai = InsertionSiteCandidates.genome_chimeric_evidence_reads_bai
        
        File? genome_abundance_plot = VirusReport.genome_abundance_plot
        File? virus_coverage_read_counts_summary = VirusReport.read_counts_summary
        File? virus_coverage_read_counts_image = VirusReport.read_counts_image
        File? virus_coverage_read_counts_log_image = VirusReport.read_counts_log_image
        Array[File]? virus_coverage_virus_images = VirusReport.virus_images
        File? igv_virus_report_html = VirusReport.html

        File? fasta_extract = ExtractChimericGenomicTargets.fasta_extract
        File? gtf_extract = ExtractChimericGenomicTargets.gtf_extract

        File? star_validate_inserts_bam = select_first([RemoveDuplicates2.bam, STAR_validate.bam, "/dev/null"])
        File? star_validate_inserts_bam_index = select_first([RemoveDuplicates2.bai, STAR_validate.bai, "/dev/null"])

        File? evidence_counts =  ChimericContigEvidenceAnalyzer.evidence_counts
        File? evidence_bam = ChimericContigEvidenceAnalyzer.evidence_bam
        File? evidence_bai = ChimericContigEvidenceAnalyzer.evidence_bai

        File? refined_counts = SummaryReport.refined_counts
        File? genome_abundance_refined_plot = SummaryReport.genome_abundance_plot
        File? igv_report_html = SummaryReport.html

    }

    
    if ( !defined(insertion_site_candidates) ) {
        call STAR_init as STAR_init_hgOnly {
            input:
                util_dir=util_dir,
                fastq1=left,
                fastq2=right,
                two_pass_mode = star_init_two_pass_mode,
                base_name=sample_id + ".hgOnly",
                star_reference=star_index_human_only,
                star_reference_dirpath = star_index_human_only_dirpath,
                extra_disk_space = star_extra_disk_space,
                disk_space_multiplier = star_fastq_disk_space_multiplier,
                memory = star_init_memory,
                use_ssd = star_use_ssd,
                cpu = star_cpu,
                autodetect_cpu = autodetect_cpu,
                docker = docker,
                preemptible = preemptible
            }

            call STAR_init as STAR_init_hgPlusVirus {
            input:
                util_dir=util_dir,
                fastq1=STAR_init_hgOnly.Unmapped_left_fq,
                fastq2=STAR_init_hgOnly.Unmapped_right_fq,
                two_pass_mode = star_init_two_pass_mode,
                base_name=sample_id + ".hgPlusVirus",
                star_reference=star_index_human_plus_virus,
                star_reference_dirpath = star_index_human_plus_virus_dirpath,
                extra_disk_space = star_extra_disk_space,
                disk_space_multiplier = star_fastq_disk_space_multiplier,
                memory = star_init_memory,
                use_ssd = star_use_ssd,
                cpu = star_cpu,
                autodetect_cpu = autodetect_cpu,
                docker = docker,
                preemptible = preemptible
        }

        call InsertionSiteCandidates {
            input:
                chimeric_junction=select_first([STAR_init_hgPlusVirus.chimeric_junction]),
                bam=select_first([STAR_init_hgPlusVirus.bam]),
                bai=select_first([STAR_init_hgPlusVirus.bai]),
                viral_fasta=viral_fasta,
                remove_duplicates=remove_duplicates,
                util_dir=util_dir,
                min_reads=min_reads,
                preemptible=preemptible,
                docker=docker,
                sample_id=sample_id
        }
        File insertion_site_candidates_output =  (if min_reads>0 then select_first([InsertionSiteCandidates.abridged_filtered]) else InsertionSiteCandidates.abridged)

        
        if(generate_reports) {
            call VirusReport {
                input:
                    bam=select_first([STAR_init_hgPlusVirus.bam]),
                    bai=select_first([STAR_init_hgPlusVirus.bai]),
                    remove_duplicates=remove_duplicates,
                    viral_fasta=viral_fasta,
                    insertion_site_candidates=insertion_site_candidates_output,
                    util_dir=util_dir,
                    memory=igv_virus_reports_memory,
                    preemptible=preemptible,
                    docker=docker,
                    sample_id=sample_id
            }
        }

    }

    if(!star_init_only) {
        File insertion_site_candidates_use = select_first([insertion_site_candidates, InsertionSiteCandidates.abridged_filtered, InsertionSiteCandidates.abridged])
        call ExtractChimericGenomicTargets {
            input:
                fasta=ref_genome_fasta,
                viral_fasta=viral_fasta,
                insertion_site_candidates_abridged=insertion_site_candidates_use,
                util_dir=util_dir,
                preemptible=preemptible,
                docker=docker,
                sample_id=sample_id
        }

        if (ExtractChimericGenomicTargets.has_chimeric_targets) {

        call STAR_validate {
            input:
                util_dir=util_dir,
                fastq1=select_first([STAR_init_hgOnly.Unmapped_left_fq, left]),
                fastq2=select_first([STAR_init_hgOnly.Unmapped_right_fq, right]),
                two_pass_mode = star_validate_two_pass_mode,
                base_name=sample_id+".validate_inserts",
                star_reference=star_index_human_plus_virus,                                                                                                                                            
                star_reference_dirpath = star_index_human_plus_virus_dirpath,
                insertions_fasta_file=ExtractChimericGenomicTargets.fasta_extract,
                extra_disk_space = star_extra_disk_space,
                disk_space_multiplier = star_fastq_disk_space_multiplier,
                memory = star_validate_memory,
                use_ssd = star_use_ssd,
                cpu = star_cpu,
                autodetect_cpu = autodetect_cpu,
                docker = docker,
                preemptible = preemptible
        }
        if(remove_duplicates) {
            call RemoveDuplicates as RemoveDuplicates2 {
                input:
                    input_bam=STAR_validate.bam,
                    input_bai=STAR_validate.bai,
                    output_bam = sample_id+".validate_inserts.sortedByCoord.out.rm.dups.bam",
                    util_dir=util_dir,
                    cpu=1,
                    preemptible=preemptible,
                    memory="2G",
                    docker=docker,
                    extra_disk_space=1,
                    disk_space_multiplier=2
            }
        }

        call ChimericContigEvidenceAnalyzer {
            input:
                bam=select_first([RemoveDuplicates2.bam, STAR_validate.bam]),
                bai=select_first([RemoveDuplicates2.bai, STAR_validate.bai]),
                gtf=ExtractChimericGenomicTargets.gtf_extract,
                min_reads=min_reads,
                util_dir=util_dir,
                preemptible=preemptible,
                docker=docker,
                sample_id=sample_id
        }

        if(generate_reports && ChimericContigEvidenceAnalyzer.insertions_recovered) {
            call SummaryReport {
                input:
                    init_counts=insertion_site_candidates_use,
                    vif_counts=ChimericContigEvidenceAnalyzer.evidence_counts,
                    alignment_bam=ChimericContigEvidenceAnalyzer.evidence_bam,
                    alignment_bai=ChimericContigEvidenceAnalyzer.evidence_bai,
                    chim_targets_gtf=ExtractChimericGenomicTargets.gtf_extract,
                    chim_targets_fasta=ExtractChimericGenomicTargets.fasta_extract,
                    gtf=ref_genome_gtf,
                    images=select_all([VirusReport.genome_abundance_plot, VirusReport.read_counts_image, VirusReport.read_counts_log_image]),
                    util_dir=util_dir,
                    memory=igv_reports_memory,
                    preemptible=preemptible,
                    docker=docker,
                    sample_id=sample_id
            }
          }
       }
    }
}



task STAR_init {
    input {
        String util_dir
        File fastq1
        File? fastq2
        File? star_reference
        String? star_reference_dirpath
        Float extra_disk_space
        Float disk_space_multiplier
        Boolean use_ssd
        Int cpu
        Int preemptible
        Float memory
        String docker
        String base_name
        String two_pass_mode
        Boolean autodetect_cpu
    }
    Int max_mate_dist = 100000
    
    command <<<
        set -e

        cpu=~{cpu}
        genomeDir="~{star_reference_dirpath}"
        if [[ "${genomeDir}" == "" ]]; then
            genomeDir="~{star_reference}"
        fi
      
        fastqs="~{fastq1} ~{fastq2}"
        readFilesCommand=""
        if [[ "~{fastq1}" == *.gz ]] ; then
            readFilesCommand="--readFilesCommand \"gunzip -c\""
        fi
      
        if [ "~{autodetect_cpu}" == "true" ]; then
            cpu=$(nproc)
        fi

        if [ -f "${genomeDir}" ] ; then
            mkdir genome_dir
            compress=""
            if [[ "${genomeDir}" == *.tar.gz ]] ; then
              compress="-I pigz"
            elif [[ "${genomeDir}" == *.tar.bz2 ]] ; then
                compress="-I pbzip2"
            fi
            tar $compress -xf ~{star_reference} -C genome_dir --strip-components 1
            genomeDir="genome_dir"
        fi

        # special case for tar of fastq files
        if [[ "~{fastq1}" == *.tar.gz ]] ; then
            mkdir fastq
            tar -I pigz -xvf ~{fastq1} -C fastq
            fastqs=$(find fastq -type f)
            readFilesCommand=""
            if [[ "$fastqs" = *.gz ]] ; then
                readFilesCommand="--readFilesCommand \"gunzip -c\""
            fi
        fi

      
      STAR \
            --runMode alignReads \
            --genomeDir $genomeDir \
            --runThreadN $cpu \
            --readFilesIn $fastqs \
            $readFilesCommand \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ~{base_name}. \
            --outSAMstrandField intronMotif \
            --outSAMunmapped Within \
            ~{"--twopassMode " + two_pass_mode} \
            --alignSJDBoverhangMin 10 \
            --genomeSuffixLengthMax 10000 \
            --limitBAMsortRAM 47271261705 \
            --alignInsertionFlush Right \
            --alignMatesGapMax ~{max_mate_dist} \
            --alignIntronMax  ~{max_mate_dist} \
            --peOverlapNbasesMin 12 \
            --peOverlapMMp 0.1 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --alignSplicedMateMapLminOverLmate 0 \
            --alignSplicedMateMapLmin 30 \
            --chimJunctionOverhangMin 12 \
             --chimOutJunctionFormat 0 \
             --chimSegmentMin 8 \
             --chimSegmentReadGapMax 3 \
             --chimScoreJunctionNonGTAG 0 \
             --chimNonchimScoreDropMin 10 \
             --chimMultimapScoreRange 10 \
             --chimMultimapNmax 20 \
             --chimOutType Junctions WithinBAM \
             --outReadsUnmapped Fastx 
      
      samtools index "~{base_name}.Aligned.sortedByCoord.out.bam"

      # always have at least the Unmapped.out.mate1 file
      touch ~{base_name}.Unmapped.out.mate1
    >>>

          
    
    output {
        File bam = "~{base_name}.Aligned.sortedByCoord.out.bam"
        File bai = "~{base_name}.Aligned.sortedByCoord.out.bam.bai"
        File output_log_final = "~{base_name}.Log.final.out"
        File output_SJ = "~{base_name}.SJ.out.tab"
        File? chimeric_junction = "~{base_name}.Chimeric.out.junction"
        File Unmapped_left_fq = "~{base_name}.Unmapped.out.mate1"
        File? Unmapped_right_fq = "~{base_name}.Unmapped.out.mate2"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(fastq1, "GB")*disk_space_multiplier + size(fastq2, "GB") * disk_space_multiplier + size(star_reference, "GB")*8 + extra_disk_space) + " " + (if use_ssd then "SSD" else "HDD")
        docker: docker
        cpu: cpu
        memory: memory + "GB"
    }
}



task STAR_validate {
    input {
        String util_dir
        File fastq1
        File? fastq2
        File? star_reference
        String? star_reference_dirpath
        Float extra_disk_space
        Float disk_space_multiplier
        Boolean use_ssd
        Int cpu
        Int preemptible
        Float memory
        String docker
        String base_name
        String two_pass_mode
        Boolean autodetect_cpu
        File insertions_fasta_file
    }
    Int max_mate_dist = 100000
    
    command <<<
        set -e

        cpu=~{cpu}
        genomeDir="~{star_reference_dirpath}"
        if [[ "${genomeDir}" == "" ]]; then
            genomeDir="~{star_reference}"
        fi

        fastqs="~{fastq1} ~{fastq2}"
        readFilesCommand=""
        if [[ "~{fastq1}" == *.gz ]] ; then
            readFilesCommand="--readFilesCommand \"gunzip -c\""
        fi
        if [ "~{autodetect_cpu}" == "true" ]; then
            cpu=$(nproc)
        fi


        if [ -f "${genomeDir}" ] ; then
            mkdir genome_dir
            compress=""
            if [[ $genomeDir == *.tar.gz ]] ; then
                compress="-I pigz"
            elif [[ $genomeDir == *.tar.bz2 ]] ; then
                compress="-I pbzip2"
            fi
            tar $compress -xf ~{star_reference} -C genome_dir --strip-components 1
            genomeDir="genome_dir"
        fi

        # special case for tar of fastq files
        if [[ "~{fastq1}" == *.tar.gz ]] ; then
            mkdir fastq
            tar -I pigz -xvf ~{fastq1} -C fastq
            fastqs=$(find fastq -type f)
            readFilesCommand=""
            if [[ "$fastqs" = *.gz ]] ; then
                readFilesCommand="--readFilesCommand \"gunzip -c\""
            fi
        fi

      STAR \
        --runMode alignReads \
        --genomeDir $genomeDir \
        --genomeFastaFiles ~{insertions_fasta_file} \
        --runThreadN $cpu \
        --readFilesIn $fastqs \
        $readFilesCommand \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ~{base_name}. \
        --outSAMstrandField intronMotif \
        --outSAMunmapped Within \
        ~{"--twopassMode " + two_pass_mode} \
        --alignSJDBoverhangMin 10 \
        --genomeSuffixLengthMax 10000 \
        --limitBAMsortRAM 47271261705 \
        --alignInsertionFlush Right \
        --outSAMfilter KeepOnlyAddedReferences \
        --alignMatesGapMax ~{max_mate_dist} \
        --alignIntronMax  ~{max_mate_dist} \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 \

        samtools index "~{base_name}.Aligned.sortedByCoord.out.bam"
    >>>

          
    
    output {
        File bam = "~{base_name}.Aligned.sortedByCoord.out.bam"
        File bai = "~{base_name}.Aligned.sortedByCoord.out.bam.bai"
        File output_log_final = "~{base_name}.Log.final.out"
        File output_SJ = "~{base_name}.SJ.out.tab"
        File? chimeric_junction = "~{base_name}.Chimeric.out.junction"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(fastq1, "GB")*disk_space_multiplier + size(fastq2, "GB") * disk_space_multiplier + size(star_reference, "GB")*8 + extra_disk_space) + " " + (if use_ssd then "SSD" else "HDD")
        docker: docker
        cpu: cpu
        memory: memory + "GB"
    }
}

task RemoveDuplicates {
    input {
        File input_bam
        File input_bai
        String util_dir
        Int cpu
        Int preemptible
        String memory
        String docker
        Float extra_disk_space
        Float disk_space_multiplier
        String output_bam
    }
    command <<<
        set -e

        ~{util_dir}/bam_mark_duplicates.py \
        -i ~{input_bam} \
        -o ~{output_bam} \
        -r

        samtools index ~{output_bam}
    >>>

    output {
        File bam = "~{output_bam}"
        File bai = "~{output_bam}.bai"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(input_bam, "GB")*disk_space_multiplier + extra_disk_space) + " HDD"
        docker: docker
        cpu: cpu
        memory: memory
    }
}

task InsertionSiteCandidates {
    input {
        File chimeric_junction
        File bam
        File bai
        File viral_fasta
        Boolean remove_duplicates
        String util_dir
        Int preemptible
        String docker
        Int min_reads
        String sample_id

    }
    # Create the prefix to add to files 
    String prefix = sample_id + ".vif.init"

    command <<<
        set -e

        ~{util_dir}/chimJ_to_virus_insertion_candidate_sites.py \
        --chimJ ~{chimeric_junction} \
        --patch_db_fasta ~{viral_fasta} \
        --output_prefix ~{prefix} \
        ~{true='--remove_duplicates' false='' remove_duplicates}

        python <<CODE
        min_reads = ~{min_reads}
        if min_reads > 0:
            import pandas as pd
            df = pd.read_csv("~{prefix}.abridged.tsv", sep='\t')
            df = df[df['total'] >= min_reads]
            df.to_csv("~{prefix}.abridged.filtered.tsv", sep='\t', index=False)

            df = pd.read_csv("~{prefix}.full.tsv", sep='\t')
            df = df[df['total'] >= min_reads]
            df.to_csv("~{prefix}.full.filtered.tsv", sep='\t', index=False)
      
        CODE

        
        # extract the chimeric read alignments:
        # first, to speed things up, extract all reads that are NOT properly paired
        samtools view -b -F 2 ~{bam} -o not_prop_pairs.bam
        samtools index not_prop_pairs.bam
        
        ~{util_dir}/extract_prelim_chimeric_genome_read_alignments.py \
           --star_bam not_prop_pairs.bam \
           --vif_full_tsv ~{prefix}.full.tsv \
           --output_bam ~{prefix}.genome_chimeric_evidence.bam
      
        # add evidence read stats
        ~{util_dir}/incorporate_read_alignment_stats.py \
          --supp_reads_bam ~{prefix}.genome_chimeric_evidence.bam \
          --vif_full_tsv ~{prefix}.full.tsv \
          --output ~{prefix}.full_read_stats.tsv

      
    >>>

    output {
        File full = "~{prefix}.full.tsv"
        File? full_filtered = "~{prefix}.full.filtered.tsv"
        File abridged = "~{prefix}.abridged.tsv"
        File? abridged_filtered = "~{prefix}.abridged.filtered.tsv"
        File abridged_detailed = "~{prefix}.abridged.detailed.tsv"
        File genome_chimeric_evidence_reads_bam = "~{prefix}.genome_chimeric_evidence.bam"
        File genome_chimeric_evidence_reads_bai = "~{prefix}.genome_chimeric_evidence.bam.bai"
        File full_read_stats = "~{prefix}.full_read_stats.tsv"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(viral_fasta, "GB") + size(chimeric_junction, "GB")*3 + size(bam, "GB")*2) + " HDD"
        docker: docker
        cpu: 1
        memory: "10GB"
    }
}

#task TopVirusCoverage {
#    input {
#        File chimeric_events
#        String util_dir
#        File bam
#        File bai
#        Int preemptible
#        String docker
#    }
#    String prefix = "vif"
#
#    command <<<
#        set -e
#
#        ~{util_dir}/plot_top_virus_coverage.Rscript \
#        --vif_report ~{chimeric_events} \
#        --bam ~{bam} \
#        --output_prefix ~{prefix}
#    >>>
#
#    output {
#        File read_counts_summary = "~{prefix}.virus_read_counts_summary.tsv"
#        File read_counts_image = "~{prefix}.virus_read_counts.png"
#        File read_counts_log_image = "~{prefix}.virus_read_counts_log.png"
#        Array[File] virus_images = glob("~{prefix}.virus_coverage_*.png")
#    }
#
#    runtime {
#        preemptible: preemptible
#        disks: "local-disk " + ceil(size(bam, "GB") + size(bai, "GB") + size(chimeric_events, "GB") + 2) + " HDD"
#        docker: docker
#        cpu: 1
#        memory: "2GB"
#    }
#}
#


#######################################
# Phase 2
#######################################

task ExtractChimericGenomicTargets {
    input {
        File fasta
        File viral_fasta
        File insertion_site_candidates_abridged
        String util_dir
        Int preemptible
        String docker
        String sample_id
    }

    # Set the file prefix 
    String prefix = sample_id + ".vif.extract"

    command <<<
        set -e

        ~{util_dir}/extract_chimeric_genomic_targets.py \
        --fasta ~{fasta} \
        --patch_db_fasta ~{viral_fasta} \
        --output_prefix ~{prefix} \
        --chim_events ~{insertion_site_candidates_abridged} \
        --pad_region_length 1000


      if [ -s ~{prefix}.fasta ]; then echo "true" > has_results; else echo "false" > has_results; fi
      
    >>>

    output {
        File fasta_extract = "~{prefix}.fasta"
        File gtf_extract = "~{prefix}.gtf"
        Boolean has_chimeric_targets = read_boolean("has_results")
    }
    
    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(viral_fasta, "GB")*2 + size(fasta, "GB")*2) + " HDD"
        docker: docker
        cpu: 1
        memory: "1GB"
    }
}

task CreateViralFasta {
    input {
        File fasta
        String picard
        Int preemptible
        String docker
    }

    command <<<
        mv ~{fasta} virus.fasta
        # we need fasta, index, and dict in same directory
        samtools faidx virus.fasta
        java -jar ~{picard} CreateSequenceDictionary -R virus.fasta
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: "1G"
        disks: "local-disk " + ceil(size(fasta, "GB")*2) + " HDD"
        preemptible: preemptible
        cpu: 1
    }

    output {
        File viral_fasta = "virus.fasta"
        File viral_fasta_index = "virus.fasta.fai"
        File viral_dict = "virus.dict"
    }
}

task ExtractViralReads {
    input {
        File bam
        File bai
        File fasta
        String picard
        Int preemptible
        String docker
    }

    command <<<
        samtools faidx ~{fasta}

        python <<CODE

        def parse_fai(path):
            values = set()
            with open(path, 'rt') as f:
                for line in f:
                    line = line.strip()
                    if line != '':
                        values.add(line.split('\t')[0])
            return values

        def to_txt(values, path):
            is_first = True
            with open(path, 'wt') as f:
                for val in values:
                    if not is_first:
                        f.write(' ')
                    f.write(val)
                    is_first = False

        to_txt(parse_fai('~{fasta}.fai'), 'extract.txt')
        CODE

        samtools view -b ~{bam} $(cat extract.txt) > tmp.bam
        samtools sort tmp.bam > virus.bam
        samtools index virus.bam
    >>>

    runtime {
        docker: docker
        bootDiskSizeGb: 12
        memory: "1G"
        disks: "local-disk " + ceil(size(bam, "GB")*2 + size(fasta, "GB")*2) + " HDD"
        preemptible: preemptible
        cpu: 1
    }

    output {
        File viral_bam = "virus.bam"
        File viral_bai = "virus.bam.bai"
    }
}

task ChimericContigEvidenceAnalyzer {
    input {
        File bam
        File bai
        File gtf

        String util_dir
        Int preemptible
        String docker
        Int min_reads
        String sample_id
    }

    String prefix = sample_id + ".vif"

    command <<<
        set -e

        ~{util_dir}/chimeric_contig_evidence_analyzer.py \
        --patch_db_bam ~{bam} \
        --patch_db_gtf ~{gtf} \
        --output_prefix ~{prefix}

        samtools index ~{prefix}.evidence.bam


      insertions_file="~{prefix}.evidence_counts.tsv"
      insertions_validated_file="insertions_validated.txt"
      num_insertions=$(cat ${insertions_file} | wc -l)
      if [[ ${num_insertions} -gt 1 ]]; then
        echo "true" > ${insertions_validated_file}
      else
        echo "false" > ${insertions_validated_file}
      fi
      
    >>>

    output {
        File evidence_counts = "~{prefix}.evidence_counts.tsv"
        File evidence_bam = "~{prefix}.evidence.bam"
        File evidence_bai = "~{prefix}.evidence.bam.bai"
        Boolean insertions_recovered = read_boolean("insertions_validated.txt")
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil( size(bam, "GB")*2 + 1) + " HDD"
        docker: docker
        cpu: 1
        memory: "1GB"
    }
}


task VirusReport {
    input {
        File bam
        File bai
        File viral_fasta
        File insertion_site_candidates
        Boolean remove_duplicates
        String util_dir
        Int preemptible
        String docker
        String memory
        String sample_id
        Int num_top_viruses = 10
    }
    Int max_coverage = 100

    String prefix = sample_id + ".vif"

    command <<<
        set -e

        ~{util_dir}/make_VIF_genome_abundance_plot.Rscript \
          --vif_report ~{insertion_site_candidates} \
          --title "Preliminary Genome Wide Abundance" \
          --output_png ~{prefix}.init.genome_plot.png

        bam=~{bam}

        ## restrict bam to only viruses
        samtools faidx ~{viral_fasta}
        awk '{printf("%s\t0\t%s\n",$1,$2);}' ~{viral_fasta}.fai  > viral_fasta.bed
        samtools view -b -L viral_fasta.bed ${bam} -o virus_only.bam
        samtools index virus_only.bam
        bam="virus_only.bam"
        
        if [ "~{remove_duplicates}" == "true" ]; then
            ~{util_dir}/bam_mark_duplicates.py -i ${bam}  -o dups.removed.bam -r
            samtools index dups.removed.bam
            bam="dups.removed.bam"
        fi

      
        # generates read_counts_summary and images
        ~{util_dir}/plot_top_virus_coverage.Rscript \
          --vif_report ~{insertion_site_candidates} \
          --bam ${bam} \
          --output_prefix ~{prefix}

        ~{util_dir}/create_insertion_site_inspector_js.py \
          --VIF_summary_tsv ~{prefix}.virus_read_counts_summary.tsv \
          --json_outfile ~{prefix}.virus.json

        # make bed for igvjs
        ~{util_dir}/create_igvjs_virus_bed.py \
            --summary ~{prefix}.virus_read_counts_summary.tsv \
            --output ~{prefix}.virus.bed \
            --num_top_viruses ~{num_top_viruses}

        # prep for making the report
        ~{util_dir}/bamsifter/bamsifter \
          -c ~{max_coverage} \
          -o ~{prefix}.virus.reads.bam \
          ${bam} 

        # IGV reports expects to find, __PREFIX__.fa, __PREFIX__.bed, __PREFIX__.reads.bam
        #ln -sf ~{viral_fasta} ~{prefix}.virus.fa
        ~{util_dir}/create_igvjs_virus_fa.py \
          ~{prefix}.virus.bed \
          ~{viral_fasta}  \
          ~{prefix}.virus.fa
      
        # generate the html
        ~{util_dir}/make_VIF_igvjs_html.py \
          --html_template ~{util_dir}/resources/igvjs_VIF.html \
          --fusions_json ~{prefix}.virus.json \
          --input_file_prefix ~{prefix}.virus \
          --html_output ~{prefix}.virus.html
    >>>

    output {
        File html = "~{prefix}.virus.html"
        File genome_abundance_plot = "~{prefix}.init.genome_plot.png"
        File virus_alignments_bam = "virus_only.bam"
        File virus_alignments_bai = "virus_only.bam.bai"
        File read_counts_summary = "~{prefix}.virus_read_counts_summary.tsv"
        File read_counts_image = "~{prefix}.virus_read_counts.png"
        File read_counts_log_image = "~{prefix}.virus_read_counts_log.png"
        Array[File] virus_images = glob("~{prefix}.virus_coverage_*.png")
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(bam, "GB")*2 + size(viral_fasta, "GB")*2 + 2) + " HDD"
        docker: docker
        cpu: 1
        memory: memory
    }
}

task SummaryReport {
    input {
        File init_counts
        File vif_counts
        File alignment_bam
        File alignment_bai
        File chim_targets_gtf
        File chim_targets_fasta
        File gtf
        Array[File] images
        String util_dir
        Int preemptible
        String docker
        String memory
        String sample_id
    }
    Int max_coverage = 100
    String prefix = sample_id + ".vif"
    String image_prefix = if(length(images)>0) then "--image " else ""
    command <<<
        set -e

        ~{util_dir}/refine_VIF_output.Rscript \
        --prelim_counts ~{init_counts} \
        --vif_counts ~{vif_counts} \
        --output ~{prefix}.refined.tsv

        ~{util_dir}/make_VIF_genome_abundance_plot.Rscript \
        --vif_report ~{prefix}.refined.tsv \
        --title "Genome Wide Abundance" \
        --output_png ~{prefix}.genome_plot.png

        ~{util_dir}/find_closest.py \
        -i ~{prefix}.refined.tsv \
        -o summary_results_tsv_with_genes.tsv \
        --gtf ~{gtf}

        ~{util_dir}/create_insertion_site_inspector_js.py \
            --VIF_summary_tsv summary_results_tsv_with_genes.tsv \
            --json_outfile igv.json

        # make bed for igvjs
        ~{util_dir}/region_gtf_to_bed.py \
            ~{chim_targets_gtf} \
            > ~{prefix}.bed

        # prep for making the report
        ~{util_dir}/bamsifter/bamsifter \
        -c ~{max_coverage} \
        -o ~{prefix}.reads.bam \
        ~{alignment_bam}

        # IGV reports expects to find, __PREFIX__.fa, __PREFIX__.bed, __PREFIX__.reads.bam
        ln -sf ~{chim_targets_fasta} ~{prefix}.fa

        ~{util_dir}/make_VIF_igvjs_html.py \
        --html_template ~{util_dir}/resources/igvjs_VIF.html \
        --fusions_json igv.json \
        --input_file_prefix ~{prefix} \
        --html_output ~{prefix}.html

        # generate the final report
        ~{util_dir}/add_to_html.py \
        --html ~{prefix}.html \
        --out ~{prefix}.html \
        --image ~{prefix}.genome_plot.png \
        ~{image_prefix}~{sep=' --image ' images}
    >>>

    output {
        File html = "~{prefix}.html"
        File refined_counts = "~{prefix}.refined.tsv"
        File genome_abundance_plot = "~{prefix}.genome_plot.png"
    }
    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil( size(images, "GB") + size(alignment_bam, "GB")*2 + size(init_counts, "GB") + size(vif_counts, "GB") + size(chim_targets_fasta,"GB")*2 + 2) + " HDD"
        docker: docker
        cpu: 1
        memory: "16GB"
    }
}





task ExtractEvidenceReads {
    input {
        File fastq1
        File? fastq2
        File orig_insertion_site_candidates
        String util_dir
        Int preemptible
        String docker
        String sample_id

    }
    String prefix = sample_id + ".vif.init"

    command <<<
        set -e
        
        # Combine the FastQ's 
        fastqs="~{fastq1} ~{fastq2}"

        # special case for tar of fastq files
        if [[ "~{fastq1}" == *.tar.gz ]] ; then
            mkdir fastq
            tar -I pigz -xvf ~{fastq1} -C fastq
            fastqs=$(find fastq -type f)
        fi

        ~{util_dir}/extract_insertion_evidence_reads.py \
            --fastqs $fastqs \
            --insertion_candidates ~{orig_insertion_site_candidates} \
            --out_prefix ~{prefix}

    >>>

    output {
        File left = "~{prefix}_1.fastq"
        File? right = "~{prefix}_2.fastq"
    }

    runtime {
        preemptible: preemptible
        disks: "local-disk " + ceil(size(fastq1, "GB") + size(orig_insertion_site_candidates, "GB")*3) + " HDD"
        docker: docker
        cpu: 1
        memory: "2GB"
    }
}

