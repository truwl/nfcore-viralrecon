version 1.0

workflow viralrecon {
	input{
		File samplesheet
		String? platform
		String? protocol
		String outdir
		String? email
		String? genome
		File? fasta
		String? gff
		String? bowtie2_index
		String? primer_bed
		String? primer_fasta
		String? primer_set
		Float? primer_set_version
		String primer_left_suffix = "_LEFT"
		String primer_right_suffix = "_RIGHT"
		Boolean? save_reference
		String? fastq_dir
		String? fast5_dir
		String? sequencing_summary
		Int min_barcode_reads = 100
		Int min_guppyplex_reads = 10
		String artic_minion_caller = "nanopolish"
		String artic_minion_aligner = "minimap2"
		String? artic_scheme
		String? artic_minion_medaka_model
		Boolean? skip_pycoqc
		Boolean? skip_nanoplot
		String? nextclade_dataset
		String? nextclade_dataset_name
		String? nextclade_dataset_reference
		String? nextclade_dataset_tag
		Int asciigenome_read_depth = 50
		Int asciigenome_window_size = 50
		String? multiqc_title
		String? multiqc_config
		String max_multiqc_email_size = "25.MB"
		Boolean? skip_mosdepth
		Boolean? skip_pangolin
		Boolean? skip_nextclade
		Boolean? skip_asciigenome
		Boolean? skip_variants_quast
		Boolean? skip_variants_long_table
		Boolean? skip_multiqc
		String kraken2_db = "s3://nf-core-awsmegatests/viralrecon/input_data/kraken2_human.tar.gz"
		String kraken2_db_name = "human"
		Boolean? kraken2_variants_host_filter
		Boolean kraken2_assembly_host_filter = true
		Boolean? save_trimmed_fail
		Boolean? skip_fastqc
		Boolean? skip_kraken2
		Boolean? skip_fastp
		Boolean? skip_cutadapt
		String? variant_caller
		String consensus_caller = "bcftools"
		Int min_mapped_reads = 1000
		Boolean? ivar_trim_noprimer
		Int? ivar_trim_offset
		Boolean? filter_duplicates
		Boolean? save_unaligned
		Boolean? save_mpileup
		Boolean? skip_ivar_trim
		Boolean skip_markduplicates = true
		Boolean? skip_picard_metrics
		Boolean? skip_snpeff
		Boolean? skip_consensus_plots
		Boolean? skip_consensus
		Boolean? skip_variants
		String assemblers = "spades"
		String spades_mode = "rnaviral"
		String? spades_hmm
		String? blast_db
		Boolean? skip_bandage
		Boolean? skip_blast
		Boolean? skip_abacas
		Boolean? skip_plasmidid
		Boolean? skip_assembly_quast
		Boolean? skip_assembly
		Boolean? help
		String publish_dir_mode = "copy"
		String? email_on_fail
		Boolean? plaintext_email
		Boolean? monochrome_logs
		String tracedir = "./results/pipeline_info"
		Boolean? enable_conda
		Boolean validate_params = true
		Boolean? show_hidden_params
		Int max_cpus = 16
		String max_memory = "128.GB"
		String max_time = "240.h"
		String custom_config_version = "master"
		String custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/master"
		String? config_profile_name
		String? config_profile_description
		String? config_profile_contact
		String? config_profile_url

	}

	call make_uuid as mkuuid {}
	call touch_uuid as thuuid {
		input:
			outputbucket = mkuuid.uuid
	}
	call run_nfcoretask as nfcoretask {
		input:
			samplesheet = samplesheet,
			platform = platform,
			protocol = protocol,
			outdir = outdir,
			email = email,
			genome = genome,
			fasta = fasta,
			gff = gff,
			bowtie2_index = bowtie2_index,
			primer_bed = primer_bed,
			primer_fasta = primer_fasta,
			primer_set = primer_set,
			primer_set_version = primer_set_version,
			primer_left_suffix = primer_left_suffix,
			primer_right_suffix = primer_right_suffix,
			save_reference = save_reference,
			fastq_dir = fastq_dir,
			fast5_dir = fast5_dir,
			sequencing_summary = sequencing_summary,
			min_barcode_reads = min_barcode_reads,
			min_guppyplex_reads = min_guppyplex_reads,
			artic_minion_caller = artic_minion_caller,
			artic_minion_aligner = artic_minion_aligner,
			artic_scheme = artic_scheme,
			artic_minion_medaka_model = artic_minion_medaka_model,
			skip_pycoqc = skip_pycoqc,
			skip_nanoplot = skip_nanoplot,
			nextclade_dataset = nextclade_dataset,
			nextclade_dataset_name = nextclade_dataset_name,
			nextclade_dataset_reference = nextclade_dataset_reference,
			nextclade_dataset_tag = nextclade_dataset_tag,
			asciigenome_read_depth = asciigenome_read_depth,
			asciigenome_window_size = asciigenome_window_size,
			multiqc_title = multiqc_title,
			multiqc_config = multiqc_config,
			max_multiqc_email_size = max_multiqc_email_size,
			skip_mosdepth = skip_mosdepth,
			skip_pangolin = skip_pangolin,
			skip_nextclade = skip_nextclade,
			skip_asciigenome = skip_asciigenome,
			skip_variants_quast = skip_variants_quast,
			skip_variants_long_table = skip_variants_long_table,
			skip_multiqc = skip_multiqc,
			kraken2_db = kraken2_db,
			kraken2_db_name = kraken2_db_name,
			kraken2_variants_host_filter = kraken2_variants_host_filter,
			kraken2_assembly_host_filter = kraken2_assembly_host_filter,
			save_trimmed_fail = save_trimmed_fail,
			skip_fastqc = skip_fastqc,
			skip_kraken2 = skip_kraken2,
			skip_fastp = skip_fastp,
			skip_cutadapt = skip_cutadapt,
			variant_caller = variant_caller,
			consensus_caller = consensus_caller,
			min_mapped_reads = min_mapped_reads,
			ivar_trim_noprimer = ivar_trim_noprimer,
			ivar_trim_offset = ivar_trim_offset,
			filter_duplicates = filter_duplicates,
			save_unaligned = save_unaligned,
			save_mpileup = save_mpileup,
			skip_ivar_trim = skip_ivar_trim,
			skip_markduplicates = skip_markduplicates,
			skip_picard_metrics = skip_picard_metrics,
			skip_snpeff = skip_snpeff,
			skip_consensus_plots = skip_consensus_plots,
			skip_consensus = skip_consensus,
			skip_variants = skip_variants,
			assemblers = assemblers,
			spades_mode = spades_mode,
			spades_hmm = spades_hmm,
			blast_db = blast_db,
			skip_bandage = skip_bandage,
			skip_blast = skip_blast,
			skip_abacas = skip_abacas,
			skip_plasmidid = skip_plasmidid,
			skip_assembly_quast = skip_assembly_quast,
			skip_assembly = skip_assembly,
			help = help,
			publish_dir_mode = publish_dir_mode,
			email_on_fail = email_on_fail,
			plaintext_email = plaintext_email,
			monochrome_logs = monochrome_logs,
			tracedir = tracedir,
			enable_conda = enable_conda,
			validate_params = validate_params,
			show_hidden_params = show_hidden_params,
			max_cpus = max_cpus,
			max_memory = max_memory,
			max_time = max_time,
			custom_config_version = custom_config_version,
			custom_config_base = custom_config_base,
			config_profile_name = config_profile_name,
			config_profile_description = config_profile_description,
			config_profile_contact = config_profile_contact,
			config_profile_url = config_profile_url,
			outputbucket = thuuid.touchedbucket
            }
		output {
			Array[File] results = nfcoretask.results
		}
	}
task make_uuid {
	meta {
		volatile: true
}

command <<<
        python <<CODE
        import uuid
        print("gs://truwl-internal-inputs/nf-viralrecon/{}".format(str(uuid.uuid4())))
        CODE
>>>

  output {
    String uuid = read_string(stdout())
  }
  
  runtime {
    docker: "python:3.8.12-buster"
  }
}

task touch_uuid {
    input {
        String outputbucket
    }

    command <<<
        echo "sentinel" > sentinelfile
        gsutil cp sentinelfile ~{outputbucket}/sentinelfile
    >>>

    output {
        String touchedbucket = outputbucket
    }

    runtime {
        docker: "google/cloud-sdk:latest"
    }
}

task fetch_results {
    input {
        String outputbucket
        File execution_trace
    }
    command <<<
        cat ~{execution_trace}
        echo ~{outputbucket}
        mkdir -p ./resultsdir
        gsutil cp -R ~{outputbucket} ./resultsdir
    >>>
    output {
        Array[File] results = glob("resultsdir/*")
    }
    runtime {
        docker: "google/cloud-sdk:latest"
    }
}

task run_nfcoretask {
    input {
        String outputbucket
		File samplesheet
		String? platform
		String? protocol
		String outdir
		String? email
		String? genome
		File? fasta
		String? gff
		String? bowtie2_index
		String? primer_bed
		String? primer_fasta
		String? primer_set
		Float? primer_set_version
		String primer_left_suffix = "_LEFT"
		String primer_right_suffix = "_RIGHT"
		Boolean? save_reference
		String? fastq_dir
		String? fast5_dir
		String? sequencing_summary
		Int min_barcode_reads = 100
		Int min_guppyplex_reads = 10
		String artic_minion_caller = "nanopolish"
		String artic_minion_aligner = "minimap2"
		String? artic_scheme
		String? artic_minion_medaka_model
		Boolean? skip_pycoqc
		Boolean? skip_nanoplot
		String? nextclade_dataset
		String? nextclade_dataset_name
		String? nextclade_dataset_reference
		String? nextclade_dataset_tag
		Int asciigenome_read_depth = 50
		Int asciigenome_window_size = 50
		String? multiqc_title
		String? multiqc_config
		String max_multiqc_email_size = "25.MB"
		Boolean? skip_mosdepth
		Boolean? skip_pangolin
		Boolean? skip_nextclade
		Boolean? skip_asciigenome
		Boolean? skip_variants_quast
		Boolean? skip_variants_long_table
		Boolean? skip_multiqc
		String kraken2_db = "s3://nf-core-awsmegatests/viralrecon/input_data/kraken2_human.tar.gz"
		String kraken2_db_name = "human"
		Boolean? kraken2_variants_host_filter
		Boolean kraken2_assembly_host_filter = true
		Boolean? save_trimmed_fail
		Boolean? skip_fastqc
		Boolean? skip_kraken2
		Boolean? skip_fastp
		Boolean? skip_cutadapt
		String? variant_caller
		String consensus_caller = "bcftools"
		Int min_mapped_reads = 1000
		Boolean? ivar_trim_noprimer
		Int? ivar_trim_offset
		Boolean? filter_duplicates
		Boolean? save_unaligned
		Boolean? save_mpileup
		Boolean? skip_ivar_trim
		Boolean skip_markduplicates = true
		Boolean? skip_picard_metrics
		Boolean? skip_snpeff
		Boolean? skip_consensus_plots
		Boolean? skip_consensus
		Boolean? skip_variants
		String assemblers = "spades"
		String spades_mode = "rnaviral"
		String? spades_hmm
		String? blast_db
		Boolean? skip_bandage
		Boolean? skip_blast
		Boolean? skip_abacas
		Boolean? skip_plasmidid
		Boolean? skip_assembly_quast
		Boolean? skip_assembly
		Boolean? help
		String publish_dir_mode = "copy"
		String? email_on_fail
		Boolean? plaintext_email
		Boolean? monochrome_logs
		String tracedir = "./results/pipeline_info"
		Boolean? enable_conda
		Boolean validate_params = true
		Boolean? show_hidden_params
		Int max_cpus = 16
		String max_memory = "128.GB"
		String max_time = "240.h"
		String custom_config_version = "master"
		String custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/master"
		String? config_profile_name
		String? config_profile_description
		String? config_profile_contact
		String? config_profile_url

	}
	command <<<
		export NXF_VER=21.10.5
		export NXF_MODE=google
		echo ~{outputbucket}
		/nextflow -c /truwl.nf.config run /viralrecon-2.4.1  -profile truwl  --input ~{samplesheet} 	~{"--samplesheet " + samplesheet}	~{"--platform " + platform}	~{"--protocol " + protocol}	~{"--outdir " + outdir}	~{"--email " + email}	~{"--genome " + genome}	~{"--fasta " + fasta}	~{"--gff " + gff}	~{"--bowtie2_index " + bowtie2_index}	~{"--primer_bed " + primer_bed}	~{"--primer_fasta " + primer_fasta}	~{"--primer_set " + primer_set}	~{"--primer_set_version " + primer_set_version}	~{"--primer_left_suffix " + primer_left_suffix}	~{"--primer_right_suffix " + primer_right_suffix}	~{true="--save_reference  " false="" save_reference}	~{"--fastq_dir " + fastq_dir}	~{"--fast5_dir " + fast5_dir}	~{"--sequencing_summary " + sequencing_summary}	~{"--min_barcode_reads " + min_barcode_reads}	~{"--min_guppyplex_reads " + min_guppyplex_reads}	~{"--artic_minion_caller " + artic_minion_caller}	~{"--artic_minion_aligner " + artic_minion_aligner}	~{"--artic_scheme " + artic_scheme}	~{"--artic_minion_medaka_model " + artic_minion_medaka_model}	~{true="--skip_pycoqc  " false="" skip_pycoqc}	~{true="--skip_nanoplot  " false="" skip_nanoplot}	~{"--nextclade_dataset " + nextclade_dataset}	~{"--nextclade_dataset_name " + nextclade_dataset_name}	~{"--nextclade_dataset_reference " + nextclade_dataset_reference}	~{"--nextclade_dataset_tag " + nextclade_dataset_tag}	~{"--asciigenome_read_depth " + asciigenome_read_depth}	~{"--asciigenome_window_size " + asciigenome_window_size}	~{"--multiqc_title " + multiqc_title}	~{"--multiqc_config " + multiqc_config}	~{"--max_multiqc_email_size " + max_multiqc_email_size}	~{true="--skip_mosdepth  " false="" skip_mosdepth}	~{true="--skip_pangolin  " false="" skip_pangolin}	~{true="--skip_nextclade  " false="" skip_nextclade}	~{true="--skip_asciigenome  " false="" skip_asciigenome}	~{true="--skip_variants_quast  " false="" skip_variants_quast}	~{true="--skip_variants_long_table  " false="" skip_variants_long_table}	~{true="--skip_multiqc  " false="" skip_multiqc}	~{"--kraken2_db " + kraken2_db}	~{"--kraken2_db_name " + kraken2_db_name}	~{true="--kraken2_variants_host_filter  " false="" kraken2_variants_host_filter}	~{true="--kraken2_assembly_host_filter  " false="" kraken2_assembly_host_filter}	~{true="--save_trimmed_fail  " false="" save_trimmed_fail}	~{true="--skip_fastqc  " false="" skip_fastqc}	~{true="--skip_kraken2  " false="" skip_kraken2}	~{true="--skip_fastp  " false="" skip_fastp}	~{true="--skip_cutadapt  " false="" skip_cutadapt}	~{"--variant_caller " + variant_caller}	~{"--consensus_caller " + consensus_caller}	~{"--min_mapped_reads " + min_mapped_reads}	~{true="--ivar_trim_noprimer  " false="" ivar_trim_noprimer}	~{"--ivar_trim_offset " + ivar_trim_offset}	~{true="--filter_duplicates  " false="" filter_duplicates}	~{true="--save_unaligned  " false="" save_unaligned}	~{true="--save_mpileup  " false="" save_mpileup}	~{true="--skip_ivar_trim  " false="" skip_ivar_trim}	~{true="--skip_markduplicates  " false="" skip_markduplicates}	~{true="--skip_picard_metrics  " false="" skip_picard_metrics}	~{true="--skip_snpeff  " false="" skip_snpeff}	~{true="--skip_consensus_plots  " false="" skip_consensus_plots}	~{true="--skip_consensus  " false="" skip_consensus}	~{true="--skip_variants  " false="" skip_variants}	~{"--assemblers " + assemblers}	~{"--spades_mode " + spades_mode}	~{"--spades_hmm " + spades_hmm}	~{"--blast_db " + blast_db}	~{true="--skip_bandage  " false="" skip_bandage}	~{true="--skip_blast  " false="" skip_blast}	~{true="--skip_abacas  " false="" skip_abacas}	~{true="--skip_plasmidid  " false="" skip_plasmidid}	~{true="--skip_assembly_quast  " false="" skip_assembly_quast}	~{true="--skip_assembly  " false="" skip_assembly}	~{true="--help  " false="" help}	~{"--publish_dir_mode " + publish_dir_mode}	~{"--email_on_fail " + email_on_fail}	~{true="--plaintext_email  " false="" plaintext_email}	~{true="--monochrome_logs  " false="" monochrome_logs}	~{"--tracedir " + tracedir}	~{true="--enable_conda  " false="" enable_conda}	~{true="--validate_params  " false="" validate_params}	~{true="--show_hidden_params  " false="" show_hidden_params}	~{"--max_cpus " + max_cpus}	~{"--max_memory " + max_memory}	~{"--max_time " + max_time}	~{"--custom_config_version " + custom_config_version}	~{"--custom_config_base " + custom_config_base}	~{"--config_profile_name " + config_profile_name}	~{"--config_profile_description " + config_profile_description}	~{"--config_profile_contact " + config_profile_contact}	~{"--config_profile_url " + config_profile_url}	-w ~{outputbucket}
	>>>
        
    output {
        File execution_trace = "pipeline_execution_trace.txt"
        Array[File] results = glob("results/*/*html")
    }
    runtime {
        docker: "truwl/nfcore-viralrecon:2.4.1_0.1.0"
        memory: "2 GB"
        cpu: 1
    }
}
    