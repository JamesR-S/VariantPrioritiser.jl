using TOML

Base.@kwdef struct ThresholdConfig
    frequency_cutoff::Float64 = 0.001
    homozygous_frequency_cutoff::Float64 = 0.01
    denovo_frequency_cutoff::Float64 = 0.0001
    spliceai_cutoff::Float64 = 0.5
    gq_cutoff::Float64 = 10.0
    gq_hom_cutoff::Float64 = 5.0
    mq_cutoff::Float64 = -10.0
    min_alt_reads::Int = 2
    min_denovo_read_fraction::Float64 = 0.1
    splice_acceptor_window::Int = 50
    splice_donor_window::Int = 10
    extended_cnv_splice_window::Int = 400
    gene_distance_cutoff::Int = 2000
    manta_gq_cutoff::Float64 = 20.0
    manta_min_total_alt_support::Int = 8
    manta_min_sr_alt_support::Int = 2
    manta_min_parent_alt_support::Int = 2
    manta_min_size::Int = 50
end

Base.@kwdef struct ResourceConfig
    gtf::Vector{String} = String[]
    omim_morbid::Vector{String} = String[]
    omim_mim2gene::Vector{String} = String[]
    omim_ensg_to_hgnc::Vector{String} = String[]
    panelapp_genes::Vector{String} = String[]
    sample_comments::Vector{String} = String[]
    imprinted_genes::Vector{String} = String[]
    excluded_variants_diag::Vector{String} = String[]
    excluded_genes_diag::Vector{String} = String[]
end

Base.@kwdef struct AppConfig
    thresholds::ThresholdConfig = ThresholdConfig()
    resources::ResourceConfig = ResourceConfig()
end

function load_config(path::Union{Nothing,String})
    config_path = isnothing(path) ? joinpath(@__DIR__, "..", "config", "defaults.toml") : path
    raw = TOML.parsefile(config_path)
    thresholds = get(raw, "thresholds", Dict{String,Any}())
    resources = get(raw, "resources", Dict{String,Any}())
    return AppConfig(
        thresholds=ThresholdConfig(
            frequency_cutoff=Float64(get(thresholds, "frequency_cutoff", 0.001)),
            homozygous_frequency_cutoff=Float64(get(thresholds, "homozygous_frequency_cutoff", 0.01)),
            denovo_frequency_cutoff=Float64(get(thresholds, "denovo_frequency_cutoff", 0.0001)),
            spliceai_cutoff=Float64(get(thresholds, "spliceai_cutoff", 0.5)),
            gq_cutoff=Float64(get(thresholds, "gq_cutoff", 10.0)),
            gq_hom_cutoff=Float64(get(thresholds, "gq_hom_cutoff", 5.0)),
            mq_cutoff=Float64(get(thresholds, "mq_cutoff", -10.0)),
            min_alt_reads=Int(get(thresholds, "min_alt_reads", 2)),
            min_denovo_read_fraction=Float64(get(thresholds, "min_denovo_read_fraction", 0.1)),
            splice_acceptor_window=Int(get(thresholds, "splice_acceptor_window", 50)),
            splice_donor_window=Int(get(thresholds, "splice_donor_window", 10)),
            extended_cnv_splice_window=Int(get(thresholds, "extended_cnv_splice_window", 400)),
            gene_distance_cutoff=Int(get(thresholds, "gene_distance_cutoff", 2000)),
            manta_gq_cutoff=Float64(get(thresholds, "manta_gq_cutoff", 20.0)),
            manta_min_total_alt_support=Int(get(thresholds, "manta_min_total_alt_support", 8)),
            manta_min_sr_alt_support=Int(get(thresholds, "manta_min_sr_alt_support", 2)),
            manta_min_parent_alt_support=Int(get(thresholds, "manta_min_parent_alt_support", 2)),
            manta_min_size=Int(get(thresholds, "manta_min_size", 50)),
        ),
        resources=ResourceConfig(
            gtf=String.(get(resources, "gtf", String[])),
            omim_morbid=String.(get(resources, "omim_morbid", String[])),
            omim_mim2gene=String.(get(resources, "omim_mim2gene", String[])),
            omim_ensg_to_hgnc=String.(get(resources, "omim_ensg_to_hgnc", String[])),
            panelapp_genes=String.(get(resources, "panelapp_genes", String[])),
            sample_comments=String.(get(resources, "sample_comments", String[])),
            imprinted_genes=String.(get(resources, "imprinted_genes", String[])),
            excluded_variants_diag=String.(get(resources, "excluded_variants_diag", String[])),
            excluded_genes_diag=String.(get(resources, "excluded_genes_diag", String[])),
        ),
    )
end
