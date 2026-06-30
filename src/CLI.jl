Base.@kwdef mutable struct RunOptions
    input::String = ""
    samples::Vector{String} = String[]
    parent1::Union{Nothing,String} = nothing
    parent2::Union{Nothing,String} = nothing
    parent1_affected::Bool = false
    parent2_affected::Bool = false
    affected::Vector{String} = String[]
    unaffected::Set{String} = Set{String}()
    extra_inputs::Vector{String} = String[]
    output::Union{Nothing,String} = nothing
    output_format::Symbol = :tsv
    pipeline_prefix::Union{Nothing,String} = nothing
    config_path::Union{Nothing,String} = nothing
    freq_cutoff::Union{Nothing,Float64} = nothing
    gq_cutoff::Union{Nothing,Float64} = nothing
    gq_hom_cutoff::Union{Nothing,Float64} = nothing
    mq_cutoff::Union{Nothing,Float64} = nothing
    min_alt_reads::Union{Nothing,Int} = nothing
    splice::Bool = false
    shared::Bool = false
    recessive_only::Bool = false
    denovo_only::Bool = false
    lowq::Bool = false
    debug::Bool = false
    discovery::Bool = false
    highlight_genes::Set{String} = Set{String}()
    only_genes::Set{String} = Set{String}()
end

Base.@kwdef struct FamilySpec
    parent1::Union{Nothing,String} = nothing
    parent2::Union{Nothing,String} = nothing
    parent1_affected::Bool = false
    parent2_affected::Bool = false
    affected::Vector{String} = String[]
    unaffected::Set{String} = Set{String}()
    shared::Bool = false
end

function parse_cli(args::Vector{String})
    options = RunOptions()
    is_unaffected_mode = false
    positional = String[]
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "-splice"
            options.splice = true
        elseif arg == "-shared"
            options.shared = true
        elseif arg == "-recessive"
            options.recessive_only = true
        elseif arg == "-denovo"
            options.denovo_only = true
        elseif arg == "-lowq"
            options.lowq = true
        elseif arg == "-debug"
            options.debug = true
        elseif arg == "-discovery"
            options.discovery = true
        elseif arg == "-x"
            options.output_format = :tsv
        elseif arg == "-s" || arg == "-a"
            options.output_format = :tsv
        elseif arg == "--html"
            options.output_format = :html
        elseif arg == "--tsv"
            options.output_format = :tsv
        elseif arg == "-f"
            i += 1
            options.freq_cutoff = parse(Float64, args[i])
        elseif arg == "-q"
            i += 1
            options.gq_cutoff = parse(Float64, args[i])
        elseif arg == "-hq"
            i += 1
            options.gq_hom_cutoff = parse(Float64, args[i])
        elseif arg == "-mq"
            i += 1
            options.mq_cutoff = parse(Float64, args[i])
        elseif arg == "-minAlt"
            i += 1
            options.min_alt_reads = parse(Int, args[i])
        elseif arg == "-extra"
            i += 1
            push!(options.extra_inputs, args[i])
        elseif arg == "-highlight"
            i += 1
            union!(options.highlight_genes, split(args[i], ','))
        elseif arg == "-onlyGenes"
            i += 1
            union!(options.only_genes, split(args[i], ','))
        elseif arg == "-unaffected"
            is_unaffected_mode = true
        elseif arg == "--config"
            i += 1
            options.config_path = args[i]
        elseif arg == "--parent1"
            i += 1
            options.parent1 = args[i]
        elseif arg == "--parent2"
            i += 1
            options.parent2 = args[i]
        elseif arg == "--affected-parent1"
            options.parent1_affected = true
        elseif arg == "--affected-parent2"
            options.parent2_affected = true
        elseif arg == "--affected"
            i += 1
            append!(options.affected, filter(!isempty, split(args[i], ',')))
        elseif arg == "--unaffected-list"
            i += 1
            union!(options.unaffected, filter(!isempty, split(args[i], ',')))
        elseif arg == "--out"
            i += 1
            options.output = args[i]
        elseif arg == "--pipeline-prefix"
            i += 1
            options.pipeline_prefix = args[i]
        elseif startswith(arg, "-")
            throw(ArgumentError("Unsupported option: $arg"))
        else
            push!(positional, arg)
            if is_unaffected_mode
                push!(options.unaffected, arg)
            end
        end
        i += 1
    end

    isempty(positional) && throw(ArgumentError("Missing input VCF/TSV path"))
    options.input = positional[1]
    if length(positional) > 1
        options.samples = positional[2:end]
    end
    if isnothing(options.pipeline_prefix)
        options.pipeline_prefix = infer_pipeline_prefix(options.input, options.output)
    end
    return options
end

function infer_pipeline_prefix(paths::Vararg{Union{Nothing,String}})
    for path in paths
        if isnothing(path)
            continue
        end
        match_obj = match(r"(r\d{2})", path)
        if match_obj !== nothing
            return match_obj.captures[1]
        end
    end
    return "r04"
end

function resolve_family(options::RunOptions)
    family = FamilySpec(shared=options.shared, unaffected=copy(options.unaffected))
    if !isnothing(options.parent1) || !isnothing(options.parent2) || !isempty(options.affected)
        return FamilySpec(
            parent1=options.parent1,
            parent2=options.parent2,
            parent1_affected=options.parent1_affected,
            parent2_affected=options.parent2_affected,
            affected=copy(options.affected),
            unaffected=copy(options.unaffected),
            shared=options.shared,
        )
    end
    if !isempty(options.samples)
        if !options.shared && length(options.samples) >= 2
            family = FamilySpec(
                parent1=options.samples[1],
                parent2=options.samples[2],
                parent1_affected=options.parent1_affected,
                parent2_affected=options.parent2_affected,
                affected=length(options.samples) > 2 ? options.samples[3:end] : String[],
                unaffected=copy(options.unaffected),
                shared=options.shared,
            )
        else
            family = FamilySpec(
                affected=copy(options.samples),
                unaffected=copy(options.unaffected),
                shared=options.shared,
            )
        end
    end
    return family
end

function resolve_output_path(options::RunOptions)
    if !isnothing(options.output)
        return options.output
    end
    stem = replace(basename(options.input), r"\.(vcf|vcf\.gz|tsv)$" => "")
    extension = options.output_format == :html ? "html" : "tsv"
    return "$(stem).prioritised.$extension"
end
