module VariantPrioritiser

include("Config.jl")
include("CLI.jl")
include("VCFLoader.jl")
include("MantaLoader.jl")
include("FamilyInference.jl")
include("Filtering.jl")
include("Report.jl")

export main

function main(args::Vector{String}=copy(ARGS))
    options = parse_cli(args)
    ensure_primary_input_exists(options.input)
    config = load_config(options.config_path)
    family = resolve_family(options)
    extra_inputs = collect_resolved_extra_inputs(options)
    headers, sample_names = scan_input_metadata(options.input, extra_inputs; debug=options.debug)
    family = infer_family(family, sample_names, options.pipeline_prefix, batch_root(options.input, options.pipeline_prefix))
    filtered_rows = prioritise_inputs(options, config, family, extra_inputs)
    output_path = resolve_output_path(options)
    write_report(output_path, filtered_rows, headers, family, options, config, sample_names)
    println(output_path)
    return nothing
end

function batch_root(input_path::AbstractString, pipeline_prefix::Union{Nothing,String})
    prefix = something(pipeline_prefix, infer_pipeline_prefix(input_path))
    input_dir = dirname(input_path)
    if basename(input_dir) == string(prefix, "_vep")
        return dirname(input_dir)
    end
    return input_dir
end

function collect_resolved_extra_inputs(options::RunOptions)
    extras = copy(options.extra_inputs)
    append!(extras, discover_manta_inputs(options.input, options.pipeline_prefix))
    return unique(extras)
end

function ensure_primary_input_exists(path::AbstractString)
    isfile(path) && return
    if occursin('*', path) || occursin('?', path)
        throw(ArgumentError("Primary input VCF not found: $path. This looks like an unexpanded shell glob; resolve it to a real file path before running Julia."))
    end
    throw(ArgumentError("Primary input VCF not found: $path"))
end

function warn_skip_missing_optional(path::AbstractString)
    @warn "Skipping missing optional input" path=path
    return nothing
end

function prioritise_inputs(options::RunOptions, config::AppConfig, family::FamilySpec, resolved_extra_inputs::Vector{String})
    context = prioritisation_context(options, config)
    filtered = Vector{Dict{String,Any}}()
    for path in [options.input; resolved_extra_inputs]
        if path != options.input && !isfile(path)
            warn_skip_missing_optional(path)
            continue
        end
        if is_manta_vcf(path)
            _, _, stream_state = prepare_manta_stream(path)
            stream_manta_rows(stream_state, row -> begin
                prioritise_manta_row!(filtered, row, family, options, context.thresholds)
            end; debug=options.debug)
        elseif endswith(lowercase(path), ".vcf") || endswith(lowercase(path), ".vcf.gz")
            _, _, stream_state = prepare_vep_stream(path)
            stream_vep_rows(stream_state, row -> begin
                prioritise_row!(filtered, row, family, options, context.thresholds, context.freq_cutoff, context.gq_cutoff, context.gq_hom_cutoff, context.mq_cutoff, context.denovocnn_calls, context.sample_sexes)
            end; debug=options.debug)
        else
            rows, _, _ = load_input(path; debug=options.debug)
            for row in rows
                prioritise_row!(filtered, row, family, options, context.thresholds, context.freq_cutoff, context.gq_cutoff, context.gq_hom_cutoff, context.mq_cutoff, context.denovocnn_calls)
            end
        end
    end
    return postprocess_prioritised_rows(filtered, family)
end

end
