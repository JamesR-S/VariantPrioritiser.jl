function discover_manta_inputs(input_path::AbstractString, pipeline_prefix::Union{Nothing,String})
    prefix = something(pipeline_prefix, "r04")
    input_stem = replace(basename(String(input_path)), r"\.(vcf|vcf\.gz|tsv)$" => "")
    search_roots = unique([
        dirname(input_path),
        dirname(dirname(input_path)),
    ])
    candidates = String[]
    for root in search_roots
        manta_dir = joinpath(root, "$(prefix)_manta")
        isdir(manta_dir) || continue
        append!(candidates, filter(path -> is_matching_annotated_manta_vcf(path, input_stem), readdir(manta_dir; join=true)))
    end
    return sort(unique(candidates))
end

function is_manta_vcf(path::AbstractString)
    lower = lowercase(path)
    return occursin(".sv.vcf", lower) || occursin("/manta/", lower) || occursin("_manta/", lower)
end

function is_annotated_manta_vcf(path::AbstractString)
    lower = lowercase(path)
    return endswith(lower, "vep_annotated.sv.vcf.gz") || endswith(lower, "vep_annotated.sv.vcf")
end

function is_matching_annotated_manta_vcf(path::AbstractString, input_stem::AbstractString)
    is_annotated_manta_vcf(path) || return false
    candidate_stem = replace(basename(String(path)), r"\.SV\.(vcf|vcf\.gz)$" => "")
    return candidate_stem == input_stem
end

function prepare_manta_stream(path::AbstractString)
    io = open_text(path)
    try
        metadata = String[]
        header_fields = String[]
        while !eof(io)
            line = chomp(readline(io))
            if startswith(line, "##")
                push!(metadata, line)
            elseif startswith(line, "#CHROM")
                header_fields = split(line[2:end], '\t')
                break
            end
        end
        isempty(header_fields) && throw(ArgumentError("VCF header not found in $path"))
        sample_names = length(header_fields) > 9 ? header_fields[10:end] : String[]
        csq_headers = try
            parse_csq_headers(metadata)
        catch err
            if err isa ArgumentError
                String[]
            else
                rethrow()
            end
        end
        headers = build_manta_headers(sample_names, csq_headers)
        return headers, sample_names, (io=io, sample_names=sample_names, csq_headers=csq_headers, path=String(path))
    catch
        close(io)
        rethrow()
    end
end

function stream_manta_rows(stream_state, callback::Function; debug::Bool=false)
    io = stream_state.io
    try
        isempty(stream_state.csq_headers) && return
        while !eof(io)
            line = chomp(readline(io))
            isempty(line) && continue
            startswith(line, "#") && continue
            fields = split(line, '\t')
            length(fields) < 10 && continue
            for row in normalise_manta_record(fields, stream_state.sample_names, stream_state.csq_headers, stream_state.path; debug=debug)
                callback(row)
            end
        end
    finally
        close(io)
    end
end

function build_manta_headers(sample_names::Vector{<:AbstractString}, csq_headers::Vector{<:AbstractString})
    headers = String[
        "candidateCategory",
        "gene",
        "geneId",
        "transcript",
        "Consequence",
        "IMPACT",
        "svType",
        "SVLEN",
        "gDNAstart",
        "gDNAend",
        "gNomen",
        "varLocation",
        "codingEffect",
        "inheritance_model",
        "Filter (VCF)",
        "Quality (VCF)",
        "MANE",
        "MANE_SELECT",
        "MANE_PLUS_CLINICAL",
    ]
    append!(headers, String.(csq_headers))
    for sample in sample_names
        append!(headers, ["GT ($sample)", "FT ($sample)", "GQ ($sample)", "PR ($sample)", "SR ($sample)"])
    end
    return unique(headers)
end

function normalise_manta_record(fields::AbstractVector{<:AbstractString}, sample_names::AbstractVector{<:AbstractString}, csq_headers::AbstractVector{<:AbstractString}, filename::AbstractString; debug::Bool=false)
    chrom, pos, _, ref, alt, qual, filter_field, info_field, format_field = fields[1:9]
    sample_fields = fields[10:end]
    info_map = parse_info_field(info_field)
    haskey(info_map, "CSQ") || return Dict{String,Any}[]
    sv_type = normalise_svtype(info_map, alt)
    end_pos = string(get(info_map, "END", pos))
    sv_len = normalise_svlen(get(info_map, "SVLEN", ""))
    format_keys = split(format_field, ':')
    rows = Dict{String,Any}[]
    for csq_entry in split(String(info_map["CSQ"]), ',')
        csq_values = split(csq_entry, '|', keepempty=true)
        length(csq_values) == length(csq_headers) || continue
        csq_map = Dict{String,String}(zip(csq_headers, csq_values))
        gene = get(csq_map, "SYMBOL", "")
        isempty(gene) && continue
        row = Dict{String,Any}()
        consequence = get(csq_map, "Consequence", "")
        var_location, coding_effect = classify_manta_consequence(consequence, sv_type)
        row["candidateCategory"] = ""
        row["source"] = "manta"
        row["chrom"] = String(chrom)
        row["inputPos"] = String(pos)
        row["inputRef"] = String(ref)
        row["inputAlt"] = String(alt)
        row["svType"] = sv_type
        row["SVLEN"] = sv_len
        row["gene"] = gene
        row["geneId"] = replace(get(csq_map, "HGNC_ID", ""), "HGNC:" => "")
        row["biotype"] = get(csq_map, "BIOTYPE", "")
        row["transcript"] = get(csq_map, "Feature", "")
        row["Consequence"] = consequence
        row["IMPACT"] = get(csq_map, "IMPACT", "")
        row["DISTANCE"] = get(csq_map, "DISTANCE", "")
        row["varLocation"] = var_location
        row["codingEffect"] = coding_effect
        row["gDNAstart"] = String(pos)
        row["gDNAend"] = end_pos
        row["gNomen"] = normalise_manta_gnomen(chrom, pos, end_pos, sv_type)
        row["cNomen"] = suffix_after_colon(get(csq_map, "HGVSc", ""))
        row["pNomen"] = suffix_after_colon(get(csq_map, "HGVSp", ""))
        row["MANE"] = get(csq_map, "MANE", "")
        row["MANE_SELECT"] = get(csq_map, "MANE_SELECT", "")
        row["MANE_PLUS_CLINICAL"] = get(csq_map, "MANE_PLUS_CLINICAL", "")
        row["Filter (VCF)"] = String(filter_field)
        row["Quality (VCF)"] = String(qual)
        row["inheritance_model"] = ""
        apply_manta_sample_fields!(row, sample_names, sample_fields, format_keys)
        push!(rows, row)
    end
    return rows
end

function normalise_svtype(info_map::Dict{String,Any}, alt::AbstractString)
    raw = uppercase(String(get(info_map, "SVTYPE", replace(String(alt), ['<', '>'] => ""))))
    raw == "DUP:TANDEM" && return "DUP"
    return raw
end

function normalise_svlen(value::Any)
    text = strip(string(value))
    isempty(text) && return ""
    first_value = split(text, ',')[1]
    parsed = tryparse(Int, first_value)
    parsed === nothing && return first_value
    return string(abs(parsed))
end

function normalise_manta_gnomen(chrom::AbstractString, start_pos::AbstractString, end_pos::AbstractString, sv_type::AbstractString)
    return string(normalise_chromosome(chrom), ":", start_pos, "-", end_pos, " ", sv_type)
end

function classify_manta_consequence(consequence::AbstractString, sv_type::AbstractString)
    terms = split(lowercase(consequence), '&')
    if any(term -> occursin("coding_sequence_variant", term), terms)
        return "exon", "coding_sequence_variant"
    elseif any(term -> occursin("splice_", term), terms)
        return "splice", "splice_variant"
    elseif any(term -> occursin("5_prime_utr_variant", term), terms)
        return "5_prime_utr", "utr_variant"
    elseif any(term -> occursin("upstream_gene_variant", term), terms)
        return "upstream", "promoter_candidate"
    elseif any(term -> occursin("transcript_ablation", term), terms)
        return "exon", sv_type == "DEL" ? "transcript_ablation" : "transcript_effect"
    elseif any(term -> occursin("transcript_amplification", term), terms)
        return "exon", "transcript_amplification"
    elseif any(term -> occursin("feature_truncation", term), terms)
        return "exon", "feature_truncation"
    elseif any(term -> occursin("feature_elongation", term), terms)
        return "exon", "feature_elongation"
    elseif any(term -> occursin("intron_variant", term), terms)
        return "intron", "intronic"
    elseif any(term -> occursin("regulatory_region", term), terms)
        return "regulatory", "regulatory_region"
    end
    return "other", ""
end

function apply_manta_sample_fields!(row::Dict{String,Any}, sample_names::AbstractVector{<:AbstractString}, sample_fields::AbstractVector{<:AbstractString}, format_keys::AbstractVector{<:AbstractString})
    for (sample, field_blob) in zip(sample_names, sample_fields)
        value_map = Dict{String,String}()
        for (key, value) in zip(format_keys, split(field_blob, ':', keepempty=true))
            value_map[key] = value
        end
        row["GT ($sample)"] = get(value_map, "GT", "")
        row["FT ($sample)"] = get(value_map, "FT", "")
        row["GQ ($sample)"] = get(value_map, "GQ", "")
        row["PR ($sample)"] = get(value_map, "PR", "")
        row["SR ($sample)"] = get(value_map, "SR", "")
    end
end
