using CodecZlib

const PRIORITY_HEADERS = [
    "candidateCategory",
    "gene",
    "transcript",
    "Consequence",
    "IMPACT",
    "codingEffect",
    "varLocation",
    "chrom",
    "inputPos",
    "inputRef",
    "inputAlt",
    "gNomen",
    "cNomen",
    "pNomen",
    "MANE",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "alleleFreq",
    "GnomAD_v4_1_AF_popmax",
    "GnomAD_v4_1_AF_all",
    "AllofUs250k_gvs_all_af",
    "AllofUs250k_gvs_max_af",
    "spliceai_max",
    "spliceai_summary",
    "Filter (VCF)",
    "Quality (VCF)",
    "MQ",
    "clinVarClinSignifs",
    "hgmdId",
]

function load_inputs(input_path::AbstractString, extra_inputs::Vector{String}; debug::Bool=false)
    rows = Vector{Dict{String,Any}}()
    headers = String[]
    sample_names = String[]
    for path in [input_path; extra_inputs]
        local_rows, local_headers, local_samples = load_input(path; debug=debug)
        append!(rows, local_rows)
        append_missing!(headers, local_headers)
        append_missing!(sample_names, local_samples)
    end
    return rows, headers, sample_names
end

function scan_input_metadata(input_path::AbstractString, extra_inputs::Vector{String}; debug::Bool=false)
    headers = String[]
    sample_names = String[]
    for path in [input_path; extra_inputs]
        if path != input_path && !isfile(path)
            warn_skip_missing_optional(path)
            continue
        end
        local_headers, local_samples = input_metadata(path; debug=debug)
        append_missing!(headers, local_headers)
        append_missing!(sample_names, local_samples)
    end
    return headers, sample_names
end

function input_metadata(path::AbstractString; debug::Bool=false)
    if endswith(lowercase(path), ".vcf") || endswith(lowercase(path), ".vcf.gz")
        if is_manta_vcf(path)
            headers, sample_names, _ = prepare_manta_stream(path)
            return headers, sample_names
        end
        return vcf_metadata(path)
    end
    _, headers, samples = load_tabular(path)
    return headers, samples
end

function load_input(path::AbstractString; debug::Bool=false)
    if endswith(lowercase(path), ".vcf") || endswith(lowercase(path), ".vcf.gz")
        if is_manta_vcf(path)
            headers, sample_names, stream_state = prepare_manta_stream(path)
            rows = Vector{Dict{String,Any}}()
            stream_manta_rows(stream_state) do row
                push!(rows, row)
            end
            return rows, headers, sample_names
        end
        return load_vep_vcf(path; debug=debug)
    end
    return load_tabular(path)
end

function load_tabular(path::AbstractString)
    io = open(path, "r")
    try
        header_line = ""
        while !eof(io)
            header_line = chomp(readline(io))
            isempty(header_line) && continue
            break
        end
        headers = split(startswith(header_line, "#") ? header_line[2:end] : header_line, '\t')
        rows = Vector{Dict{String,Any}}()
        while !eof(io)
            line = chomp(readline(io))
            isempty(line) && continue
            values = split(line, '\t')
            row = Dict{String,Any}()
            for (header, value) in zip(headers, values)
                row[header] = value
            end
            push!(rows, row)
        end
        samples = [header[5:end-1] for header in headers if startswith(header, "GT (")]
        return rows, headers, samples
    finally
        close(io)
    end
end

function load_vep_vcf(path::AbstractString; debug::Bool=false)
    headers, sample_names, stream_state = prepare_vep_stream(path)
    rows = Vector{Dict{String,Any}}()
    stream_vep_rows(stream_state) do row
        push!(rows, row)
    end
    return rows, headers, sample_names
end

function vcf_metadata(path::AbstractString)
    headers, sample_names, _ = prepare_vep_stream(path)
    return headers, sample_names
end

function prepare_vep_stream(path::AbstractString)
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
        assembly = parse_assembly(metadata)
        csq_headers = parse_csq_headers(metadata)
        headers = build_headers(sample_names, csq_headers)
        return headers, sample_names, (io=io, sample_names=sample_names, csq_headers=csq_headers, assembly=assembly, path=String(path))
    catch
        close(io)
        rethrow()
    end
end

function stream_vep_rows(stream_state, callback::Function; debug::Bool=false)
    io = stream_state.io
    try
        while !eof(io)
            line = chomp(readline(io))
            isempty(line) && continue
            startswith(line, "#") && continue
            split_line = split(line, '\t')
            length(split_line) < 8 && continue
            for row in normalise_vcf_record(split_line, stream_state.sample_names, stream_state.csq_headers, stream_state.assembly, stream_state.path; debug=debug)
                callback(row)
            end
        end
    finally
        close(io)
    end
end

function open_text(path::AbstractString)
    if endswith(lowercase(path), ".gz")
        return GzipDecompressorStream(open(path, "r"))
    end
    return open(path, "r")
end

function parse_assembly(metadata::Vector{String})
    for line in metadata
        if startswith(line, "##VEP-command-line=")
            match_obj = match(r"--assembly\s+([A-Za-z0-9._-]+)", line)
            match_obj !== nothing && return match_obj.captures[1]
        end
    end
    return "GRCh38"
end

function parse_csq_headers(metadata::Vector{String})
    for line in metadata
        if occursin("ID=CSQ", line)
            parts = split(line, "Format: "; limit=2)
            if length(parts) == 2
                format_text = strip(parts[2])
                format_text = replace(format_text, "\">" => "")
                format_text = replace(format_text, "\"" => "")
                format_text = replace(format_text, ">" => "")
                return split(strip(format_text), '|')
            end
        end
    end
    preview = join(first(metadata, min(length(metadata), 5)), "\n")
    throw(ArgumentError("CSQ header not found in VCF metadata. First metadata lines:\n" * preview))
end

function build_headers(sample_names::Vector{<:AbstractString}, csq_headers::Vector{<:AbstractString})
    headers = copy(PRIORITY_HEADERS)
    append!(headers, ["geneId", "varType", "distNearestSS", "protein", "assembly", "nearestSSChange", "localSpliceEffect", "hgmdSubCategory", "hgmdPhenotype", "clinVarReviewStatus", "clinVarPhenotypes", "gDNAstart", "gDNAend", "wtCodon", "varCodon", "branchPointChange", "uorfEffect", "shortenedIntron", "Consequence"])
    for sample in sample_names
        append!(headers, ["GT ($sample)", "AD ($sample)", "DP ($sample)", "GQ ($sample)", "PL ($sample)"])
    end
    return headers
end

function normalise_vcf_record(fields::AbstractVector{<:AbstractString}, sample_names::AbstractVector{<:AbstractString}, csq_headers::AbstractVector{<:AbstractString}, assembly::AbstractString, filename::AbstractString; debug::Bool=false)
    chrom, pos, _, ref, alt, qual, filter_field, info_field = fields[1:8]
    format_keys = length(fields) >= 9 ? split(fields[9], ':') : String[]
    sample_fields = length(fields) >= 10 ? fields[10:end] : String[]
    info_map = parse_info_field(info_field)
    haskey(info_map, "CSQ") || return Dict{String,Any}[]
    alts = split(alt, ',')
    csq_entries = split(String(info_map["CSQ"]), ',')
    rows = Vector{Dict{String,Any}}()
    for csq_entry in csq_entries
        csq_values = split(csq_entry, '|', keepempty=true)
        length(csq_values) == length(csq_headers) || continue
        csq_map = Dict{String,String}(zip(csq_headers, csq_values))
        get(csq_map, "Feature_type", "") == "Transcript" || continue
        isempty(get(csq_map, "SYMBOL", "")) && continue
        alt_number = match_vep_allele(ref, alts, get(csq_map, "Allele", ""))
        alt_number == 0 && continue
        input_alt = alts[alt_number]
        row = Dict{String,Any}()
        consequence = get(csq_map, "Consequence", "")
        var_location, coding_effect, var_type = classify_variant(consequence, ref, input_alt)
        row["candidateCategory"] = ""
        row["chrom"] = chrom
        row["inputPos"] = pos
        row["inputRef"] = ref
        row["inputAlt"] = input_alt
        row["varType"] = var_type
        row["varLocation"] = var_location
        row["codingEffect"] = coding_effect
        row["IMPACT"] = get(csq_map, "IMPACT", "")
        row["nearestSSChange"] = ""
        row["localSpliceEffect"] = ""
        row["hgmdSubCategory"] = get(csq_map, "HGMD_CLASS", "")
        row["hgmdPhenotype"] = get(csq_map, "HGMD_PHEN", "")
        row["hgmdId"] = get(csq_map, "HGMD", "")
        row["clinVarClinSignifs"] = get(csq_map, "ClinVar_CLNSIG", "")
        row["clinVarReviewStatus"] = get(csq_map, "ClinVar_CLNREVSTAT", "")
        row["clinVarPhenotypes"] = get(csq_map, "ClinVar_CLNDN", "")
        row["geneId"] = replace(get(csq_map, "HGNC_ID", ""), "HGNC:" => "")
        row["gene"] = get(csq_map, "SYMBOL", "")
        row["distNearestSS"] = normalise_splice_distance(get(csq_map, "SpliceDistance", ""))
        row["protein"] = get(csq_map, "ENSP", "")
        row["transcript"] = get(csq_map, "Feature_type", "") == "Transcript" ? get(csq_map, "Feature", "") : ""
        row["MANE"] = get(csq_map, "MANE", "")
        row["MANE_SELECT"] = get(csq_map, "MANE_SELECT", "")
        row["MANE_PLUS_CLINICAL"] = get(csq_map, "MANE_PLUS_CLINICAL", "")
        row["GnomAD_v4_1_AF_all"] = first_nonempty(csq_map, ["GnomAD_v4_1_AF_all", "gnomADg_AF", "gnomadAltFreq_all"])
        row["GnomAD_v4_1_AF_popmax"] = first_nonempty(csq_map, ["GnomAD_v4_1_AF_popmax"])
        row["AllofUs250k_gvs_all_af"] = get(csq_map, "AllofUs250k_gvs_all_af", "")
        row["AllofUs250k_gvs_max_af"] = first_nonempty([
            alt_specific_value(info_map, "AllofUs250k_gvs_max_af", alt_number),
            get(csq_map, "AllofUs250k_gvs_max_af", ""),
        ])
        row["GnomAD_v4_1_AC_all"] = get(csq_map, "GnomAD_v4_1_AC_all", "")
        row["GnomAD_v4_1_N_Hom_all"] = get(csq_map, "GnomAD_v4_1_N_Hom_all", "")
        spliceai_values = spliceai_scores(csq_map)
        row["spliceai_max"] = string(maximum(spliceai_values))
        row["spliceai_summary"] = spliceai_summary(spliceai_values)
        row["gDNAstart"] = pos
        row["gDNAend"] = string(parse(Int, pos) + max(length(ref), length(input_alt)) - 1)
        wt_codon, var_codon = parse_codons(get(csq_map, "Codons", ""))
        row["wtCodon"] = wt_codon
        row["varCodon"] = var_codon
        row["gNomen"] = string(chrom, ":", pos, ref, ">", input_alt)
        row["cNomen"] = suffix_after_colon(get(csq_map, "HGVSc", ""))
        row["pNomen"] = suffix_after_colon(get(csq_map, "HGVSp", ""))
        row["branchPointChange"] = get(csq_map, "BranchPointDistance", "")
        row["uorfEffect"] = get(csq_map, "5UTR_consequence", "")
        row["shortenedIntron"] = get(csq_map, "shortenedIntron", "")
        row["assembly"] = assembly
        row["gnomadFilter"] = get(info_map, "gnomAD_FILTER", "")
        row["Quality (VCF)"] = qual
        row["Filter (VCF)"] = filter_field
        row["MQ"] = get(info_map, "MQ", "")
        row["alleleFreq"] = alt_specific_value(info_map, "AF", alt_number)
        row["Consequence"] = consequence
        apply_sample_fields!(row, sample_names, sample_fields, format_keys, alt_number)
        push!(rows, row)
    end
    return rows
end

function parse_info_field(text::AbstractString)
    info = Dict{String,Any}()
    for part in split(text, ';')
        isempty(part) && continue
        if occursin('=', part)
            key, value = split(part, '='; limit=2)
            info[key] = value
        else
            info[part] = true
        end
    end
    return info
end

function apply_sample_fields!(row::Dict{String,Any}, sample_names::AbstractVector{<:AbstractString}, sample_fields::AbstractVector{<:AbstractString}, format_keys::AbstractVector{<:AbstractString}, alt_number::Int)
    for (sample, field_blob) in zip(sample_names, sample_fields)
        value_map = Dict{String,String}()
        for (key, value) in zip(format_keys, split(field_blob, ':', keepempty=true))
            value_map[key] = value
        end
        row["GT ($sample)"] = get(value_map, "GT", "")
        row["AD ($sample)"] = alt_specific_ad(get(value_map, "AD", ""), alt_number)
        row["DP ($sample)"] = get(value_map, "DP", "")
        row["GQ ($sample)"] = get(value_map, "GQ", "")
        row["PL ($sample)"] = get(value_map, "PL", "")
    end
end

function match_vep_allele(ref::AbstractString, alts::AbstractVector{<:AbstractString}, vep_allele::AbstractString)
    allele = vep_allele == "-" ? "" : vep_allele
    missing_bases = 1
    for alt in alts
        alt == "*" && continue
        if !isempty(ref) && !isempty(alt) && first(ref) != first(alt)
            missing_bases = 0
        end
    end
    for (index, alt) in enumerate(alts)
        alt == "*" && continue
        if missing_bases < length(alt) && allele == alt[(missing_bases + 1):end]
            return index
        end
        if allele == alt
            return index
        end
        if isempty(allele) && length(ref) > length(alt)
            return index
        end
        if length(alt) > length(ref)
            inserted = alt[(length(ref) + 1):end]
            allele == inserted && return index
        end
    end
    return 0
end

function classify_variant(consequence::AbstractString, ref::AbstractString, alt::AbstractString)
    var_location = if occursin("splice_", consequence)
        "splice site"
    elseif any(occursin(token, consequence) for token in [
        "coding_sequence_variant",
        "protein_altering_variant",
        "frameshift_variant",
        "inframe_deletion",
        "inframe_insertion",
        "missense_variant",
        "non_coding_transcript_exon_variant",
        "start_lost",
        "stop_gained",
        "stop_lost",
        "stop_retained_variant",
        "synonymous_variant",
    ])
        "exon"
    elseif occursin("3_prime_UTR_variant", consequence)
        "3'UTR"
    elseif occursin("5_prime_UTR_variant", consequence)
        "5'UTR"
    elseif occursin("intron_variant", consequence)
        "intron"
    elseif occursin("upstream_gene_variant", consequence)
        "upstream"
    elseif occursin("downstream_gene_variant", consequence)
        "downstream"
    else
        consequence
    end
    coding_effect = if occursin("start_lost", consequence)
        "start loss"
    elseif occursin("stop_gained", consequence)
        "stop gain"
    elseif occursin("stop_lost", consequence)
        "stop loss"
    elseif occursin("frameshift_variant", consequence)
        "frameshift"
    elseif occursin("inframe_deletion", consequence) || occursin("inframe_insertion", consequence)
        "in-frame"
    elseif occursin("missense_variant", consequence)
        "missense"
    elseif occursin("protein_altering_variant", consequence) || occursin("coding_sequence_variant", consequence)
        "protein-altering"
    elseif occursin("synonymous_variant", consequence) || occursin("stop_retained_variant", consequence)
        "synonymous"
    else
        ""
    end
    var_type = if length(ref) == 1 && length(alt) == 1
        "substitution"
    elseif length(ref) > length(alt)
        "deletion"
    elseif length(ref) < length(alt)
        "insertion"
    else
        consequence
    end
    return var_location, coding_effect, var_type
end

function normalise_splice_distance(text::AbstractString)
    isempty(text) && return ""
    value = split(text, '&')[1]
    if startswith(value, "acceptor")
        return value[9:end]
    elseif startswith(value, "donor")
        return value[6:end]
    end
    return value
end

function parse_codons(text::AbstractString)
    length(text) == 7 || return "", ""
    return text[1:3], text[5:7]
end

function suffix_after_colon(text::AbstractString)
    index = findfirst(':', text)
    index === nothing && return ""
    return text[(index + 1):end]
end

function alt_specific_value(info::Dict{String,Any}, key::AbstractString, alt_number::Int)
    raw = get(info, key, "")
    raw isa Bool && return ""
    values = split(String(raw), ',')
    if alt_number <= length(values)
        return values[alt_number]
    end
    return values[1]
end

function alt_specific_ad(ad_text::AbstractString, alt_number::Int)
    isempty(ad_text) && return ""
    values = split(ad_text, ',')
    length(values) < 2 && return ad_text
    alt_index = min(alt_number + 1, length(values))
    return string(values[1], ",", values[alt_index])
end

function first_nonempty(values::Dict{String,String}, keys::AbstractVector{<:AbstractString})
    for key in keys
        value = get(values, key, "")
        !isempty(value) && return value
    end
    return ""
end

function first_nonempty(values::AbstractVector{<:AbstractString})
    for value in values
        !isempty(value) && return String(value)
    end
    return ""
end

function spliceai_scores(csq_map::Dict{String,String})
    return [
        parse_float_string(get(csq_map, "SpliceAI_pred_DS_AG", "")),
        parse_float_string(get(csq_map, "SpliceAI_pred_DS_AL", "")),
        parse_float_string(get(csq_map, "SpliceAI_pred_DS_DG", "")),
        parse_float_string(get(csq_map, "SpliceAI_pred_DS_DL", "")),
    ]
end

function spliceai_summary(scores::Vector{Float64})
    labels = ["AG", "AL", "DG", "DL"]
    parts = String[]
    for (label, score) in zip(labels, scores)
        score > 0 || continue
        push!(parts, string(label, "=", round(score, digits=3)))
    end
    return join(parts, ", ")
end

function parse_float_string(text::AbstractString)
    isempty(text) && return 0.0
    parsed = tryparse(Float64, text)
    return parsed === nothing ? 0.0 : parsed
end

function append_missing!(target::Vector{String}, values::AbstractVector{<:AbstractString})
    for value in values
        string_value = String(value)
        string_value in target || push!(target, string_value)
    end
    return target
end
