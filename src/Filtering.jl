function prioritise_rows(rows::Vector{Dict{String,Any}}, family::FamilySpec, options::RunOptions, config::AppConfig)
    thresholds = config.thresholds
    freq_cutoff = something(options.freq_cutoff, thresholds.frequency_cutoff)
    gq_cutoff = something(options.gq_cutoff, thresholds.gq_cutoff)
    gq_hom_cutoff = something(options.gq_hom_cutoff, thresholds.gq_hom_cutoff)
    mq_cutoff = something(options.mq_cutoff, thresholds.mq_cutoff)
    root_dir = batch_root(options.input, options.pipeline_prefix)
    denovocnn_calls = load_denovocnn_calls(options.pipeline_prefix, root_dir, denovocnn_proband(family))
    sample_sexes = load_sample_sexes(family, options.pipeline_prefix, root_dir)

    filtered = Vector{Dict{String,Any}}()
    for row in rows
        prioritise_row!(filtered, row, family, options, thresholds, freq_cutoff, gq_cutoff, gq_hom_cutoff, mq_cutoff, denovocnn_calls, sample_sexes)
    end
    return postprocess_prioritised_rows(filtered, family)
end

function prioritise_row!(filtered::Vector{Dict{String,Any}}, row::Dict{String,Any}, family::FamilySpec, options::RunOptions, thresholds::ThresholdConfig, freq_cutoff::Float64, gq_cutoff::Float64, gq_hom_cutoff::Float64, mq_cutoff::Float64, denovocnn_calls::Set{String}, sample_sexes::Dict{String,String})
    is_denovocnn_hit = row_variant_key(row) in denovocnn_calls
    row["__is_denovocnn_hit__"] = is_denovocnn_hit
    passes_gene_filters(row, options, thresholds) || return
    category = classify_family_model(row, family, gq_cutoff, gq_hom_cutoff, denovocnn_calls, sample_sexes)
    isempty(category) && return
    if !is_denovocnn_hit || category == "recessive_homozygous_candidate"
        passes_quality_filter(row, options.lowq, mq_cutoff) || return
        is_interesting_variant(row, options.splice, thresholds) || return
    end
    passes_frequency_filter(row, category, freq_cutoff, thresholds) || return
    if options.recessive_only && !occursin("recessive", category)
        return
    end
    if options.denovo_only && !occursin("de_novo", category)
        return
    end
    row["candidateCategory"] = category
    push!(filtered, row)
end

function prioritise_manta_row!(filtered::Vector{Dict{String,Any}}, row::Dict{String,Any}, family::FamilySpec, options::RunOptions, thresholds::ThresholdConfig)
    row["analysis_mode"] = (!isnothing(singleton_sample(family)) && isnothing(family.parent1) && isnothing(family.parent2)) ? "singleton" : "family"
    passes_manta_gene_filters(row, thresholds) || return
    passes_manta_quality_filter(row, family, thresholds) || return
    passes_manta_consequence_filter(row, thresholds) || return
    category = classify_manta_category(row)
    isempty(category) && return
    row["candidateCategory"] = category
    row["inheritance_model"] = classify_manta_inheritance(row, family)
    push!(filtered, row)
end

function prioritisation_context(options::RunOptions, config::AppConfig)
    thresholds = config.thresholds
    root = batch_root(options.input, options.pipeline_prefix)
    return (
        thresholds=thresholds,
        freq_cutoff=something(options.freq_cutoff, thresholds.frequency_cutoff),
        gq_cutoff=something(options.gq_cutoff, thresholds.gq_cutoff),
        gq_hom_cutoff=something(options.gq_hom_cutoff, thresholds.gq_hom_cutoff),
        mq_cutoff=something(options.mq_cutoff, thresholds.mq_cutoff),
        denovocnn_calls=load_denovocnn_calls(options.pipeline_prefix, root, denovocnn_proband(resolve_family(options))),
        sample_sexes=load_sample_sexes(resolve_family(options), options.pipeline_prefix, root),
    )
end

function passes_gene_filters(row::Dict{String,Any}, options::RunOptions, thresholds::ThresholdConfig)
    gene = string(get(row, "gene", ""))
    if !isempty(options.only_genes) && !(gene in options.only_genes)
        return false
    end
    if startswith(gene, "MUC") || startswith(gene, "HLA")
        return false
    end
    var_location = string(get(row, "varLocation", ""))
    distance = parse_float(get(row, "distNearestSS", ""))
    if (var_location == "upstream" || var_location == "downstream") && !isnan(distance) && abs(distance) > thresholds.gene_distance_cutoff
        return false
    end
    return true
end

function passes_manta_gene_filters(row::Dict{String,Any}, thresholds::ThresholdConfig)
    gene = string(get(row, "gene", ""))
    isempty(gene) && return false
    startswith(gene, "MUC") && return false
    startswith(gene, "HLA") && return false
    return true
end

function passes_quality_filter(row::Dict{String,Any}, lowq::Bool, mq_cutoff::Float64)
    lowq && return true
    filter_field = uppercase(string(get(row, "Filter (VCF)", "")))
    var_type = string(get(row, "varType", ""))
    allowed = if var_type == "substitution"
        Set(["PASS", ".", "MQ40"])
    else
        Set(["PASS", ".", "QD2", "MQ40", "MQ40;QD2"])
    end
    mq = parse_float(get(row, "MQ", ""))
    if filter_field in allowed
        return true
    end
    return !isnan(mq) && mq > mq_cutoff
end

function passes_manta_quality_filter(row::Dict{String,Any}, family::FamilySpec, thresholds::ThresholdConfig)
    uppercase(string(get(row, "Filter (VCF)", ""))) == "PASS" || return false
    isempty(family.affected) && return false
    proband = family.affected[1]
    child_state = genotype_state(row, proband)
    child_state in (:het, :hom_alt) || return false
    uppercase(string(get(row, "FT ($proband)", ""))) == "PASS" || return false
    parse_float(get(row, "GQ ($proband)", "")) >= thresholds.manta_gq_cutoff || return false
    sv_size = parse_sv_size(row)
    sv_size >= thresholds.manta_min_size || return false
    pr_alt = sample_alt_support(row, proband, "PR")
    sr_alt = sample_alt_support(row, proband, "SR")
    sr_alt >= thresholds.manta_min_sr_alt_support || return false
    (pr_alt + sr_alt) >= thresholds.manta_min_total_alt_support || return false
    for parent in filter(!isnothing, [family.parent1, family.parent2])
        passes_manta_parent_qc(row, parent, thresholds) || return false
    end
    if !isnothing(family.parent1) && !isnothing(family.parent2)
        p1 = genotype_state(row, family.parent1)
        p2 = genotype_state(row, family.parent2)
        p1 == :hom_alt && return false
        p2 == :hom_alt && return false
        if child_state == :het
            p1 == :ref || return false
            p2 == :ref || return false
        elseif child_state == :hom_alt
            p1 == :het || return false
            p2 == :het || return false
        end
    end
    return true
end

function passes_manta_parent_qc(row::Dict{String,Any}, parent::String, thresholds::ThresholdConfig)
    state = genotype_state(row, parent)
    state == :missing && return false
    uppercase(string(get(row, "FT ($parent)", ""))) == "PASS" || return false
    parse_float(get(row, "GQ ($parent)", "")) >= thresholds.manta_gq_cutoff || return false
    if state in (:het, :hom_alt)
        total_alt = sample_alt_support(row, parent, "PR") + sample_alt_support(row, parent, "SR")
        total_alt >= thresholds.manta_min_parent_alt_support || return false
    end
    return true
end

function is_interesting_variant(row::Dict{String,Any}, include_extended_splice::Bool, thresholds::ThresholdConfig)
    consequence = lowercase(string(get(row, "Consequence", "")))
    impact = uppercase(string(get(row, "IMPACT", "")))
    var_location = string(get(row, "varLocation", ""))
    coding_effect = string(get(row, "codingEffect", ""))
    if impact in ("HIGH", "MODERATE")
        return true
    elseif impact in ("LOW", "MODIFIER")
        return has_spliceai_support(row, thresholds)
    elseif !isempty(impact)
        return false
    end
    occursin("splice_", consequence) && return true
    if is_protein_altering(coding_effect)
        return true
    end
    if coding_effect == "synonymous" || var_location == "intron"
        return has_spliceai_support(row, thresholds)
    end
    if var_location == "exon" && isempty(coding_effect)
        return false
    end
    return false
end

function passes_manta_consequence_filter(row::Dict{String,Any}, thresholds::ThresholdConfig)
    sv_type = uppercase(string(get(row, "svType", "")))
    consequence = lowercase(string(get(row, "Consequence", "")))
    impact = uppercase(string(get(row, "IMPACT", "")))
    protein_coding = manta_is_protein_coding(row)
    impact_hit = impact in ("HIGH", "MODERATE")
    coding_hit = protein_coding && (impact_hit || manta_coding_consequence(consequence))
    is_singleton_mode = string(get(row, "analysis_mode", "")) == "singleton"
    if is_singleton_mode
        if sv_type in ("DEL", "INS")
            return protein_coding && (impact_hit || (isempty(impact) && manta_coding_consequence(consequence)))
        end
        return false
    end
    if sv_type in ("DUP", "INS")
        return protein_coding && (impact_hit || (isempty(impact) && manta_coding_consequence(consequence)))
    elseif sv_type == "DEL"
        coding_hit && return true
        protein_coding && occursin("5_prime_utr_variant", consequence) && return true
        if occursin("upstream_gene_variant", consequence)
            distance = parse_float(get(row, "DISTANCE", get(row, "distNearestSS", "")))
            return protein_coding && !isnan(distance) && abs(distance) <= thresholds.gene_distance_cutoff
        end
        return false
    end
    return false
end

function passes_frequency_filter(row::Dict{String,Any}, category::String, freq_cutoff::Float64, thresholds::ThresholdConfig)
    cutoff = if category == "de_novo_candidate"
        thresholds.denovo_frequency_cutoff
    elseif category == "recessive_homozygous_candidate"
        thresholds.homozygous_frequency_cutoff
    else
        freq_cutoff
    end
    return row_max_frequency(row) <= cutoff
end

function row_max_frequency(row::Dict{String,Any})
    values = [
        parse_float(get(row, "GnomAD_v4_1_AF_popmax", "")),
        parse_float(get(row, "GnomAD_v4_1_AF_all", "")),
        parse_float(get(row, "gnomadAltFreq_all", "")),
        parse_float(get(row, "AllofUs250k_gvs_max_af", "")),
    ]
    filtered = [value for value in values if !isnan(value)]
    isempty(filtered) && return 0.0
    return maximum(filtered)
end

function parse_sv_size(row::Dict{String,Any})
    raw = parse_float(get(row, "SVLEN", ""))
    isnan(raw) && return 0
    return Int(round(abs(raw)))
end

function sample_alt_support(row::Dict{String,Any}, sample::String, field::String)
    raw = string(get(row, "$field ($sample)", ""))
    isempty(raw) && return 0
    raw == "." && return 0
    parts = split(raw, ',')
    length(parts) < 2 && return 0
    parsed = tryparse(Int, parts[2])
    return parsed === nothing ? 0 : parsed
end

function manta_coding_consequence(consequence::AbstractString)
    terms = split(lowercase(consequence), '&')
    coding_terms = (
        "coding_sequence_variant",
        "transcript_ablation",
        "transcript_amplification",
        "feature_truncation",
        "feature_elongation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "splice_donor_5th_base_variant",
        "stop_gained",
        "stop_lost",
        "start_lost",
        "frameshift_variant",
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant",
    )
    return any(term -> term in coding_terms, terms)
end

function manta_is_protein_coding(row::Dict{String,Any})
    return lowercase(string(get(row, "biotype", ""))) == "protein_coding"
end

function classify_manta_category(row::Dict{String,Any})
    sv_type = uppercase(string(get(row, "svType", "")))
    if string(get(row, "analysis_mode", "")) == "singleton"
        sv_type == "DEL" && return "manta_deletion_candidate"
        sv_type == "INS" && return "manta_insertion_candidate"
        return ""
    end
    sv_type == "DEL" && return "manta_deletion_candidate"
    sv_type == "DUP" && return "manta_duplication_candidate"
    sv_type == "INS" && return "manta_insertion_candidate"
    return ""
end

function classify_manta_inheritance(row::Dict{String,Any}, family::FamilySpec)
    isempty(family.affected) && return ""
    child = family.affected[1]
    child_state = genotype_state(row, child)
    child_state in (:het, :hom_alt) || return ""
    if isnothing(family.parent1) || isnothing(family.parent2)
        return "present_in_proband"
    end
    p1 = genotype_state(row, family.parent1)
    p2 = genotype_state(row, family.parent2)
    p1_ft = uppercase(string(get(row, "FT ($(family.parent1))", "")))
    p2_ft = uppercase(string(get(row, "FT ($(family.parent2))", "")))
    p1_pass = p1_ft == "PASS"
    p2_pass = p2_ft == "PASS"
    if p1 == :ref && p2 == :ref && p1_pass && p2_pass
        return "apparent_de_novo"
    elseif child_state == :het
        return "filtered_non_denovo_het"
    elseif (p1 in (:het, :hom_alt)) && p2 == :ref && p1_pass
        return "maternal"
    elseif (p2 in (:het, :hom_alt)) && p1 == :ref && p2_pass
        return "paternal"
    elseif (p1 in (:het, :hom_alt)) && (p2 in (:het, :hom_alt)) && p1_pass && p2_pass
        return "biparental"
    end
    return "uncertain"
end

function classify_family_model(row::Dict{String,Any}, family::FamilySpec, gq_cutoff::Float64, gq_hom_cutoff::Float64, denovocnn_calls::Set{String}, sample_sexes::Dict{String,String})
    family.shared && return classify_shared_variant(row, family)
    x_linked = classify_x_linked_variant(row, family, sample_sexes, gq_hom_cutoff)
    isempty(x_linked) || return x_linked
    if !isnothing(family.parent1) && !isnothing(family.parent2) && !isempty(family.affected)
        length(family.affected) == 1 && return classify_trio_variant(row, family, gq_cutoff, gq_hom_cutoff, denovocnn_calls)
        return classify_family_variant(row, family)
    elseif known_parent_count(family) == 1 && length(family.affected) == 1
        return classify_one_parent_variant(row, family, gq_cutoff, gq_hom_cutoff)
    elseif !isempty(family.affected) && (length(family.affected) > 1 || !isempty(family.unaffected))
        return classify_segregating_variant(row, family)
    elseif !isnothing(singleton_sample(family))
        return classify_singleton_variant(row, singleton_sample(family), gq_cutoff, gq_hom_cutoff)
    elseif !isempty(family.affected)
        return classify_shared_variant(row, family)
    elseif !isnothing(family.parent1) && !isnothing(family.parent2)
        p1 = genotype_state(row, family.parent1)
        p2 = genotype_state(row, family.parent2)
        return (p1 != :ref || p2 != :ref) ? "parental_variant" : ""
    end
    return ""
end

function known_parent_count(family::FamilySpec)
    return (!isnothing(family.parent1) ? 1 : 0) + (!isnothing(family.parent2) ? 1 : 0)
end

function known_parent(family::FamilySpec)
    return isnothing(family.parent1) ? family.parent2 : family.parent1
end

function singleton_sample(family::FamilySpec)
    return length(family.affected) == 1 ? family.affected[1] : nothing
end

function classify_singleton_variant(row::Dict{String,Any}, sample::String, gq_cutoff::Float64, gq_hom_cutoff::Float64)
    autosomal_or_x(row) || return ""
    state = genotype_state(row, sample)
    gq = parse_float(get(row, "GQ ($sample)", ""))
    if state == :hom_alt && gq >= gq_hom_cutoff
        return "recessive_homozygous_candidate"
    elseif state == :het && gq >= gq_cutoff && singleton_comphet_eligible(row)
        return "singleton_possible_compound_heterozygous_component"
    end
    return ""
end

function classify_one_parent_variant(row::Dict{String,Any}, family::FamilySpec, gq_cutoff::Float64, gq_hom_cutoff::Float64)
    autosomal_or_x(row) || return ""
    child = family.affected[1]
    parent = known_parent(family)
    parent === nothing && return ""
    child_state = genotype_state(row, child)
    parent_state = genotype_state(row, parent)
    child_gq = parse_float(get(row, "GQ ($child)", ""))
    if child_state == :hom_alt && child_gq >= gq_hom_cutoff && parent_state in (:het, :hom_alt)
        return "recessive_homozygous_candidate"
    elseif child_state == :het && child_gq >= gq_cutoff
        if parent_state in (:het, :hom_alt)
            return "one_parent_inherited_candidate"
        elseif parent_state == :ref
            return "one_parent_other_parent_candidate"
        end
    end
    return ""
end

function classify_segregating_variant(row::Dict{String,Any}, family::FamilySpec)
    autosomal_or_x(row) || return ""
    affected_states = [genotype_state(row, sample) for sample in family.affected]
    all(state != :ref && state != :missing for state in affected_states) || return ""
    unaffected_states = [genotype_state(row, sample) for sample in family.unaffected]
    if all(state == :hom_alt for state in affected_states) && all(state != :hom_alt for state in unaffected_states)
        return "recessive_family_candidate"
    end
    all(state == :ref || state == :missing for state in unaffected_states) || return ""
    return "segregating_family_candidate"
end

function singleton_comphet_eligible(row::Dict{String,Any})
    consequence = lowercase(string(get(row, "Consequence", "")))
    impact = uppercase(string(get(row, "IMPACT", "")))
    coding_effect = string(get(row, "codingEffect", ""))
    if impact in ("HIGH", "MODERATE")
        return true
    elseif impact in ("LOW", "MODIFIER")
        return false
    elseif !isempty(impact)
        return false
    end
    if is_protein_altering(coding_effect)
        return true
    end
    return occursin("splice_acceptor_variant", consequence) || occursin("splice_donor_variant", consequence)
end

function classify_trio_variant(row::Dict{String,Any}, family::FamilySpec, gq_cutoff::Float64, gq_hom_cutoff::Float64, denovocnn_calls::Set{String})
    autosomal_or_x(row) || return ""
    child = family.affected[1]
    p1 = genotype_state(row, family.parent1)
    p2 = genotype_state(row, family.parent2)
    child_state = genotype_state(row, child)
    child_gq = parse_float(get(row, "GQ ($child)", ""))

    if p1 == :ref && p2 == :ref && (child_state == :het || child_state == :hom_alt) && row_variant_key(row) in denovocnn_calls
        return "de_novo_candidate"
    end
    if p1 == :het && p2 == :het && child_state == :hom_alt && child_gq >= gq_hom_cutoff
        return "recessive_homozygous_candidate"
    end
    if child_state == :het && child_gq >= gq_cutoff
        if p1 == :het || p1 == :hom_alt || p2 == :het || p2 == :hom_alt
            return "inherited_candidate"
        end
    end
    return ""
end

function classify_family_variant(row::Dict{String,Any}, family::FamilySpec)
    autosomal_or_x(row) || return ""
    proband = family.affected[1]
    proband_state = genotype_state(row, proband)
    p1 = genotype_state(row, family.parent1)
    p2 = genotype_state(row, family.parent2)
    if p1 == :ref && p2 == :ref && (proband_state == :het || proband_state == :hom_alt) && get(row, "__is_denovocnn_hit__", false) == true
        return "de_novo_candidate"
    end
    affected_states = [genotype_state(row, sample) for sample in family.affected]
    affected_gqs = [parse_float(get(row, "GQ ($sample)", "")) for sample in family.affected]
    if all(state == :hom_alt for state in affected_states) && all(gq >= 5.0 for gq in affected_gqs)
        if p1 == :het && p2 == :het
            return "recessive_homozygous_candidate"
        end
        return ""
    end
    if all(state == :het for state in affected_states)
        origin = inherited_parent(row, family)
        if origin == :maternal
            return "family_inherited_maternal_candidate"
        elseif origin == :paternal
            return "family_inherited_paternal_candidate"
        end
    end
    return ""
end

function classify_shared_variant(row::Dict{String,Any}, family::FamilySpec)
    for sample in family.affected
        genotype_state(row, sample) == :ref && return ""
    end
    for sample in family.unaffected
        genotype_state(row, sample) == :hom_alt && return ""
    end
    return "shared_candidate"
end

function genotype_state(row::Dict{String,Any}, sample::String)
    gt = replace(string(get(row, "GT ($sample)", "")), '|' => '/')
    gt in ("", ".", "./.", ".//.") && return :missing
    gt == "0" && return :ref
    gt == "1" && return :hom_alt
    gt == "0/0" && return :ref
    gt in ("0/1", "1/0", "0/2", "2/0", "1/2", "2/1") && return :het
    if occursin(r"^([1-9][0-9]*)/\1$", gt)
        return :hom_alt
    end
    return :other
end

function autosomal_or_x(row::Dict{String,Any})
    chrom = uppercase(normalise_chromosome(string(get(row, "chrom", ""))))
    return !(chrom in ("M", "MT"))
end

function variant_sort_key(row::Dict{String,Any})
    category_rank = Dict(
        "de_novo_candidate" => 1,
        "recessive_homozygous_candidate" => 2,
        "x_linked_recessive_candidate" => 3,
        "compound_heterozygous_candidate" => 4,
        "imprinted_gene_flag" => 5,
    )
    category = string(get(row, "candidateCategory", ""))
    rank = get(category_rank, category, 99)
    frequency = row_max_frequency(row)
    gene = string(get(row, "gene", ""))
    position = something(tryparse(Int, string(get(row, "inputPos", "0"))), 0)
    return (rank, frequency, gene, position)
end

function is_protein_altering(coding_effect::AbstractString)
    return coding_effect in ("start loss", "stop gain", "stop loss", "frameshift", "in-frame", "missense", "protein-altering")
end

function has_spliceai_support(row::Dict{String,Any}, thresholds::ThresholdConfig)
    spliceai_max = parse_float(get(row, "spliceai_max", ""))
    return !isnan(spliceai_max) && spliceai_max >= thresholds.spliceai_cutoff
end

function postprocess_prioritised_rows(rows::Vector{Dict{String,Any}}, family::FamilySpec)
    collapsed = collapse_transcript_rows(rows)
    reclassified = Vector{Dict{String,Any}}()
    if !isnothing(singleton_sample(family)) && isnothing(family.parent1) && isnothing(family.parent2)
        singleton_comphet_genes = identify_singleton_possible_compound_heterozygous_genes(collapsed, family)
        for row in collapsed
            category = string(get(row, "candidateCategory", ""))
            if category == "singleton_possible_compound_heterozygous_component"
                string(get(row, "gene", "")) in singleton_comphet_genes || continue
                row["candidateCategory"] = "singleton_possible_compound_heterozygous_candidate"
            elseif !(category in ("recessive_homozygous_candidate", "x_linked_recessive_candidate", "manta_deletion_candidate", "manta_insertion_candidate"))
                continue
            end
            push!(reclassified, row)
        end
        sort!(reclassified, by=variant_sort_key)
        return reclassified
    end
    if known_parent_count(family) == 1 && length(family.affected) == 1
        comphet_keys = identify_one_parent_compound_heterozygous_variants(collapsed, family)
        for row in collapsed
            category = string(get(row, "candidateCategory", ""))
            key = row_variant_key(row)
            if category in ("one_parent_inherited_candidate", "one_parent_other_parent_candidate")
                key in comphet_keys || continue
                row["candidateCategory"] = "compound_heterozygous_candidate"
            elseif !(category in ("recessive_homozygous_candidate", "x_linked_recessive_candidate", "manta_deletion_candidate", "manta_duplication_candidate", "manta_insertion_candidate"))
                continue
            end
            push!(reclassified, row)
        end
        sort!(reclassified, by=variant_sort_key)
        return reclassified
    end
    if !isempty(family.affected) && (length(family.affected) > 1 || !isempty(family.unaffected)) && known_parent_count(family) < 2
        for row in collapsed
            category = string(get(row, "candidateCategory", ""))
            if !(category in ("recessive_family_candidate", "segregating_family_candidate", "x_linked_recessive_candidate", "manta_deletion_candidate", "manta_duplication_candidate", "manta_insertion_candidate"))
                continue
            end
            push!(reclassified, row)
        end
        sort!(reclassified, by=variant_sort_key)
        return reclassified
    end
    if known_parent_count(family) == 2 && length(family.affected) > 1
        comphet_keys = identify_family_compound_heterozygous_variants(collapsed, family)
        for row in collapsed
            category = string(get(row, "candidateCategory", ""))
            key = row_variant_key(row)
            if category in ("family_inherited_maternal_candidate", "family_inherited_paternal_candidate")
                key in comphet_keys || continue
                row["candidateCategory"] = "compound_heterozygous_candidate"
            elseif !(category in ("de_novo_candidate", "recessive_homozygous_candidate", "x_linked_recessive_candidate", "manta_deletion_candidate", "manta_duplication_candidate", "manta_insertion_candidate"))
                continue
            end
            push!(reclassified, row)
        end
        sort!(reclassified, by=variant_sort_key)
        return reclassified
    end
    comphet_keys = identify_compound_heterozygous_variants(collapsed, family)
    for row in collapsed
        category = string(get(row, "candidateCategory", ""))
        key = row_variant_key(row)
        if category == "inherited_candidate"
            key in comphet_keys || continue
            row["candidateCategory"] = "compound_heterozygous_candidate"
        elseif !(category in ("de_novo_candidate", "recessive_homozygous_candidate", "x_linked_recessive_candidate", "manta_deletion_candidate", "manta_duplication_candidate", "manta_insertion_candidate"))
            continue
        end
        push!(reclassified, row)
    end
    sort!(reclassified, by=variant_sort_key)
    return reclassified
end

function identify_family_compound_heterozygous_variants(rows::Vector{Dict{String,Any}}, family::FamilySpec)
    (known_parent_count(family) == 2 && length(family.affected) > 1) || return Set{String}()
    grouped = Dict{String,Dict{Symbol,Set{String}}}()
    for row in rows
        category = string(get(row, "candidateCategory", ""))
        category in ("family_inherited_maternal_candidate", "family_inherited_paternal_candidate") || continue
        gene = string(get(row, "gene", ""))
        isempty(gene) && continue
        group = get!(grouped, gene, Dict(:maternal => Set{String}(), :paternal => Set{String}()))
        if category == "family_inherited_maternal_candidate"
            push!(group[:maternal], row_variant_key(row))
        else
            push!(group[:paternal], row_variant_key(row))
        end
    end
    selected = Set{String}()
    for gene_groups in values(grouped)
        if !isempty(gene_groups[:maternal]) && !isempty(gene_groups[:paternal])
            union!(selected, gene_groups[:maternal])
            union!(selected, gene_groups[:paternal])
        end
    end
    return selected
end

function identify_one_parent_compound_heterozygous_variants(rows::Vector{Dict{String,Any}}, family::FamilySpec)
    parent = known_parent(family)
    child = singleton_sample(family)
    (parent === nothing || child === nothing) && return Set{String}()
    grouped = Dict{String,Dict{String,Set{String}}}()
    for row in rows
        category = string(get(row, "candidateCategory", ""))
        category in ("one_parent_inherited_candidate", "one_parent_other_parent_candidate") || continue
        genotype_state(row, child) == :het || continue
        gene = string(get(row, "gene", ""))
        isempty(gene) && continue
        group = get!(grouped, gene, Dict("present_in_parent" => Set{String}(), "absent_from_parent" => Set{String}()))
        bucket = category == "one_parent_inherited_candidate" ? "present_in_parent" : "absent_from_parent"
        push!(group[bucket], row_variant_key(row))
    end
    selected = Set{String}()
    for gene_groups in values(grouped)
        if !isempty(gene_groups["present_in_parent"]) && !isempty(gene_groups["absent_from_parent"])
            union!(selected, gene_groups["present_in_parent"])
            union!(selected, gene_groups["absent_from_parent"])
        end
    end
    return selected
end

function load_sample_sexes(family::FamilySpec, pipeline_prefix::Union{Nothing,String}, root_dir::AbstractString)
    isnothing(pipeline_prefix) && return Dict{String,String}()
    metrics_path = joinpath(root_dir, string(pipeline_prefix, "_metrics"), "XY_coverage")
    isfile(metrics_path) || return Dict{String,String}()
    relevant = family_relevant_samples_for_filtering(family)
    sexes = Dict{String,String}()
    open(metrics_path, "r") do io
        for line in eachline(io)
            split_line = split(chomp(line), '\t')
            length(split_line) < 4 && continue
            sample = split_line[1]
            sample in relevant || continue
            sexes[sample] = split_line[4]
        end
    end
    return sexes
end

function family_relevant_samples_for_filtering(family::FamilySpec)
    samples = Set{String}()
    !isnothing(family.parent1) && push!(samples, family.parent1)
    !isnothing(family.parent2) && push!(samples, family.parent2)
    union!(samples, family.affected)
    union!(samples, family.unaffected)
    return samples
end

function classify_x_linked_variant(row::Dict{String,Any}, family::FamilySpec, sample_sexes::Dict{String,String}, gq_hom_cutoff::Float64)
    chrom = uppercase(normalise_chromosome(string(get(row, "chrom", ""))))
    chrom == "X" || return ""
    isempty(family.affected) && return ""
    for sample in family.affected
        sex = get(sample_sexes, sample, "Unknown")
        state = genotype_state(row, sample)
        gq = parse_float(get(row, "GQ ($sample)", ""))
        sex in ("Male", "Female") || return ""
        state == :hom_alt || return ""
        gq >= gq_hom_cutoff || return ""
    end
    for sample in family.unaffected
        genotype_state(row, sample) == :hom_alt && return ""
    end
    mother = x_parent_by_sex(family, sample_sexes, "Female")
    father = x_parent_by_sex(family, sample_sexes, "Male")
    if !isnothing(mother)
        mother_state = genotype_state(row, mother)
        any(get(sample_sexes, sample, "Unknown") == "Male" for sample in family.affected) && !(mother_state in (:het, :hom_alt)) && return ""
        any(get(sample_sexes, sample, "Unknown") == "Female" for sample in family.affected) && !(mother_state in (:het, :hom_alt)) && return ""
    end
    if !isnothing(father) && any(get(sample_sexes, sample, "Unknown") == "Female" for sample in family.affected)
        father_state = genotype_state(row, father)
        father_state in (:het, :hom_alt) || return ""
    end
    return "x_linked_recessive_candidate"
end

function x_parent_by_sex(family::FamilySpec, sample_sexes::Dict{String,String}, sex::String)
    for parent in filter(!isnothing, [family.parent1, family.parent2])
        get(sample_sexes, parent, "Unknown") == sex && return parent
    end
    return nothing
end

function identify_singleton_possible_compound_heterozygous_genes(rows::Vector{Dict{String,Any}}, family::FamilySpec)
    sample = singleton_sample(family)
    sample === nothing && return Set{String}()
    grouped = Dict{String,Set{String}}()
    for row in rows
        string(get(row, "candidateCategory", "")) == "singleton_possible_compound_heterozygous_component" || continue
        genotype_state(row, sample) == :het || continue
        gene = string(get(row, "gene", ""))
        isempty(gene) && continue
        push!(get!(grouped, gene, Set{String}()), row_variant_key(row))
    end
    return Set([gene for (gene, keys) in grouped if length(keys) >= 2])
end

function collapse_transcript_rows(rows::Vector{Dict{String,Any}})
    grouped = Dict{Tuple{String,String},Dict{String,Any}}()
    for row in rows
        key = (string(get(row, "gene", "")), row_variant_key(row))
        current = get(grouped, key, nothing)
        if current === nothing || transcript_priority(row) > transcript_priority(current)
            grouped[key] = row
        elseif transcript_priority(row) == transcript_priority(current) && consequence_rank(row) < consequence_rank(current)
            grouped[key] = row
        end
    end
    return collect(values(grouped))
end

function transcript_priority(row::Dict{String,Any})
    !isempty(string(get(row, "MANE_PLUS_CLINICAL", ""))) && return 4
    !isempty(string(get(row, "MANE_SELECT", ""))) && return 3
    !isempty(string(get(row, "MANE", ""))) && return 2
    return 1
end

function consequence_rank(row::Dict{String,Any})
    coding_effect = string(get(row, "codingEffect", ""))
    consequence = lowercase(string(get(row, "Consequence", "")))
    if occursin("stop_gained", consequence) || coding_effect == "frameshift"
        return 1
    elseif coding_effect == "start loss" || coding_effect == "stop loss"
        return 2
    elseif coding_effect == "missense" || coding_effect == "in-frame"
        return 3
    elseif occursin("splice_", consequence) || has_spliceai_support(row, ThresholdConfig())
        return 4
    elseif coding_effect == "synonymous"
        return 5
    elseif string(get(row, "varLocation", "")) == "intron"
        return 6
    end
    return 7
end

function identify_compound_heterozygous_variants(rows::Vector{Dict{String,Any}}, family::FamilySpec)
    (!isnothing(family.parent1) && !isnothing(family.parent2) && length(family.affected) == 1) || return Set{String}()
    child = family.affected[1]
    grouped = Dict{String,Vector{Dict{String,Any}}}()
    for row in rows
        string(get(row, "candidateCategory", "")) == "inherited_candidate" || continue
        genotype_state(row, child) == :het || continue
        push!(get!(grouped, string(get(row, "gene", "")), Vector{Dict{String,Any}}()), row)
    end
    selected = Set{String}()
    for gene_rows in values(grouped)
        maternal = Dict{String,Dict{String,Any}}()
        paternal = Dict{String,Dict{String,Any}}()
        for row in gene_rows
            origin = inherited_parent(row, family)
            key = row_variant_key(row)
            if origin == :maternal
                maternal[key] = row
            elseif origin == :paternal
                paternal[key] = row
            end
        end
        if !isempty(maternal) && !isempty(paternal)
            union!(selected, keys(maternal))
            union!(selected, keys(paternal))
        end
    end
    return selected
end

function inherited_parent(row::Dict{String,Any}, family::FamilySpec)
    p1_state = genotype_state(row, family.parent1)
    p2_state = genotype_state(row, family.parent2)
    if (p1_state == :het || p1_state == :hom_alt) && p2_state == :ref
        return :maternal
    elseif (p2_state == :het || p2_state == :hom_alt) && p1_state == :ref
        return :paternal
    end
    return :unknown
end

function parse_float(value::Any)
    value isa Real && return Float64(value)
    text = strip(string(value))
    isempty(text) && return NaN
    parsed = tryparse(Float64, text)
    return parsed === nothing ? NaN : parsed
end

function row_variant_key(row::Dict{String,Any})
    chrom = normalise_chromosome(string(get(row, "chrom", "")))
    pos = string(get(row, "inputPos", ""))
    ref = uppercase(string(get(row, "inputRef", "")))
    alt = uppercase(string(get(row, "inputAlt", "")))
    return string(chrom, ":", pos, ":", ref, ":", alt)
end

function normalise_chromosome(chrom::AbstractString)
    startswith(chrom, "chr") ? chrom[4:end] : String(chrom)
end

function denovocnn_proband(family::FamilySpec)
    return isempty(family.affected) ? nothing : family.affected[1]
end

function load_denovocnn_calls(pipeline_prefix::Union{Nothing,String}, root_dir::AbstractString=pwd(), proband::Union{Nothing,String}=nothing)
    isnothing(pipeline_prefix) && return Set{String}()
    directory = joinpath(root_dir, string(pipeline_prefix, "_denovocnn"))
    isdir(directory) || return Set{String}()
    calls = Set{String}()
    for file in readdir(directory; join=true)
        basename_matches = isnothing(proband) || occursin(proband, basename(file))
        basename_matches || continue
        if endswith(file, "_denovos.filtered.txt") || endswith(file, "_denovos") || endswith(file, ".vcf")
            union!(calls, parse_denovocnn_file(file))
        end
    end
    return calls
end

function parse_denovocnn_file(path::AbstractString)
    calls = Set{String}()
    open(path, "r") do io
        for line in eachline(io)
            stripped = strip(line)
            isempty(stripped) && continue
            startswith(stripped, "#") && continue
            fields = split(stripped, '\t')
            key = denovocnn_variant_key(fields)
            key === nothing || push!(calls, key)
        end
    end
    return calls
end

function denovocnn_variant_key(fields::Vector{SubString{String}})
    return denovocnn_variant_key(String.(fields))
end

function denovocnn_variant_key(fields::Vector{String})
    length(fields) < 4 && return nothing
    if occursin(r"^\d+$", fields[2]) && is_dna_allele(fields[3]) && is_dna_allele(fields[4])
        return string(normalise_chromosome(fields[1]), ":", fields[2], ":", uppercase(fields[3]), ":", uppercase(fields[4]))
    end
    if length(fields) >= 5 && occursin(r"^\d+$", fields[2]) && is_dna_allele(fields[4]) && is_dna_allele(fields[5])
        return string(normalise_chromosome(fields[1]), ":", fields[2], ":", uppercase(fields[4]), ":", uppercase(fields[5]))
    end
    return nothing
end

function is_dna_allele(value::AbstractString)
    text = uppercase(strip(value))
    isempty(text) && return false
    text == "." && return false
    return occursin(r"^[ACGTN]+$", text)
end
