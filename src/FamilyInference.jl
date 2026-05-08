Base.@kwdef mutable struct RelatednessRecord
    pair::String = ""
    shared_homozygosity::Float64 = NaN
    potential_homozygosity::Float64 = NaN
    opposite_homozygous_fraction::Float64 = NaN
    relatedness_score::Float64 = NaN
    category::String = ""
end

function infer_family(family::FamilySpec, sample_names::Vector{String}, pipeline_prefix::Union{Nothing,String}, root_dir::AbstractString=pwd())
    if family.shared
        affected = isempty(family.affected) ? copy(sample_names) : copy(family.affected)
        filter!(sample -> !(sample in family.unaffected), affected)
        return FamilySpec(affected=affected, unaffected=copy(family.unaffected), shared=true)
    end

    if !isempty(family.affected) || !isnothing(family.parent1) || !isnothing(family.parent2)
        return family
    end

    if !isnothing(family.parent1) && !isnothing(family.parent2)
        affected = isempty(family.affected) ? [sample for sample in sample_names if sample != family.parent1 && sample != family.parent2 && !(sample in family.unaffected)] : family.affected
        return FamilySpec(parent1=family.parent1, parent2=family.parent2, affected=affected, unaffected=copy(family.unaffected), shared=false)
    end

    control_family = infer_family_from_control(sample_names, root_dir)
    control_family !== nothing && return control_family

    if length(sample_names) == 2
        return FamilySpec(parent1=sample_names[1], parent2=sample_names[2], unaffected=copy(family.unaffected), shared=false)
    end

    relatedness = load_relatedness(sample_names, something(pipeline_prefix, "r04"), root_dir)
    inferred = infer_parents_from_relatedness(sample_names, relatedness)
    if inferred !== nothing
        parent1, parent2 = inferred
        affected = [sample for sample in sample_names if sample != parent1 && sample != parent2 && !(sample in family.unaffected)]
        return FamilySpec(parent1=parent1, parent2=parent2, affected=affected, unaffected=copy(family.unaffected), shared=false)
    end

    if length(sample_names) >= 2
        return FamilySpec(parent1=sample_names[1], parent2=sample_names[2], affected=sample_names[3:end], unaffected=copy(family.unaffected), shared=false)
    end
    return FamilySpec(affected=copy(sample_names), unaffected=copy(family.unaffected), shared=false)
end

function infer_family_from_control(sample_names::Vector{String}, root_dir::AbstractString=pwd())
    control_path = joinpath(root_dir, "control")
    isfile(control_path) || return nothing
    sample_set = Set(sample_names)
    open(control_path, "r") do io
        for line in eachline(io)
            stripped = strip(line)
            startswith(stripped, "TRIO ") || continue
            parts = split(stripped)
            length(parts) == 4 || continue
            trio_set = Set(parts[2:4])
            trio_set == sample_set || continue
            proband, father, mother = parts[2], parts[3], parts[4]
            return FamilySpec(parent1=mother, parent2=father, affected=[proband], shared=false)
        end
    end
    return nothing
end

function load_relatedness(sample_names::Vector{String}, pipeline_prefix::String, root_dir::AbstractString=pwd())
    records = Dict{String,RelatednessRecord}()
    homozygosity_path = joinpath(root_dir, string(pipeline_prefix, "_metrics"), "homozygosity.csv")
    if isfile(homozygosity_path)
        open(homozygosity_path, "r") do io
            for line in eachline(io)
                split_line = split(chomp(line), '\t')
                length(split_line) < 5 && continue
                pair = split_line[1]
                space_pos = findfirst(' ', pair)
                space_pos === nothing && continue
                left = pair[1:(space_pos - 1)]
                right = pair[(space_pos + 1):end]
                (left in sample_names && right in sample_names) || continue
                opposite_hom = tryparse(Float64, split_line[4])
                total_vars = tryparse(Float64, split_line[5])
                opposite_fraction = if opposite_hom === nothing || total_vars === nothing || total_vars == 0.0
                    NaN
                else
                    100.0 * opposite_hom / total_vars
                end
                potential = something(tryparse(Float64, split_line[3]), NaN)
                category = classify_relatedness(opposite_fraction, potential)
                records[pair] = RelatednessRecord(
                    pair=pair,
                    shared_homozygosity=something(tryparse(Float64, split_line[2]), NaN),
                    potential_homozygosity=something(tryparse(Float64, split_line[3]), NaN),
                    opposite_homozygous_fraction=opposite_fraction,
                    category=category,
                )
            end
        end
    end
    relatedness_path = joinpath(root_dir, string(pipeline_prefix, "_metrics"), "relatedness2.csv")
    if isfile(relatedness_path)
        open(relatedness_path, "r") do io
            for line in eachline(io)
                split_line = split(chomp(line), '\t')
                length(split_line) == 7 || continue
                pair = string(split_line[1], " ", split_line[2])
                if haskey(records, pair)
                    value = tryparse(Float64, split_line[7])
                    value === nothing && continue
                    records[pair].relatedness_score = 2.0 * value
                    if records[pair].relatedness_score > 0.75
                        records[pair].category = "Samples from the same person or identical twins"
                    elseif isnan(records[pair].opposite_homozygous_fraction) || (records[pair].shared_homozygosity == 0.0 && records[pair].potential_homozygosity == 0.0)
                        records[pair].category = classify_relatedness_from_score(records[pair].relatedness_score)
                    end
                end
            end
        end
    end
    return records
end

function classify_relatedness(opposite_fraction::Float64, potential_homozygosity::Float64)
    if !isnan(opposite_fraction) && opposite_fraction < 3.0
        return "Parent-child"
    elseif !isnan(potential_homozygosity) && potential_homozygosity < 10.0
        return "Unrelated"
    elseif !isnan(potential_homozygosity) && potential_homozygosity < 15.0
        return "Possibly related"
    elseif !isnan(potential_homozygosity) && potential_homozygosity >= 60.0
        return "Siblings"
    end
    return "Related"
end

function classify_relatedness_from_score(score::Float64)
    isnan(score) && return "Unknown"
    if score >= 0.35
        return "Parent-child or siblings"
    elseif score >= 0.08
        return "Related"
    elseif score > -0.05
        return "Unrelated"
    end
    return "Unrelated"
end

function infer_parents_from_relatedness(sample_names::Vector{String}, relatedness::Dict{String,RelatednessRecord})
    length(sample_names) < 3 && return nothing
    best_pair = nothing
    best_score = -1
    for i in 1:(length(sample_names) - 1)
        for j in (i + 1):length(sample_names)
            sample1 = sample_names[i]
            sample2 = sample_names[j]
            pair_key = sample1 < sample2 ? string(sample1, " ", sample2) : string(sample2, " ", sample1)
            pair_record = get(relatedness, pair_key, nothing)
            pair_record === nothing && continue
            is_partner_pair = startswith(pair_record.category, "Unrelated") || startswith(pair_record.category, "Possibly related") || startswith(pair_record.category, "Related")
            is_partner_pair || continue
            child_links = 0
            for sample in sample_names
                sample == sample1 && continue
                sample == sample2 && continue
                key1 = sample < sample1 ? string(sample, " ", sample1) : string(sample1, " ", sample)
                key2 = sample < sample2 ? string(sample, " ", sample2) : string(sample2, " ", sample)
                record1 = get(relatedness, key1, nothing)
                record2 = get(relatedness, key2, nothing)
                if record1 !== nothing && record2 !== nothing && startswith(record1.category, "Parent-child") && startswith(record2.category, "Parent-child")
                    child_links += 1
                end
            end
            if child_links > best_score
                best_pair = (sample1, sample2)
                best_score = child_links
            end
        end
    end
    return best_pair
end
