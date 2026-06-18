function write_report(path::String, rows::Vector{Dict{String,Any}}, headers::Vector{String}, family::FamilySpec, options::RunOptions, config::AppConfig, sample_names::Vector{String})
    if options.output_format == :html || endswith(lowercase(path), ".html")
        return write_html_report(path, rows, headers, family, options, config, sample_names)
    end
    return write_tsv_report(path, rows, headers)
end

function write_tsv_report(path::String, rows::Vector{Dict{String,Any}}, headers::Vector{String})
    final_headers = ordered_headers(rows, headers)
    open(path, "w") do io
        println(io, join(final_headers, '\t'))
        for row in rows
            values = [escape_tsv(get(row, header, "")) for header in final_headers]
            println(io, join(values, '\t'))
        end
    end
    return path
end

function write_html_report(path::String, rows::Vector{Dict{String,Any}}, headers::Vector{String}, family::FamilySpec, options::RunOptions, config::AppConfig, sample_names::Vector{String})
    pipeline_prefix = something(options.pipeline_prefix, "r04")
    root_dir = batch_root(options.input, options.pipeline_prefix)
    ordered_samples = report_sample_names(sample_names, family)
    sample_details = load_sample_details(config.resources.sample_comments)
    sample_qc = load_sample_qc(ordered_samples, pipeline_prefix, root_dir)
    relatedness = load_relatedness(ordered_samples, pipeline_prefix, root_dir)
    gene_annotations = load_gene_annotations(config.resources)
    annotate_rows!(rows, gene_annotations, sample_qc, family)
    sections = report_sections(rows, family)
    open(path, "w") do io
        println(io, "<!doctype html>")
        println(io, "<html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">")
        println(io, "<title>Variant Prioritisation Report</title>")
        println(io, html_styles())
        println(io, "</head><body>")
        println(io, "<main class=\"page\">")
        println(io, "<section class=\"hero\">")
        println(io, "<div><p class=\"eyebrow\">r04 Variant Prioritisation</p><h1>Family Report</h1><p class=\"sub\">Input: ", html_escape(options.input), "</p></div>")
        println(io, "<div class=\"hero-stats\">")
        println(io, stat_tile("Variants kept", string(length(rows))))
        println(io, stat_tile("Categories", string(length(sections))))
        println(io, stat_tile("Samples", string(length(ordered_samples))))
        println(io, "</div></section>")
        println(io, family_summary_html(family))
        println(io, validation_panels_html(ordered_samples, sample_details, sample_qc, relatedness, family))
        println(io, category_summary_html(sections))
        for section in sections
            println(io, variant_section_html(section.id, section.rows, headers, family, section.title))
        end
        println(io, html_scripts())
        println(io, "</main></body></html>")
    end
    return path
end

function ordered_headers(rows::Vector{Dict{String,Any}}, discovered_headers::Vector{String})
    headers = copy(PRIORITY_HEADERS)
    for row in rows
        for key in keys(row)
            key in headers || push!(headers, key)
        end
    end
    for header in discovered_headers
        header in headers || push!(headers, header)
    end
    return headers
end

function escape_tsv(value::Any)
    text = replace(string(value), '\t' => ' ')
    return replace(text, '\n' => ' ')
end

function load_sample_details(paths::Vector{String})
    details = Dict{String,Vector{String}}()
    for path in paths
        isfile(path) || continue
        open(path, "r") do io
            for line in eachline(io)
                stripped = chomp(line)
                isempty(stripped) && continue
                parts = split(stripped, '\t'; limit=2)
                length(parts) == 2 || continue
                push!(get!(details, parts[1], String[]), parts[2])
            end
        end
        break
    end
    return details
end

struct ReportSection
    id::String
    title::String
    rows::Vector{Dict{String,Any}}
end

function report_sections(rows::Vector{Dict{String,Any}}, family::FamilySpec)
    sections = ReportSection[]
    push_section!(sections, "recessive_homozygous_candidate", "Homozygous Variants", [row for row in rows if string(get(row, "candidateCategory", "")) == "recessive_homozygous_candidate"])
    push_section!(sections, "x_linked_recessive_candidate", "X-Linked Variants", [row for row in rows if string(get(row, "candidateCategory", "")) == "x_linked_recessive_candidate"])
    singleton_mode = !isnothing(singleton_sample(family)) && isnothing(family.parent1) && isnothing(family.parent2)
    one_parent_mode = known_parent_count(family) == 1 && length(family.affected) == 1
    comphet_title = singleton_mode || one_parent_mode ? "Variants Which Could Be Compound Heterozygous" : "Compound Heterozygous Variants"
    comphet_rows = singleton_mode ?
        [row for row in rows if string(get(row, "candidateCategory", "")) == "singleton_possible_compound_heterozygous_candidate"] :
        [row for row in rows if string(get(row, "candidateCategory", "")) == "compound_heterozygous_candidate"]
    push_section!(sections, singleton_mode ? "singleton_possible_compound_heterozygous_candidate" : "compound_heterozygous_candidate", comphet_title, comphet_rows)
    push_section!(sections, "de_novo_candidate", "De Novo Variants", [row for row in rows if string(get(row, "candidateCategory", "")) == "de_novo_candidate"])
    segregating_rows = [row for row in rows if string(get(row, "candidateCategory", "")) in ("segregating_family_candidate", "recessive_family_candidate")]
    push_section!(sections, "segregating_family_candidate", "Segregating Variants", segregating_rows)
    push_section!(sections, "imprinted_gene_flag", "Variants In Imprinted Genes", [row for row in rows if get(row, "is_imprinted_gene", false) == true])
    push_section!(sections, "manta_deletion_candidate", "Manta Deletions", [row for row in rows if string(get(row, "candidateCategory", "")) == "manta_deletion_candidate"])
    push_section!(sections, "manta_duplication_candidate", "Manta Duplications", [row for row in rows if string(get(row, "candidateCategory", "")) == "manta_duplication_candidate"])
    push_section!(sections, "manta_insertion_candidate", "Manta Insertions", [row for row in rows if string(get(row, "candidateCategory", "")) == "manta_insertion_candidate"])
    return sections
end

function push_section!(sections::Vector{ReportSection}, id::String, title::String, rows::Vector{Dict{String,Any}})
    isempty(rows) && return
    sort!(rows, by=report_row_sort_key)
    push!(sections, ReportSection(id, title, rows))
end

function family_summary_html(family::FamilySpec)
    lines = String[]
    !isempty(family.affected) && push!(lines, "<strong>Proband:</strong> " * html_escape(join(family.affected, ", ")))
    !isnothing(family.parent1) && push!(lines, "<strong>Parent 1:</strong> " * html_escape(family.parent1))
    !isnothing(family.parent2) && push!(lines, "<strong>Parent 2:</strong> " * html_escape(family.parent2))
    !isempty(family.unaffected) && push!(lines, "<strong>Unaffected:</strong> " * html_escape(join(collect(family.unaffected), ", ")))
    family.shared && push!(lines, "<strong>Mode:</strong> Shared-variant analysis")
    return "<section class=\"panel summary\"><h2>Family Summary</h2><div class=\"summary-grid\">" * join(["<div class=\"summary-item\">$line</div>" for line in lines], "") * "</div></section>"
end

function sample_cards_html(sample_names::Vector{String}, sample_details::Dict{String,Vector{String}}, family::FamilySpec)
    cards = String[]
    for sample in sample_names
        role = sample_role(sample, family)
        details = get(sample_details, sample, ["No sample details available"])
        body = join(["<li>" * html_escape(detail) * "</li>" for detail in details], "")
        push!(cards, "<article class=\"sample-card\"><div class=\"sample-head\"><h3>" * html_escape(sample) * "</h3><span class=\"role\">" * html_escape(role) * "</span></div><ul>$body</ul></article>")
    end
    return "<section class=\"panel\"><h2>Sample Details</h2><div class=\"sample-grid\">" * join(cards, "") * "</div></section>"
end

function category_summary_html(sections)
    cards = String[]
    for section in sections
        push!(cards, "<a class=\"category-card\" href=\"#" * html_escape(section.id) * "\"><span class=\"category-name\">" * html_escape(section.title) * "</span><strong>" * string(length(section.rows)) * "</strong></a>")
    end
    return "<section class=\"panel\"><h2>Variant Classes</h2><div class=\"category-grid\">" * join(cards, "") * "</div></section>"
end

function variant_section_html(category::String, rows::Vector{Dict{String,Any}}, headers::Vector{String}, family::FamilySpec, title::Union{Nothing,String}=nothing)
    table_headers = variant_table_headers(headers, family, rows)
    section_title = isnothing(title) ? pretty_category(category) : title
    table_id = category * "_table"
    impact_filter = "IMPACT" in table_headers
    controls = impact_filter ?
        "<label class=\"impact-filter\"><span>Impact</span><select class=\"impact-filter-select\" data-table-id=\"" * html_escape(table_id) * "\"><option value=\"\">All</option><option value=\"HIGH\">HIGH</option><option value=\"MODERATE\">MODERATE</option><option value=\"LOW\">LOW</option><option value=\"MODIFIER\">MODIFIER</option><option value=\"__blank__\">Blank</option></select></label>" :
        ""
    table = ["<section class=\"panel category-section\" id=\"" * html_escape(category) * "\"><div class=\"section-head\"><div><h2>" * html_escape(section_title) * "</h2><span class=\"pill\" data-table-count=\"" * html_escape(table_id) * "\">" * string(length(rows)) * " variants</span></div><div class=\"section-controls\">" * controls * "<button type=\"button\" class=\"copy-button\" data-table-id=\"" * html_escape(table_id) * "\">Copy Table</button></div></div>"]
    push!(table, "<div class=\"table-wrap\"><table id=\"" * html_escape(table_id) * "\"><thead><tr>")
    push!(table, "<th>Assessed</th>")
    for header in table_headers
        push!(table, "<th>" * html_escape(header) * "</th>")
    end
    push!(table, "</tr></thead><tbody>")
    for row in rows
        row_id = assessed_row_id(category, row)
        row_impact = uppercase(string(get(row, "IMPACT", "")))
        push!(table, "<tr data-row-id=\"" * html_escape(row_id) * "\" data-impact=\"" * html_escape(row_impact) * "\">")
        push!(table, "<td><input type=\"checkbox\" class=\"assessed-checkbox\" data-row-id=\"" * html_escape(row_id) * "\" aria-label=\"Mark variant assessed\"></td>")
        for header in table_headers
            push!(table, "<td>" * html_escape(display_value(row, header)) * "</td>")
        end
        push!(table, "</tr>")
    end
    push!(table, "</tbody></table></div></section>")
    return join(table, "")
end

function variant_table_headers(headers::Vector{String}, family::FamilySpec, rows::Vector{Dict{String,Any}})
    if all(row -> haskey(row, "svType"), rows)
        return sv_variant_table_headers(headers, family, rows)
    end
    return small_variant_table_headers(headers, family)
end

function small_variant_table_headers(headers::Vector{String}, family::FamilySpec)
    selected = ["gene", "omim_annotations", "panelapp_status", "panelapp_phenotypes", "transcript", "IMPACT", "Consequence", "gNomen", "cNomen", "pNomen", "MANE_SELECT", "GnomAD_v4_1_AF_all", "GnomAD_v4_1_AF_popmax", "spliceai_summary", "clinVarClinSignifs", "imprinting_status", "roh_overlap"]
    sample_headers = String[]
    for sample in report_sample_names(headers_to_samples(headers), family)
        append!(sample_headers, ["GT ($sample)", "AD ($sample)", "GQ ($sample)"])
    end
    for header in sample_headers
        header in headers || continue
        header in selected || push!(selected, header)
    end
    return selected
end

function sv_variant_table_headers(headers::Vector{String}, family::FamilySpec, rows::Vector{Dict{String,Any}})
    selected = ["gene", "omim_annotations", "panelapp_status", "panelapp_phenotypes", "svType", "SVLEN", "inheritance_model", "transcript", "IMPACT", "Consequence", "gNomen", "cNomen", "pNomen", "MANE_SELECT", "imprinting_status"]
    for sample in report_sample_names(headers_to_samples(headers), family)
        append!(selected, ["GT ($sample)", "FT ($sample)", "GQ ($sample)", "PR ($sample)", "SR ($sample)"])
    end
    return [header for header in selected if header in headers || any(row -> haskey(row, header), rows)]
end

function headers_to_samples(headers::Vector{String})
    return [header[5:end-1] for header in headers if startswith(header, "GT (")]
end

function display_value(row::Dict{String,Any}, header::String)
    value = get(row, header, "")
    if header == "GnomAD_v4_1_AF_all"
        parsed = tryparse(Float64, string(value))
        parsed === nothing && return string(value)
        return parsed == 0 ? "0" : string(round(parsed, sigdigits=3))
    end
    return string(value)
end

Base.@kwdef mutable struct SampleQc
    sample::String = ""
    sex::String = "Unknown"
    abnormal_sex::Bool = false
    read_count::Union{Nothing,Int} = nothing
    read_count_is_raw::Bool = false
    mean_depth::Union{Nothing,Float64} = nothing
    coverage20::Union{Nothing,Float64} = nothing
    vdp::Union{Nothing,Float64} = nothing
    exomedepth_correlation::Union{Nothing,Float64} = nothing
    cnv_count::Union{Nothing,Int} = nothing
    contamination::Union{Nothing,Float64} = nothing
    homozygosity::Union{Nothing,Float64} = nothing
    checkfastq::String = ""
    homozygous_regions::Vector{NamedTuple{(:chrom, :start, :stop),Tuple{String,Int,Int}}} = NamedTuple{(:chrom, :start, :stop),Tuple{String,Int,Int}}[]
end

Base.@kwdef struct GeneAnnotations
    omim::Dict{String,Vector{String}} = Dict{String,Vector{String}}()
    panelapp::Dict{String,Vector{String}} = Dict{String,Vector{String}}()
    imprinted::Dict{String,String} = Dict{String,String}()
end

function validation_panels_html(sample_names::Vector{String}, sample_details::Dict{String,Vector{String}}, sample_qc::Dict{String,SampleQc}, relatedness::Dict{String,RelatednessRecord}, family::FamilySpec)
    return "<section class=\"panel\"><h2>Sample Details</h2><div class=\"sample-grid\">" *
        join([sample_card_html(sample, sample_details, sample_qc, family) for sample in sample_names], "") *
        "</div></section>" *
        relatedness_panel_html(relatedness, sample_names)
end

function sample_display_order(sample_names::Vector{String}, family::FamilySpec)
    ordered = String[]
    append_if_present!(ordered, family.affected)
    append_if_present!(ordered, filter(!isnothing, [family.parent1, family.parent2]))
    append_if_present!(ordered, [sample for sample in sample_names if !(sample in ordered)])
    return ordered
end

function report_sample_names(sample_names::Vector{String}, family::FamilySpec)
    relevant = family_relevant_samples(family)
    filtered = isempty(relevant) ? sample_names : [sample for sample in sample_names if sample in relevant]
    return sample_display_order(filtered, family)
end

function family_relevant_samples(family::FamilySpec)
    samples = Set{String}()
    !isnothing(family.parent1) && push!(samples, family.parent1)
    !isnothing(family.parent2) && push!(samples, family.parent2)
    union!(samples, family.affected)
    union!(samples, family.unaffected)
    return samples
end

function append_if_present!(target::Vector{String}, values)
    for value in values
        sample = string(value)
        isempty(sample) && continue
        sample in target || push!(target, sample)
    end
    return target
end

function sample_card_html(sample::String, sample_details::Dict{String,Vector{String}}, sample_qc::Dict{String,SampleQc}, family::FamilySpec)
    role = sample_role(sample, family)
    details = get(sample_details, sample, ["No sample details available"])
    qc = get(sample_qc, sample, SampleQc(sample=sample))
    detail_list = join(["<li>" * html_escape(detail) * "</li>" for detail in details], "")
    qc_list = join(filter(!isempty, [
        qc.sex == "Unknown" ? "" : "<li><strong>Sex:</strong> " * html_escape(qc.sex * (qc.abnormal_sex ? " (Abnormal)" : "")) * "</li>",
        isnothing(qc.read_count) ? "" : "<li><strong>Read count:</strong> " * string(qc.read_count) * (qc.read_count_is_raw ? " raw" : " dedup") * "</li>",
        isnothing(qc.mean_depth) ? "" : "<li><strong>Mean depth:</strong> " * string(round(qc.mean_depth, digits=1)) * "</li>",
        isnothing(qc.coverage20) ? "" : "<li><strong>Coverage @20X:</strong> " * string(round(qc.coverage20, digits=1)) * "%</li>",
        isnothing(qc.vdp) ? "" : "<li><strong>VDP:</strong> " * string(round(qc.vdp, digits=1)) * "%</li>",
        isnothing(qc.exomedepth_correlation) ? "" : "<li><strong>ExomeDepth correlation:</strong> " * string(round(qc.exomedepth_correlation * 100, digits=2)) * "</li>",
        isnothing(qc.cnv_count) ? "" : "<li><strong>CNV count:</strong> " * string(qc.cnv_count) * "</li>",
        isnothing(qc.contamination) ? "" : "<li><strong>Contamination:</strong> " * string(round(qc.contamination * 100, digits=2)) * "%</li>",
        isnothing(qc.homozygosity) ? "" : "<li><strong>Homozygosity:</strong> " * string(round(qc.homozygosity, digits=2)) * "%</li>",
        isempty(qc.checkfastq) ? "" : "<li><strong>FASTQ checks:</strong> " * html_escape(qc.checkfastq) * "</li>",
    ]))
    return "<article class=\"sample-card\"><div class=\"sample-head\"><h3>" * html_escape(sample) * "</h3><span class=\"role\">" * html_escape(role) * "</span></div><ul>$detail_list$qc_list</ul></article>"
end

function relatedness_panel_html(relatedness::Dict{String,RelatednessRecord}, sample_names::Vector{String})
    isempty(relatedness) && return ""
    rows = String[]
    for pair in sort(collect(keys(relatedness)))
        record = relatedness[pair]
        parts = split(pair, ' ')
        length(parts) == 2 || continue
        (parts[1] in sample_names && parts[2] in sample_names) || continue
        push!(rows,
            "<tr><td>" * html_escape(parts[1]) * "</td><td>" * html_escape(parts[2]) * "</td><td>" * html_escape(isnan(record.relatedness_score) ? "" : string(round(record.relatedness_score, digits=3))) *
            "</td><td>" * html_escape(isnan(record.shared_homozygosity) ? "" : string(round(record.shared_homozygosity, digits=2))) *
            "</td><td>" * html_escape(isnan(record.potential_homozygosity) ? "" : string(round(record.potential_homozygosity, digits=2))) *
            "</td><td>" * html_escape(isnan(record.opposite_homozygous_fraction) ? "" : string(round(record.opposite_homozygous_fraction, digits=2))) *
            "</td><td>" * html_escape(record.category) * "</td></tr>")
    end
    isempty(rows) && return ""
    return "<section class=\"panel\"><h2>Relatedness And Parent-Child Validation</h2><div class=\"table-wrap\"><table><thead><tr><th>Sample 1</th><th>Sample 2</th><th>Relatedness</th><th>Shared Homozygosity</th><th>Potential Homozygosity</th><th>Opposite Homozygous %</th><th>Likely Relationship</th></tr></thead><tbody>" * join(rows, "") * "</tbody></table></div></section>"
end

function load_sample_qc(sample_names::Vector{String}, pipeline_prefix::String, root_dir::AbstractString=pwd())
    qc = Dict(name => SampleQc(sample=name) for name in sample_names)
    metrics_dir = joinpath(root_dir, string(pipeline_prefix, "_metrics"))
    load_xy_coverage!(qc, joinpath(metrics_dir, "XY_coverage"))
    load_coverage!(qc, joinpath(metrics_dir, "coverage_report"))
    load_exomedepth!(qc, joinpath(metrics_dir, "ExomeDepth_stats.csv"))
    load_contamination!(qc, metrics_dir)
    load_checkfastq!(qc, metrics_dir)
    load_homozygosity_summary!(qc, joinpath(metrics_dir, "homozygosity.csv"))
    load_homozygous_regions!(qc, metrics_dir)
    return qc
end

function load_xy_coverage!(qc::Dict{String,SampleQc}, path::String)
    isfile(path) || return
    open(path, "r") do io
        for line in eachline(io)
            split_line = split(chomp(line), '\t')
            length(split_line) < 3 && continue
            sample = get(qc, split_line[1], nothing)
            sample === nothing && continue
            if length(split_line) > 3
                sample.sex = split_line[4]
            end
            if length(split_line) > 4
                sample.abnormal_sex = split_line[5] == "Abnormal"
            end
        end
    end
end

function load_coverage!(qc::Dict{String,SampleQc}, path::String)
    isfile(path) || return
    open(path, "r") do io
        for line in eachline(io)
            split_line = split(chomp(line), '\t')
            length(split_line) < 3 && continue
            sample = get(qc, split_line[1], nothing)
            sample === nothing && continue
            sample.mean_depth = tryparse(Float64, split_line[2])
            sample.coverage20 = tryparse(Float64, split_line[3])
            if length(split_line) >= 4
                sample.vdp = tryparse(Float64, split_line[4])
            end
        end
    end
end

function load_exomedepth!(qc::Dict{String,SampleQc}, path::String)
    isfile(path) || return
    open(path, "r") do io
        for line in eachline(io)
            split_line = split(chomp(line), '\t')
            length(split_line) < 3 && continue
            sample = get(qc, split_line[1], nothing)
            sample === nothing && continue
            sample.exomedepth_correlation = tryparse(Float64, split_line[2])
            sample.cnv_count = tryparse(Int, split_line[3])
        end
    end
end

function load_contamination!(qc::Dict{String,SampleQc}, metrics_dir::AbstractString)
    for sample in keys(qc)
        path = joinpath(metrics_dir, string(sample, "_cleanCall.csv"))
        isfile(path) || continue
        lines = readlines(path)
        length(lines) < 2 && continue
        split_line = split(chomp(lines[2]), '\t')
        length(split_line) >= 7 || continue
        qc[sample].contamination = tryparse(Float64, split_line[7])
    end
end

function load_checkfastq!(qc::Dict{String,SampleQc}, metrics_dir::AbstractString)
    for sample in keys(qc)
        path = joinpath(metrics_dir, string(sample, "_checkFastq.txt"))
        isfile(path) || continue
        notes = String[]
        open(path, "r") do io
            for line in eachline(io)
                stripped = chomp(line)
                parsed = tryparse(Int, stripped)
                if parsed === nothing
                    push!(notes, stripped)
                else
                    qc[sample].read_count = parsed
                    qc[sample].read_count_is_raw = true
                end
            end
        end
        qc[sample].checkfastq = join(notes, "; ")
    end
end

function load_homozygosity_summary!(qc::Dict{String,SampleQc}, path::String)
    isfile(path) || return
    open(path, "r") do io
        for line in eachline(io)
            split_line = split(chomp(line), '\t')
            isempty(split_line) && continue
            sample = get(qc, split_line[1], nothing)
            sample === nothing && continue
            length(split_line) > 1 && (sample.homozygosity = tryparse(Float64, split_line[2]))
        end
    end
end

function load_homozygous_regions!(qc::Dict{String,SampleQc}, metrics_dir::AbstractString)
    for sample in keys(qc)
        path = joinpath(metrics_dir, string(sample, "_homozygosity.csv"))
        isfile(path) || continue
        open(path, "r") do io
            for line in eachline(io)
                split_line = split(chomp(line), '\t')
                length(split_line) == 4 || continue
                start_pos = tryparse(Int, split_line[2])
                stop_pos = tryparse(Int, split_line[3])
                (start_pos === nothing || stop_pos === nothing) && continue
                push!(qc[sample].homozygous_regions, (chrom=normalise_chromosome(split_line[1]), start=start_pos, stop=stop_pos))
            end
        end
    end
end

function load_gene_annotations(resources::ResourceConfig)
    return GeneAnnotations(
        omim=load_omim_annotations(resources),
        panelapp=load_panelapp_annotations(resources.panelapp_genes),
        imprinted=load_imprinted_genes(resources.imprinted_genes),
    )
end

function load_omim_annotations(resources::ResourceConfig)
    annotations = Dict{String,Vector{String}}()
    ensg_to_hgnc = load_ensg_to_hgnc(resources.omim_ensg_to_hgnc)
    mim_to_ensg = load_mim_to_ensg(resources.omim_mim2gene)
    for path in resources.omim_morbid
        isfile(path) || continue
        open(path, "r") do io
            for line in eachline(io)
                startswith(line, "#") && continue
                split_line = split(chomp(line), '\t')
                length(split_line) == 4 || continue
                hgnc = get(ensg_to_hgnc, get(mim_to_ensg, split_line[3], ""), nothing)
                hgnc === nothing && continue
                push!(get!(annotations, hgnc, String[]), "OMIM: " * split_line[1])
            end
        end
        break
    end
    return annotations
end

function load_panelapp_annotations(paths::Vector{String})
    annotations = Dict{String,Vector{String}}()
    for path in paths
        isfile(path) || continue
        open(path, "r") do io
            for line in eachline(io)
                split_line = split(chomp(line), '\t')
                length(split_line) >= 5 || continue
                split_line[3] == "3" || continue
                annotation = "Green; " * split_line[5]
                for key in (split_line[1], split_line[2], replace(split_line[1], "HGNC:" => ""))
                    push!(get!(annotations, key, String[]), annotation)
                end
            end
        end
        break
    end
    return annotations
end

function load_ensg_to_hgnc(paths::Vector{String})
    map = Dict{String,String}()
    for path in paths
        isfile(path) || continue
        open(path, "r") do io
            for line in eachline(io)
                split_line = split(chomp(line), '\t')
                length(split_line) == 2 || continue
                startswith(split_line[2], "HGNC:") || continue
                map[split_line[1]] = replace(split_line[2], "HGNC:" => "")
            end
        end
        break
    end
    return map
end

function load_mim_to_ensg(paths::Vector{String})
    map = Dict{String,String}()
    for path in paths
        isfile(path) || continue
        open(path, "r") do io
            for line in eachline(io)
                split_line = split(chomp(line), '\t')
                length(split_line) == 5 || continue
                occursin("gene", split_line[2]) || continue
                map[split_line[1]] = split_line[5]
            end
        end
        break
    end
    return map
end

function load_imprinted_genes(paths::Vector{String})
    annotations = Dict{String,String}()
    for path in paths
        isfile(path) || continue
        open(path, "r") do io
            for line in eachline(io)
                split_line = split(chomp(line), '\t')
                length(split_line) >= 2 || continue
                annotations[split_line[1]] = split_line[2]
            end
        end
        break
    end
    return annotations
end

function annotate_rows!(rows::Vector{Dict{String,Any}}, gene_annotations::GeneAnnotations, sample_qc::Dict{String,SampleQc}, family::FamilySpec)
    affected_samples = !isempty(family.affected) ? family.affected : String[]
    for row in rows
        gene = string(get(row, "gene", ""))
        gene_id = string(get(row, "geneId", ""))
        omim = resolve_gene_annotations(gene_annotations.omim, gene, gene_id)
        panelapp = resolve_gene_annotations(gene_annotations.panelapp, gene, gene_id)
        row["omim_annotations"] = isempty(omim) ? "" : join(unique(omim), "; ")
        row["panelapp_status"] = isempty(panelapp) ? "" : "PanelApp Green"
        row["panelapp_phenotypes"] = isempty(panelapp) ? "" : join(unique(panelapp), "; ")
        combined = vcat(omim, panelapp)
        row["gene_annotations"] = isempty(combined) ? "" : join(unique(combined), "; ")
        imprint = get(gene_annotations.imprinted, gene, nothing)
        row["is_imprinted_gene"] = imprint !== nothing
        row["imprinting_status"] = imprint === nothing ? "" : imprint
        overlaps = String[]
        for sample in affected_samples
            genotype_state(row, sample) == :hom_alt || continue
            overlaps_roh(row, get(sample_qc, sample, SampleQc(sample=sample))) || continue
            push!(overlaps, sample)
        end
        row["roh_overlap"] = isempty(overlaps) ? "" : "Within ROH: " * join(overlaps, ", ")
    end
end

function resolve_gene_annotations(annotations::Dict{String,Vector{String}}, gene::String, gene_id::String)
    values = String[]
    for key in (gene, gene_id, replace(gene_id, "HGNC:" => ""))
        isempty(key) && continue
        append!(values, get(annotations, key, String[]))
    end
    return unique(values)
end

function overlaps_roh(row::Dict{String,Any}, sample_qc::SampleQc)
    chrom = normalise_chromosome(string(get(row, "chrom", "")))
    pos = something(tryparse(Int, string(get(row, "inputPos", ""))), 0)
    pos == 0 && return false
    for region in sample_qc.homozygous_regions
        if region.chrom == chrom && region.start <= pos <= region.stop
            return true
        end
    end
    return false
end

function sample_role(sample::String, family::FamilySpec)
    sample == family.parent1 && return "Parent 1"
    sample == family.parent2 && return "Parent 2"
    sample in family.affected && return "Affected"
    sample in family.unaffected && return "Unaffected"
    return "Sample"
end

function pretty_category(category::String)
    return replace(category, "_" => " ")
end

function stat_tile(label::String, value::String)
    return "<div class=\"stat\"><span>" * html_escape(label) * "</span><strong>" * html_escape(value) * "</strong></div>"
end

function report_row_sort_key(row::Dict{String,Any})
    chrom = normalise_chromosome(string(get(row, "chrom", "")))
    chrom_rank = chromosome_sort_rank(chrom)
    pos = something(tryparse(Int, string(get(row, "inputPos", "0"))), 0)
    ref = string(get(row, "inputRef", ""))
    alt = string(get(row, "inputAlt", ""))
    gene = string(get(row, "gene", ""))
    return (chrom_rank, pos, chrom, ref, alt, gene)
end

function chromosome_sort_rank(chrom::AbstractString)
    text = uppercase(strip(String(chrom)))
    parsed = tryparse(Int, text)
    parsed !== nothing && return parsed
    text == "X" && return 23
    text == "Y" && return 24
    text == "M" && return 25
    text == "MT" && return 25
    return 99
end

function assessed_row_id(category::String, row::Dict{String,Any})
    return category * "|" * row_variant_key(row) * "|" * string(get(row, "gene", ""))
end

function html_escape(value::AbstractString)
    text = replace(value, "&" => "&amp;")
    text = replace(text, "<" => "&lt;")
    text = replace(text, ">" => "&gt;")
    text = replace(text, "\"" => "&quot;")
    return text
end

html_escape(value::Any) = html_escape(string(value))

function html_styles()
    return """
<style>
:root{--bg:#f4f0e8;--ink:#1e1c1a;--muted:#6f685f;--panel:#fffdf8;--line:#d8cfbf;--accent:#0f766e;--accent-2:#c2410c;--shadow:0 12px 30px rgba(30,28,26,.08)}
*{box-sizing:border-box} body{margin:0;background:radial-gradient(circle at top left,#f8f3ea, #efe7d7 55%, #e9deca);color:var(--ink);font-family:Georgia, 'Times New Roman', serif}
.page{max-width:1400px;margin:0 auto;padding:32px 20px 64px}
.hero{display:grid;grid-template-columns:1.6fr 1fr;gap:20px;align-items:end;margin-bottom:24px}
.eyebrow{text-transform:uppercase;letter-spacing:.12em;font:600 12px/1.2 Helvetica, Arial, sans-serif;color:var(--accent)}
h1,h2,h3{margin:0 0 10px} h1{font-size:48px;line-height:1} h2{font-size:28px} h3{font-size:20px}
.sub{color:var(--muted);margin:0}
.hero-stats,.sample-grid,.category-grid,.summary-grid{display:grid;gap:14px}
.hero-stats{grid-template-columns:repeat(3,1fr)}
.panel{background:rgba(255,253,248,.92);backdrop-filter:blur(6px);border:1px solid var(--line);border-radius:22px;box-shadow:var(--shadow);padding:20px;margin-bottom:18px}
.stat,.sample-card,.category-card,.summary-item{background:linear-gradient(180deg,#fff, #f9f3e8);border:1px solid var(--line);border-radius:18px;padding:16px}
.stat span,.role{display:block;font:600 12px/1.2 Helvetica, Arial, sans-serif;text-transform:uppercase;letter-spacing:.08em;color:var(--muted)}
.stat strong{display:block;font-size:32px;margin-top:10px}
.summary-grid{grid-template-columns:repeat(auto-fit,minmax(220px,1fr))}
.sample-grid{grid-template-columns:repeat(auto-fit,minmax(260px,1fr))}
.sample-head{display:flex;justify-content:space-between;gap:12px;align-items:flex-start}
.sample-card ul{margin:12px 0 0 18px;padding:0}
.category-grid{grid-template-columns:repeat(auto-fit,minmax(180px,1fr))}
.category-card{text-decoration:none;color:inherit;display:flex;justify-content:space-between;align-items:center}
.category-name{text-transform:capitalize}
.section-head{display:flex;justify-content:space-between;align-items:center;gap:12px;margin-bottom:14px}
.section-controls{display:flex;align-items:center;gap:10px;flex-wrap:wrap}
.copy-button{border:1px solid var(--line);background:linear-gradient(180deg,#fff,#f6ecdd);color:var(--ink);border-radius:999px;padding:10px 14px;font:600 12px/1.2 Helvetica, Arial, sans-serif;letter-spacing:.04em;text-transform:uppercase;cursor:pointer}
.copy-button:hover{border-color:var(--accent);color:var(--accent)}
.impact-filter{display:flex;align-items:center;gap:8px;font:600 12px/1.2 Helvetica, Arial, sans-serif;letter-spacing:.04em;text-transform:uppercase;color:var(--muted)}
.impact-filter-select{border:1px solid var(--line);background:#fffdfa;border-radius:999px;padding:9px 12px;color:var(--ink);font:600 12px/1.2 Helvetica, Arial, sans-serif}
.pill{font:600 12px/1.2 Helvetica, Arial, sans-serif;padding:8px 10px;border-radius:999px;background:#e6f4f1;color:var(--accent)}
.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:16px}
table{width:100%;border-collapse:collapse;background:#fffdfa}
th,td{padding:10px 12px;border-bottom:1px solid #ece4d8;vertical-align:top}
th{position:sticky;top:0;background:#f8f1e4;font:600 12px/1.2 Helvetica, Arial, sans-serif;letter-spacing:.04em;text-transform:uppercase;text-align:left}
tr:nth-child(even) td{background:#fffcf6}
.assessed-checkbox{width:18px;height:18px;accent-color:var(--accent)}
tr.assessed td{background:#dff3ec !important}
tr.assessed:hover td{background:#d7efe7 !important}
@media (max-width: 900px){.hero{grid-template-columns:1fr}.hero-stats{grid-template-columns:1fr 1fr}.page{padding:20px 14px 48px}h1{font-size:36px}}
</style>
"""
end

function html_scripts()
    return """
<script>
const assessedStorageKey = "variant-prioritiser-assessed";

function loadAssessedRows() {
  try {
    return JSON.parse(localStorage.getItem(assessedStorageKey) || "{}");
  } catch (_) {
    return {};
  }
}

function saveAssessedRows(state) {
  localStorage.setItem(assessedStorageKey, JSON.stringify(state));
}

function applyAssessedState() {
  const state = loadAssessedRows();
  document.querySelectorAll(".assessed-checkbox").forEach((checkbox) => {
    const rowId = checkbox.dataset.rowId;
    const checked = !!state[rowId];
    checkbox.checked = checked;
    const row = checkbox.closest("tr");
    if (row) {
      row.classList.toggle("assessed", checked);
    }
    checkbox.addEventListener("change", () => {
      const next = loadAssessedRows();
      next[rowId] = checkbox.checked;
      saveAssessedRows(next);
      if (row) {
        row.classList.toggle("assessed", checkbox.checked);
      }
    });
  });
}

function setCopyButtonState(button, label) {
  if (!button) return;
  const original = button.dataset.originalLabel || button.textContent;
  button.dataset.originalLabel = original;
  button.textContent = label;
  setTimeout(() => {
    button.textContent = button.dataset.originalLabel || original;
  }, 1200);
}

function fallbackCopyText(text) {
  const textarea = document.createElement("textarea");
  textarea.value = text;
  textarea.setAttribute("readonly", "");
  textarea.style.position = "fixed";
  textarea.style.top = "-1000px";
  textarea.style.left = "-1000px";
  document.body.appendChild(textarea);
  textarea.focus();
  textarea.select();
  textarea.setSelectionRange(0, textarea.value.length);
  try {
    return document.execCommand("copy");
  } finally {
    document.body.removeChild(textarea);
  }
}

function copyTable(tableId, button) {
  const table = document.getElementById(tableId);
  if (!table) return;
  const rows = [...table.querySelectorAll("tr")].map((tr) => {
    const cells = [...tr.querySelectorAll("th,td")];
    return cells
      .filter((cell, index) => index !== 0)
      .map((cell) => cell.innerText.replace(/\\t/g, " ").replace(/\\n/g, " ").trim())
      .join("\\t");
  });
  const text = rows.join("\\n");
  if (navigator.clipboard && window.isSecureContext) {
    navigator.clipboard.writeText(text).then(() => {
      setCopyButtonState(button, "Copied");
    }).catch(() => {
      if (fallbackCopyText(text)) {
        setCopyButtonState(button, "Copied");
      } else {
        setCopyButtonState(button, "Press Ctrl+C");
      }
    });
    return;
  }
  if (fallbackCopyText(text)) {
    setCopyButtonState(button, "Copied");
  } else {
    setCopyButtonState(button, "Press Ctrl+C");
  }
}

function bindCopyButtons() {
  document.querySelectorAll(".copy-button").forEach((button) => {
    button.addEventListener("click", () => copyTable(button.dataset.tableId, button));
  });
}

function updateVisibleCount(tableId) {
  const table = document.getElementById(tableId);
  const pill = document.querySelector('[data-table-count="' + tableId + '"]');
  if (!table || !pill) return;
  const visible = [...table.querySelectorAll("tbody tr")].filter((row) => row.style.display !== "none").length;
  pill.textContent = visible + " variants";
}

function applyImpactFilter(tableId, selectedImpact) {
  const table = document.getElementById(tableId);
  if (!table) return;
  table.querySelectorAll("tbody tr").forEach((row) => {
    const rowImpact = (row.dataset.impact || "").toUpperCase();
    const show =
      selectedImpact === "" ||
      (selectedImpact === "__blank__" ? rowImpact === "" : rowImpact === selectedImpact);
    row.style.display = show ? "" : "none";
  });
  updateVisibleCount(tableId);
}

function bindImpactFilters() {
  document.querySelectorAll(".impact-filter-select").forEach((select) => {
    applyImpactFilter(select.dataset.tableId, select.value);
    select.addEventListener("change", () => applyImpactFilter(select.dataset.tableId, select.value));
  });
}

document.addEventListener("DOMContentLoaded", () => {
  applyAssessedState();
  bindCopyButtons();
  bindImpactFilters();
});
</script>
"""
end
