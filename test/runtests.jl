using Test
include(joinpath(@__DIR__, "..", "src", "VariantPrioritiser.jl"))
const VP = VariantPrioritiser

@testset "Singleton heterozygous inclusion" begin
    family = VP.FamilySpec(affected=["child"])
    row = Dict{String,Any}("chrom" => "1", "inputPos" => "100", "inputRef" => "A",
        "inputAlt" => "G", "gene" => "GENE1", "GT (child)" => "0/1",
        "GQ (child)" => "30", "IMPACT" => "MODERATE", "Filter (VCF)" => "PASS")
    options = VP.RunOptions(input="/tmp/singleton-test.vcf")
    enabled = VP.RunOptions(input=options.input, include_singleton_hets=true)
    config = VP.AppConfig()
    run(rows, opts=options, cfg=config, fam=family) = VP.prioritise_rows(deepcopy(rows), fam, opts, cfg)
    @test isempty(run([row]))
    result = run([row], enabled)
    @test only(result)["candidateCategory"] == "singleton_heterozygous_candidate"
    @test any(section -> section.title == "Heterozygous Variants", VP.report_sections(result, family))
    @test length(run([row], options, VP.AppConfig(thresholds=VP.ThresholdConfig(include_singleton_hets=true)))) == 1
    @test VP.parse_cli(["--include-singleton-hets", "input.vcf", "child"]).include_singleton_hets
    @test !VP.load_config(nothing).thresholds.include_singleton_hets
    mktemp() do path, io
        write(io, "[thresholds]\ninclude_singleton_hets = true\n")
        close(io)
        @test VP.load_config(path).thresholds.include_singleton_hets
    end
    for change in [Dict("GQ (child)" => "1"), Dict("GnomAD_v4_1_AF_all" => "0.1"),
                   Dict("IMPACT" => "LOW"), Dict("GT (child)" => "0/0")]
        @test isempty(run([merge(row, change)], enabled))
    end
    pair = [row, merge(row, Dict("inputPos" => "200"))]
    @test all(r -> r["candidateCategory"] == "singleton_possible_compound_heterozygous_candidate", run(pair, enabled))
    @test length(run(pair)) == 2
    @test only(run([merge(row, Dict("GT (child)" => "1/1"))], enabled))["candidateCategory"] == "recessive_homozygous_candidate"
    @test isempty(run([row], VP.RunOptions(input=options.input, include_singleton_hets=true, recessive_only=true)))
    @test isempty(run([row], VP.RunOptions(input=options.input, include_singleton_hets=true, denovo_only=true)))
    trio = VP.FamilySpec(parent1="mother", parent2="father", affected=["child"])
    trio_row = merge(row, Dict("GT (mother)" => "0/1", "GT (father)" => "0/0"))
    @test run([trio_row], enabled, config, trio) == run([trio_row], options, config, trio)
end

@testset "CSQ cohort frequency and AlphaGenome annotations" begin
    thresholds = VP.ThresholdConfig()
    for config in (VP.AppConfig(), VP.load_config(nothing))
        @test config.thresholds.exeter_genomes_joint_af_cutoff == 0.01
        @test config.thresholds.spliceai_cutoff == 0.05
        @test config.thresholds.alphagenome_splicing_cutoff == 1.0
    end
    mktemp() do path, io
        write(io, "[thresholds]\nexeter_genomes_joint_af_cutoff = 0.02\nspliceai_cutoff = 0.2\nalphagenome_splicing_cutoff = 2.0\n")
        close(io)
        custom = VP.load_config(path).thresholds
        @test custom.exeter_genomes_joint_af_cutoff == 0.02
        @test custom.spliceai_cutoff == 0.2
        @test custom.alphagenome_splicing_cutoff == 2.0
    end
    for category in ("inherited_candidate", "de_novo_candidate", "recessive_homozygous_candidate")
        for (af, expected) in (("0.009", true), ("0.01", false), ("0.02", false), ("", true), (".", true))
            row = Dict{String,Any}("Exeter_Genomes_Joint_AF" => af)
            @test VP.passes_frequency_filter(row, category, 0.001, thresholds) == expected
        end
        @test !VP.passes_frequency_filter(Dict{String,Any}("Exeter_Genomes_Joint_AF" => "0.001", "GnomAD_v4_1_AF_all" => "0.1"), category, 0.001, thresholds)
    end
    for impact in ("MODIFIER", "LOW", "")
        for (spliceai, alphagenome, expected) in (("0.051", "", true), ("0.05", "1", false),
                ("0", "1.01", true), ("", "1.01", true), (".", ".", false), ("0.6", "0", true))
            row = Dict{String,Any}("IMPACT" => impact, "varLocation" => "intron", "Consequence" => "intron_variant",
                "spliceai_max" => spliceai, "AlphaGenome_splicing" => alphagenome)
            @test VP.is_interesting_variant(row, false, thresholds) == expected
        end
    end
    @test !VP.is_interesting_variant(Dict{String,Any}("IMPACT" => "MODIFIER", "varLocation" => "intron", "AlphaGenome_PHRED" => "30"), false, thresholds)

    csq_headers = ["Allele", "Feature_type", "Feature", "SYMBOL", "Consequence", "IMPACT",
        "Exeter_Genomes_Joint_AF", "AlphaGenome_PHRED", "AlphaGenome_splicing"]
    fields = ["1", "100", ".", "A", "G", "50", "PASS",
        "CSQ=G|Transcript|ENST1|GENE1|intron_variant|MODIFIER|0.009|25|1.1", "GT:GQ", "0/1:30"]
    row = only(VP.normalise_vcf_record(fields, ["child"], csq_headers, "GRCh38", "test.vcf"))
    for (field, value) in (("Exeter_Genomes_Joint_AF", "0.009"), ("AlphaGenome_PHRED", "25"), ("AlphaGenome_splicing", "1.1"))
        @test row[field] == value
        @test field in VP.PRIORITY_HEADERS
        @test field in VP.small_variant_table_headers(VP.PRIORITY_HEADERS, VP.FamilySpec(affected=["child"]))
    end
    options = VP.RunOptions(input="test.vcf", include_singleton_hets=true)
    run(r) = VP.prioritise_rows([deepcopy(r)], VP.FamilySpec(affected=["child"]), options, VP.AppConfig())
    @test length(run(row)) == 1
    @test isempty(run(merge(row, Dict("Exeter_Genomes_Joint_AF" => "0.01"))))
    @test isempty(run(merge(row, Dict("AlphaGenome_splicing" => "1"))))
    @test length(run(merge(row, Dict("AlphaGenome_splicing" => "1", "spliceai_max" => "0.051")))) == 1
end
