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
