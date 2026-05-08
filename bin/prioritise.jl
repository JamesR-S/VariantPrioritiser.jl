#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "VariantPrioritiser.jl"))

using .VariantPrioritiser

VariantPrioritiser.main()
