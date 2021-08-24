using CSV
using Dates
using DataFrames
using Printf

using MAT
using MATLAB

using JLD

using Parameters

using DelimitedFiles

using Random
using StatsBase
using PyPlot
using Distributions

using DifferentialEquations

using LSODA

using SparseArrays
using LinearAlgebra

using Roots
using QuadGK

using EMpht

import Base.tryparse

include("funcs/structs.jl")
include("funcs/common.jl")
include("funcs/fluid.jl")
include("funcs/import_line.jl")
include("funcs/generate_sys.jl")
include("funcs/plotting.jl")

