include("format_data.jl")
include("longest_path.jl")

baseMVA = 100
bus_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_bus.csv"))
branch_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_branch.csv"))
gen_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_gen.csv"))
gencost_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_gencost.csv"))

case118 = input_matpower_data(bus_data_case118, branch_data_case118, gen_data_case118, gencost_data_case118, baseMVA)

num_branch_considered = 70
long_paths_case118 = Vector{Float64}(undef, num_branch_considered)
for (i, n1 ,n2) in case118.branches[186-num_branch_considered+1:186]
    path_length, path_length_ub = MTZ_longest_path(case118, n1, n2, 600, 0.005)
    long_paths_case118[i-(186-num_branch_considered)] = path_length_ub
end
lp_118_df = DataFrame(longest_path_length = long_paths_case118)
CSV.write("MP_power_dispatch/OBBT_DC_OTS/Blumsack_118_longest_pathsA.csv", lp_118_df)
