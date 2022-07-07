
#include("kSP_bigM_method.jl")
include("data_proc.jl")

# An example of computing the k values for big M values in DC OTS problem
#-----------------------------------------------------------------------------------------------------------------
#=
load_factors = [95,100,105]
idx_load_factors = [2]
instances = 1:20
branch_data = DataFrame(CSV.File("Data/118bus_data/Branch_Data_IEEE118_merged.csv"))
#branch_data = DataFrame(CSV.File("Data/300bus_data/IEEE300_branch_merged.csv"))
branch_index = Dict{Tuple{Int64,Int64},Int64}() # a dictionary that mapps buses into branch index
for i in branch_data.index
    fbus = branch_data.node1[i]
    tbus = branch_data.node2[i]
    branch_index[fbus,tbus] = i
    branch_index[tbus,fbus] = i
end
gen_data = DataFrame(CSV.File("Data/118bus_data/Generator_Data_IEEE118.csv"))
#gen_data = DataFrame(CSV.File("Data/300bus_data/IEEE300_gen.csv"))
ksp_paths = DataFrame(CSV.File("results_118/ksp_paths.csv"))
max_pathS = [25]
max_edgeS = [3,5,7]
compute_ksp_k_val(1,false,load_factors,idx_load_factors,instances,branch_data,gen_data,ksp_paths,max_pathS,max_edgeS)
=#
#------------------------------------------------------------------------------------------------------------------



branch_data =  DataFrame(CSV.File("Data/300bus_data/IEEE300_branch_merged.csv"))
gen_data = DataFrame(CSV.File("Data/300bus_data/IEEE300_gen.csv"))
ksp_weights = DataFrame(CSV.File("results_300/ksp_weights.csv"))

load_factors = [90,95,100,105,110]
idx_load_factors = [3]
instances = 1:20
card_limits = [45]
max_pathS = [15,20,25,30]
max_edgeS = [5,7,9]
conv_levelS = [5]
output_file = "results_300/OTS_ksp_para_tuning_new.csv"
output_df = OTS_ksp_bigM(2, load_factors, idx_load_factors, instances, card_limits, branch_data, gen_data, ksp_weights, max_pathS, max_edgeS, conv_levelS)
CSV.write(output_file, output_df)


# analyze optimality gaps
#=
input_df = DataFrame(CSV.File("results_118/OTS_ksp_card.csv"))
dropmissing!(input_df)
num_instances = 20
card_limits = [10,15,20,25]
output_file = "results_118/analyze_opt_gap_ksp_118.csv"
analyze_opt_gap(card_limits, num_instances, input_df, output_file)
=#

# Fuller heuristic
#=
fuller_input = DataFrame(CSV.File("results_300/OTS_fuller_heuristic_300.csv"))
fuller_input = fuller_input[fuller_input.size_candi .== 2,:]
lwp_input = DataFrame(CSV.File("results_300/OTS_lwp_paraTune_300.csv"))
load_factors = [0.95,1.0,1.05]
lwp_load_factors = [95,100,105]
output_df = analyze_fuller(2, fuller_input, lwp_input, load_factors, lwp_load_factors, 20)
output_file = "results_300/fuller_compare_ksp_300.csv"
CSV.write(output_file, output_df)
=#