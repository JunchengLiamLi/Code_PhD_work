# compute longest paths in IEEE118(Blumsack) network 
#= ---------------------------------------------------------------------------------------------------------
include("format_data.jl")
include("longest_path.jl")

baseMVA = 100
bus_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_bus.csv"))
branch_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_branch.csv"))
gen_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_gen.csv"))
gencost_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_gencost.csv"))

case118 = input_matpower_data(bus_data_case118, branch_data_case118, gen_data_case118, gencost_data_case118, baseMVA)

num_branch = 186
long_paths_case118 = Vector{Float64}(undef, num_branch)
for (i, n1 ,n2) in case118.branches
    path_length, path_length_ub = MTZ_longest_path(case118, n1, n2, 600, 0.005)
    long_paths_case118[i] = path_length_ub
end
lp_118_df = DataFrame(longest_path_length = long_paths_case118)
CSV.write("MP_power_dispatch/OBBT_DC_OTS/Blumsack_118_longest_paths.csv", lp_118_df)
------------------------------------------------------------------------------------------------------------------- =# 


# run optimality-based bound tightening method 
#= ----------------------------------------------------------------------------------------------------------------
include("tighten_bound_OTS.jl")

function run_OBBT_DCOTS(network, max_ang_diff)
    @assert network == "118" || network == "300" "network must be 118-bus case or 300-bus case"
    if network == "118"
        baseMVA = 100
        bus_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_bus.csv"))
        branch_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_branch.csv"))
        gen_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_gen.csv"))
        gencost_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_gencost.csv"))

        case118 = input_matpower_data(bus_data_case118, branch_data_case118, gen_data_case118, gencost_data_case118, baseMVA)
        
        long_paths_case118_lengths = DataFrame(CSV.File("MP_power_dispatch/OBBT_DC_OTS/Blumsack_118_longest_paths.csv"))[:,1]
        case118_bigM = similar(long_paths_case118_lengths)
        case118_b = compute_line_susceptance_dc(case118)
        for i in 1:length(case118_bigM)
            case118_bigM[i] = abs(long_paths_case118_lengths[i] * max_ang_diff * case118_b[i])
        end
        num_demand = 30
        num_branch = 186
        case118_random_demands = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_random_demands.csv"))
        sBigM_plus_mx = Matrix{Float64}(undef, num_branch, num_demand)
        sBigM_minus_mx = Matrix{Float64}(undef, num_branch, num_demand)
        for i in 1:num_demand
            for j in 1:118
                case118.bus_real_loads[j] = case118_random_demands[i,j]
            end  
            streng_bigM_plus, streng_bigM_minus = branch_indep_strengthen_bigM_DCOTS(case118, case118_bigM, max_ang_diff)
            sBigM_plus_mx[:,i] = streng_bigM_plus
            sBigM_minus_mx[:,i] = streng_bigM_minus
        end
    end
    return sBigM_plus_mx, sBigM_minus_mx
end

max_ang_diff = pi/6
sBigM_plus_mx, sBigM_minus_mx = run_OBBT_DCOTS("118", max_ang_diff)

sBigM_plus_df = DataFrame(sBigM_plus_mx,:auto)
sBigM_minus_df = DataFrame(sBigM_minus_mx,:auto)
CSV.write("MP_power_dispatch/OBBT_DC_OTS/BI_sBigM_plus_118_random_demands.csv", sBigM_plus_df)
CSV.write("MP_power_dispatch/OBBT_DC_OTS/BI_sBigM_minus_118_random_demands.csv", sBigM_minus_df)
# ----------------------------------------------------------------------------------------------------------------- =#

# run DC OTS 
#= ---------------------------------------------------------------------------------------------------------------
include("run_power_models.jl")

baseMVA = 100
bus_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_bus.csv"))
branch_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_branch.csv"))
gen_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_gen.csv"))
gencost_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_gencost.csv"))

case118 = input_matpower_data(bus_data_case118, branch_data_case118, gen_data_case118, gencost_data_case118, baseMVA)

max_ang_diff = pi/6

long_paths_case118_lengths = DataFrame(CSV.File("MP_power_dispatch/OBBT_DC_OTS/Blumsack_118_longest_paths.csv"))[:,1]
case118_bigM = similar(long_paths_case118_lengths)
case118_b = compute_line_susceptance_dc(case118)
for i in 1:length(case118_bigM)
    case118_bigM[i] = abs(long_paths_case118_lengths[i] * max_ang_diff * case118_b[i])
end
case118_sBigM_plus_df = DataFrame(CSV.File("MP_power_dispatch/OBBT_DC_OTS/BI_sBigM_plus_118_random_demands.csv"))
case118_sBigM_minus_df = DataFrame(CSV.File("MP_power_dispatch/OBBT_DC_OTS/BI_sBigM_minus_118_random_demands.csv"))

num_demand = 30
num_branch = 186

case118_random_demands = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_random_demands.csv"))
#=
obj_val_lpBigM_vt = Vector{Float64}()
gap_lpBigM_vt = Vector{Float64}()
time_lpBigM_vt = Vector{Float64}()
nodes_lpBigM_vt = Vector{Int64}()
term_lpBigM_vt = Vector()
=#

obj_val_sBigM_vt = Vector{Float64}()
gap_sBigM_vt = Vector{Float64}()
time_sBigM_vt = Vector{Float64}()
nodes_sBigM_vt = Vector{Int64}()
term_sBigM_vt = Vector()
ub_vt = Vector()

for i in 1:num_demand
    for j in 1:118
        case118.bus_real_loads[j] = case118_random_demands[i,j]
    end
    
    #mdcots_lpBigM = build_DCOTS(case118, max_ang_diff, case118_bigM, -case118_bigM)
    #obj_val_lpBigM, gap_lpBigM, time_lpBigM, nodes_lpBigM, term_lpBigM = run_DC_OTS(case118, mdcots_lpBigM, Gurobi.Optimizer, 15, 3600)
    #push!(obj_val_lpBigM_vt, obj_val_lpBigM)
    #push!(gap_lpBigM_vt, gap_lpBigM)
    #push!(time_lpBigM_vt, time_lpBigM)
    #push!(nodes_lpBigM_vt, nodes_lpBigM)
    #push!(term_lpBigM_vt, term_lpBigM)

    case118_sBigM_plus = case118_sBigM_plus_df[:,i]
    case118_sBigM_minus = case118_sBigM_minus_df[:,i]
    mdcots_sBigM = build_DCOTS(case118, max_ang_diff, case118_sBigM_plus, case118_sBigM_minus)
    obj_val_sBigM, gap_sBigM, time_sBigM, nodes_sBigM, term_sBigM = run_DC_OTS(case118, mdcots_sBigM, Gurobi.Optimizer, 15, 3600)
    push!(obj_val_sBigM_vt, obj_val_sBigM)
    push!(gap_sBigM_vt, gap_sBigM)
    push!(time_sBigM_vt, time_sBigM)
    push!(nodes_sBigM_vt, nodes_sBigM)
    push!(term_sBigM_vt, term_sBigM)

    mdcopf = build_DCOPF(case118, max_ang_diff)
    set_optimizer(mdcopf, Gurobi.Optimizer)
    optimize!(mdcopf)
    push!(ub_vt, has_values(mdcopf))
end
#DCOTS_lpBigM_df = DataFrame(obj_val = obj_val_lpBigM_vt, gap = gap_lpBigM_vt, time = time_lpBigM_vt, nodes = nodes_lpBigM_vt, status = term_lpBigM_vt)
DCOTS_sBigM_df = DataFrame(obj_val = obj_val_sBigM_vt, gap = gap_sBigM_vt, time = time_sBigM_vt, nodes = nodes_sBigM_vt, status = term_sBigM_vt, ub = ub_vt)
#CSV.write("MP_power_dispatch/OBBT_DC_OTS/DCOTS_118_lpBigM.csv", DCOTS_lpBigM_df)
CSV.write("MP_power_dispatch/OBBT_DC_OTS/DCOTS_118_sBigM_BI.csv", DCOTS_sBigM_df)
#------------------------------------------------------------------------------------------------------------------------ =#


include("run_power_models.jl")

baseMVA = 100
bus_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_bus.csv"))
branch_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_branch.csv"))
gen_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_gen.csv"))
gencost_data_case118 = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_gencost.csv"))

case118 = input_matpower_data(bus_data_case118, branch_data_case118, gen_data_case118, gencost_data_case118, baseMVA)

max_ang_diff = pi/6

long_paths_case118_lengths = DataFrame(CSV.File("MP_power_dispatch/OBBT_DC_OTS/Blumsack_118_longest_paths.csv"))[:,1]
case118_bigM = similar(long_paths_case118_lengths)
case118_b = compute_line_susceptance_dc(case118)
for i in 1:length(case118_bigM)
    case118_bigM[i] = 100 * abs(long_paths_case118_lengths[i] * max_ang_diff * case118_b[i])
end

load_index = 24

case118_random_demands = DataFrame(CSV.File("MP_power_dispatch/DC_Data/Blumsack118_random_demands.csv"))
for j in 1:118
    case118.bus_real_loads[j] = case118_random_demands[load_index,j]
end

case118_dcots_lpBigM = build_DCOTS(case118, max_ang_diff, case118_bigM, -case118_bigM)
result_lpBigM = run_DC_OTS(case118, case118_dcots_lpBigM, Gurobi.Optimizer, 15, 120)

case118_sBigM_plus_df = DataFrame(CSV.File("MP_power_dispatch/OBBT_DC_OTS/BI_sBigM_plus_118_random_demands.csv"))
case118_sBigM_minus_df = DataFrame(CSV.File("MP_power_dispatch/OBBT_DC_OTS/BI_sBigM_minus_118_random_demands.csv"))

case118_sBigM_plus = case118_sBigM_plus_df[:,load_index]
case118_sBigM_minus = case118_sBigM_minus_df[:,load_index]
mdcots_sBigM = build_DCOTS(case118, max_ang_diff, case118_sBigM_plus, case118_sBigM_minus)
result_sBigM = run_DC_OTS(case118, mdcots_sBigM, Gurobi.Optimizer, 15, 120)
