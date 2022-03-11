using Gurobi

include("build_power_models.jl")

function run_DC_OTS(data::power_system_data, model, solver, limit_switch_off::Int64, time_limit)
    br_status = model[:branch_status]
    num_branch = length(data.branches)
    @constraint(model, sum(br_status[(i,fbus,tbus)] for (i,fbus,tbus) in data.branches) >= num_branch - limit_switch_off)
    set_optimizer(model, solver)
    set_optimizer_attribute(model, "TimeLimit", time_limit)
    optimize!(model)
    
    return objective_value(model), relative_gap(model), solve_time(model), node_count(model), termination_status(model)
    # return value.(model[:branch_status]), objective_value(model), relative_gap(model), solve_time(model), node_count(model), termination_status(model)
end