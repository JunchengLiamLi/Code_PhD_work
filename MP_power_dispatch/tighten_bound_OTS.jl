using Gurobi

include("build_power_models.jl")

function build_max_ang_diff_DCOTS(data::power_system_data, bigM_plus, bigM_minus, max_ang_diff, cost_ub)
    b = compute_line_susceptance_dc(data)
    model_mad = Model();
    model_mad = add_var_dc(model_mad, data)
    model_mad = add_var_switch(model_mad, data)
    model_mad = add_constr_real_gen_limit(model_mad, data)
    model_mad = add_constr_dc_nodal_flow_balance(model_mad, data)
    model_mad = add_constr_dc_KVL_switch(model_mad,data,b,bigM_plus,bigM_minus)
    model_mad = add_constr_dc_line_flow_switch(model_mad,data,b,max_ang_diff)
    production = model_mad[:real_prod]
    @constraint(model_mad, 100*sum(production[g]*data.gen_c1[g] for g in data.gen) <= cost_ub)

    return model_mad
end

function compute_max_ang_diff_ub(model_mad, data::power_system_data, idx::Int64)
    volt_ang = model_mad[:volt_ang]
    n1 = data.branches[idx][2]
    n2 = data.branches[idx][3]
    @objective(model_mad, Max, volt_ang[n1]- volt_ang[n2])
    
    set_optimizer(model_mad, Gurobi.Optimizer)
    set_silent(model_mad)
    set_optimizer_attribute(model_mad, "NodeLimit", 0)
    optimize!(model_mad)

    return objective_bound(model_mad)
end


function compute_min_ang_diff_lb(model_mad, data::power_system_data, idx::Int64)
    volt_ang = model_mad[:volt_ang]
    n1 = data.branches[idx][2]
    n2 = data.branches[idx][3]
    @objective(model_mad, Min, volt_ang[n1]- volt_ang[n2])
    
    set_optimizer(model_mad, Gurobi.Optimizer)
    set_silent(model_mad)
    set_optimizer_attribute(model_mad, "NodeLimit", 0)
    optimize!(model_mad)

    return objective_bound(model_mad)
end

function strengthen_bigM_DCOTS(data::power_system_data, init_bigM, max_ang_diff; stop_tol=0.05, max_iter = 20)
    mdcopf = build_DCOPF(data, max_ang_diff)
    set_optimizer(mdcopf, Gurobi.Optimizer)
    set_silent(mdcopf)
    optimize!(mdcopf)
    if has_values(mdcopf)
        cost_ub = objective_value(mdcopf)
    else
        cost_ub = 100000000
    end
    
    num_branches = length(data.branches)
    streng_bigM_plus = zeros(num_branches) 
    streng_bigM_minus = zeros(num_branches)
    bigM_plus = copy(init_bigM)
    bigM_minus = -init_bigM
    
    b = compute_line_susceptance_dc(data)
    for i in 1:num_branches
        iter = 0
        while max(abs((bigM_plus[i]-streng_bigM_plus[i])/ bigM_plus[i]), abs((bigM_minus[i]-streng_bigM_minus[i])/ bigM_minus[i])) > stop_tol && iter < max_iter
            if iter > 0
                bigM_plus[i] = streng_bigM_plus[i]
                bigM_minus[i] = streng_bigM_minus[i]
            end
            mdcots_mad = build_max_ang_diff_DCOTS(data, bigM_plus, bigM_minus, max_ang_diff, cost_ub)
            mad_ub = compute_max_ang_diff_ub(mdcots_mad, data, i)
            mad_lb = compute_min_ang_diff_lb(mdcots_mad, data, i)
            streng_bigM_plus[i] = -mad_ub*b[i]
            streng_bigM_minus[i] = -mad_lb*b[i]
            iter += 1
        end
    end
    return streng_bigM_plus, streng_bigM_minus
end

function branch_indep_strengthen_bigM_DCOTS(data::power_system_data, init_bigM, max_ang_diff; stop_tol=0.02, max_iter = 20)
    mdcopf = build_DCOPF(data, max_ang_diff)
    set_optimizer(mdcopf, Gurobi.Optimizer)
    set_silent(mdcopf)
    optimize!(mdcopf)
    if has_values(mdcopf)
        cost_ub = objective_value(mdcopf)
    else
        cost_ub = 100000000
    end
    
    num_branches = length(data.branches)
    streng_bigM_plus = zeros(num_branches) 
    streng_bigM_minus = zeros(num_branches)
    
    b = compute_line_susceptance_dc(data)
    for i in 1:num_branches
        bigM_plus = copy(init_bigM)
        bigM_minus = -init_bigM
        iter = 0
        while max(abs((bigM_plus[i]-streng_bigM_plus[i])/ bigM_plus[i]), abs((bigM_minus[i]-streng_bigM_minus[i])/ bigM_minus[i])) > stop_tol && iter < max_iter
            if iter > 0
                bigM_plus[i] = streng_bigM_plus[i]
                bigM_minus[i] = streng_bigM_minus[i]
            end
            mdcots_mad = build_max_ang_diff_DCOTS(data, bigM_plus, bigM_minus, max_ang_diff, cost_ub)
            br_status = mdcots_mad[:branch_status]
            @constraint(mdcots_mad, br_status[data.branches[i]] == 0)
            
            mad_ub = compute_max_ang_diff_ub(mdcots_mad, data, i)
            mad_lb = compute_min_ang_diff_lb(mdcots_mad, data, i)
            streng_bigM_plus[i] = -mad_ub*b[i]
            streng_bigM_minus[i] = -mad_lb*b[i]
            iter += 1
        end
    end
    return streng_bigM_plus, streng_bigM_minus
end