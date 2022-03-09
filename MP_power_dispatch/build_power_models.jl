using JuMP

function compute_line_susceptance_dc(data)
    b = Dict{Int64,Float64}()
    for (idx,n1,n2) in data.branches
        b[idx] = -data.branch_react[idx]/(data.branch_react[idx]^2 + data.branch_resist[idx]^2)
    end
    return b
end

function compute_branch_admittance(data::power_system_data)
    g = Dict{Tuple{Int64,Int64,Int64}, Float64}()
    b = Dict{Tuple{Int64,Int64,Int64}, Float64}()
    G = Dict{Tuple{Int64,Int64,Int64}, Float64}()
    B = Dict{Tuple{Int64,Int64,Int64}, Float64}()

    for (i,n1,n2) in data.branches
        g[i,n1,n2] = data.branch_resist[i] / (data.branch_resist[i]^2 + data.branch_react[i]^2)
        G[i,n1,n2] = -g[i,n1,n2]/data.tap_ratio[i]
        G[i,n2,n1] = G[i,n1,n2]
        b[i,n1,n2] = -data.branch_react[i] / (data.branch_resist[i]^2 + data.branch_react[i]^2)
        B[i,n1,n2] = -b[i,n1,n2]/data.tap_ratio[i]
        B[i,n2,n1] = B[i,n1,n2]
    end
    return G,B
end

# Add variables to the model
#--------------------------------------------------------------------------------------------
function add_var_dc(model, data)
    @variable(model, real_prod[data.gen], lower_bound = 0)
    @variable(model, volt_ang[data.buses])
    @variable(model, flow[data.branches])
    return model
end

function add_var_ac(model, data, asym_branches)
    @variable(model, real_prod[data.gen])
    @variable(model, imag_prod[data.gen])
    @variable(model, volt_ang[data.buses], start = 0)
    @variable(model, volt_mag[data.buses], start = 1)
    @variable(model, real_flow[asym_branches])
    @variable(model, imag_flow[asym_branches])
    @variable(model, volt_ang_diff[asym_branches])
    return model
end

function add_var_switch(model, data)
    @variable(model, branch_status[data.branches], binary = true)
    return model
end
#------------------------------------------------------------------------------------------------------------

# Add constraints to the model
#------------------------------------------------------------------------------------------------------------
function add_constr_dc_nodal_flow_balance(model, data)
    flow = model[:flow]
    production = model[:real_prod]
    for b in data.buses 
        if b in data.gen_bus
            @constraint(model, sum(flow[(idx,n1,n2)] for (idx,n1,n2) in data.branches if n2 == b)
            - sum(flow[(idx,n1,n2)] for (idx,n1,n2) in data.branches if n1 == b)
            + sum(production[g] for g in data.gen if data.gen_bus[g] == b)
            == data.bus_real_loads[b]
            )
        else
            @constraint(model, sum(flow[(idx,n1,n2)] for (idx,n1,n2) in data.branches if n2 == b)
            - sum(flow[(idx,n1,n2)] for (idx,n1,n2) in data.branches if n1 == b)
            == data.bus_real_loads[b]
            )
        end
    end
    return model
end


function add_constr_dc_KVL(model,data, b)
    flow = model[:flow]
    volt_ang = model[:volt_ang]
    @constraint(model, [(idx,n1,n2) in data.branches], flow[(idx,n1,n2)] + data.baseMVA*b[idx]*(volt_ang[n1] - volt_ang[n2]) == 0)
    return model
end

function add_constr_dc_KVL_switch(model,data,b,bigM)
    flow = model[:flow]
    volt_ang = model[:volt_ang]
    branch_status = model[:branch_status]
    for (idx, n1, n2) in data.branches
        if bigM[idx] > 0
            @constraint(model, flow[(idx,n1,n2)] + data.baseMVA*b[idx]*(volt_ang[n1] - volt_ang[n2]) <= bigM[idx]*(1-branch_status[(idx,n1,n2)]) ) # real power flow across each branch
        else 
            @constraint(model, flow[(idx,n1,n2)] <= 0)
        end
    end
    @constraint(model, [(idx,n1,n2) in data.branches], flow[(idx,n1,n2)] + data.baseMVA*b[idx]*(volt_ang[n1] - volt_ang[n2]) >= -bigM[idx]*(1-branch_status[(idx,n1,n2)]) )
    
    return model
end

function add_constr_dc_line_flow(model,data,b,ang_diff_limit)
    flow = model[:flow]
    for (idx,n1,n2) in data.branches
        if data.branch_current_limit[idx] > 0.01
            @constraint(model, flow[(idx,n1,n2)] <= data.branch_current_limit[idx])
            @constraint(model, flow[(idx,n1,n2)] >= -data.branch_current_limit[idx])
        end
    end
    @constraint(model, [(idx,n1,n2) in data.branches], flow[(idx,n1,n2)] <= -data.baseMVA*b[idx]*ang_diff_limit)
    @constraint(model, [(idx,n1,n2) in data.branches], flow[(idx,n1,n2)] >= data.baseMVA*b[idx]*ang_diff_limit)
    return model
end

function add_constr_dc_line_flow_switch(model,data,b,ang_diff_limit)
    flow = model[:flow]
    branch_status = model[:branch_status]
    for (idx,n1,n2) in data.branches
        if data.branch_current_limit[idx] > 0.01
            @constraint(model, flow[(idx,n1,n2)] <= branch_status[(idx,n1,n2)] * data.branch_current_limit[idx])
            @constraint(model, flow[(idx,n1,n2)] >= -branch_status[(idx,n1,n2)] * data.branch_current_limit[idx])
        end
    end
    @constraint(model, [(idx,n1,n2) in data.branches], flow[(idx,n1,n2)] <= -data.baseMVA*b[idx]*ang_diff_limit*branch_status[(idx,n1,n2)])
    @constraint(model, [(idx,n1,n2) in data.branches], flow[(idx,n1,n2)] >= data.baseMVA*b[idx]*ang_diff_limit*branch_status[(idx,n1,n2)])
    return model
end

function add_constr_real_gen_limit(model, data)
    real_prod = model[:real_prod]
    @constraint(model, [g in data.gen], real_prod[g] <= data.gen_pmax[g]/data.baseMVA)
    @constraint(model, [g in data.gen], real_prod[g] >= data.gen_pmin[g]/data.baseMVA)
    return model
end

function add_constr_imag_gen_limit(model, data)
    imag_prod = model[:imag_prod]
    @constraint(model, [g in data.gen], imag_prod[g] <= data.gen_qmax[g]/data.baseMVA) 
    @constraint(model, [g in data.gen], imag_prod[g] >= data.gen_qmin[g]/data.baseMVA)
    return model
end

function is_end_node(node, edge)
    if node == edge[2]
        return true
    elseif node == edge[3]
        return true
    else
        return false
    end
end

function add_constr_ac_nodal_flow_balance(model, data, asym_branches)
    real_prod = model[:real_prod]
    imag_prod = model[:imag_prod]
    volt_mag = model[:volt_mag]
    real_flow = model[:real_flow]
    imag_flow = model[:imag_flow]
    
    adjacent_nodes = Dict{Int64, Vector{Int64}}()
    for b in data.buses
        edge_cover = findall(is_end_node.(b, data.branches))
        adjacent_node(node, edge) = node == edge[2] ? edge[3] : edge[2]
        adjacent_nodes[b] = [adjacent_node(b, edge) for edge in data.branches[edge_cover]]
    end
    
    @NLconstraint(model,[b in data.buses], data.bus_shunt_conduc[b] * volt_mag[b]^2 
                    + sum(real_flow[(i,n1,n2)] for (i,n1,n2) in asym_branches if (n1 == b && n2 in adjacent_nodes[b]))
                    == sum(real_prod[g] for g in data.gen if data.gen_bus[g] == b) - data.bus_real_loads[b]/data.baseMVA)
    @NLconstraint(model,[b in data.buses], -data.bus_shunt_suscep[b] * volt_mag[b]^2 
                    + sum(imag_flow[(i,n1,n2)] for (i,n1,n2) in asym_branches if (n1 == b && n2 in adjacent_nodes[b]))
                    == sum(imag_prod[g] for g in data.gen if data.gen_bus[g] == b) - data.bus_imag_loads[b]/data.baseMVA) 
    return model
end

function add_constr_ac_KVL_polar(model, data, asym_branches, G, B)
    volt_ang_diff = model[:volt_ang_diff]
    volt_mag = model[:volt_mag]
    real_flow = model[:real_flow]
    imag_flow = model[:imag_flow]
    @NLconstraint(model, [(i,n1,n2) in asym_branches], real_flow[(i,n1,n2)] == -G[i,n1,n2]*volt_mag[n1]^2 
                    + G[i,n1,n2]*volt_mag[n1]*volt_mag[n2]*cos(volt_ang_diff[(i,n1,n2)])
                    + B[i,n1,n2]*volt_mag[n1]*volt_mag[n2]*sin(volt_ang_diff[(i,n1,n2)]) )
    @NLconstraint(model, [(i,n1,n2) in asym_branches], imag_flow[(i,n1,n2)] == 
                    -(-B[i,n1,n2] + 0.5*data.branch_shunt_suscep[i])*volt_mag[n1]^2
                    - B[i,n1,n2]*volt_mag[n1]*volt_mag[n2]*cos(volt_ang_diff[(i,n1,n2)])
                    + G[i,n1,n2]*volt_mag[n1]*volt_mag[n2]*sin(volt_ang_diff[(i,n1,n2)]) )
    return model
end

function add_constr_ac_line_flow(model, data, asym_branches)
    real_flow = model[:real_flow]
    imag_flow = model[:imag_flow]
    for (i,n1,n2) in data.branches
        if data.branch_current_limit[i] > 0.01
            @NLconstraint(model, real_flow[(i,n1,n2)]^2 + imag_flow[(i,n1,n2)]^2 <= (data.branch_current_limit[i]/data.baseMVA)^2)
            @NLconstraint(model, real_flow[(i,n2,n1)]^2 + imag_flow[(i,n2,n1)]^2 <= (data.branch_current_limit[i]/data.baseMVA)^2)
        else
            @NLconstraint(model, real_flow[(i,n1,n2)]^2 + imag_flow[(i,n1,n2)]^2 <= 50^2)
            @NLconstraint(model, real_flow[(i,n2,n1)]^2 + imag_flow[(i,n2,n1)]^2 <= 50^2)
        end
    end
    return model
end

function add_constr_ac_volt_limit_polar(model, data, asym_branches, max_ang_diff)
    volt_ang = model[:volt_ang]
    volt_ang_diff = model[:volt_ang_diff]
    volt_mag = model[:volt_mag]
    @constraint(model, [(i,n1,n2) in asym_branches], volt_ang_diff[(i,n1,n2)] == volt_ang[n1] - volt_ang[n2])
    @constraint(model, [(i,n1,n2) in asym_branches], volt_ang_diff[(i,n1,n2)] <= max_ang_diff)
    @constraint(model, [(i,n1,n2) in asym_branches], volt_ang_diff[(i,n1,n2)] >= -max_ang_diff)
    @constraint(model, [b in data.buses], volt_mag[b] <= data.bus_volt_mag_max[b])
    @constraint(model, [b in data.buses], volt_mag[b] >= data.bus_volt_mag_min[b])
    @constraint(model, volt_ang[4] == 0) # fix voltage angle at reference bus
    return model
end
#----------------------------------------------------------------------------------------------------------------------

# Add objective function to the model
#--------------------------------------------------------------------------------------------------------------------
function add_obj_linear_prod_cost(model,data)
    production = model[:real_prod]
    @objective(model, Min, sum(production[g]*data.gen_c1[g] for g in data.gen))
    return model
end

function add_obj_quad_prod_cost(model,data)
    prod = model[:real_prod]
    @objective(model, Min, sum(prod[g]*data.gen_c1[g] + prod[g]^2*data.gen_c2[g] for g in data.gen))
    return model
end
#------------------------------------------------------------------------------------------------------------------


function build_DCOPF(data::power_system_data, ang_diff_limit)
    b = compute_line_susceptance_dc(data)
    model_dcopf = Model();
    model_dcopf = add_var_dc(model_dcopf, data)
    model_dcopf = add_constr_real_gen_limit(model_dcopf, data)
    model_dcopf = add_constr_dc_nodal_flow_balance(model_dcopf, data)
    model_dcopf = add_constr_dc_KVL(model_dcopf,data,b)
    model_dcopf = add_constr_dc_line_flow(model_dcopf,data,b,ang_diff_limit)
    model_dcopf = add_obj_linear_prod_cost(model_dcopf,data)

    return model_dcopf
end

function build_DCOTS(data::power_system_data, ang_diff_limit, bigM)
    b = compute_line_susceptance_dc(data)
    model_dcots = Model();
    model_dcots = add_var_dc(model_dcots, data)
    model_dcots = add_var_switch(model_dcots, data)
    model_dcots = add_constr_real_gen_limit(model_dcots, data)
    model_dcots = add_constr_dc_nodal_flow_balance(model_dcots, data)
    model_dcots = add_constr_dc_KVL_switch(model_dcots,data,b,bigM)
    model_dcots = add_constr_dc_line_flow_switch(model_dcots,data,b,ang_diff_limit)
    model_dcots = add_obj_linear_prod_cost(model_dcots,data)

    return model_dcots
end

function build_ACOPF_polar(data::power_system_data, max_ang_diff)
    asym_branches = Vector{Tuple{Int64,Int64,Int64}}()
    for (i,n1,n2) in data.branches
        push!(asym_branches, (i,n1,n2))
        push!(asym_branches, (i,n2,n1))
    end
    G,B = compute_branch_admittance(data)
    acopf_polar = Model();
    acopf_polar = add_var_ac(acopf_polar, data, asym_branches)
    acopf_polar = add_constr_real_gen_limit(acopf_polar, data)
    acopf_polar = add_constr_imag_gen_limit(acopf_polar, data)
    acopf_polar = add_constr_ac_nodal_flow_balance(acopf_polar, data, asym_branches)
    acopf_polar = add_constr_ac_KVL_polar(acopf_polar, data, asym_branches, G, B)
    acopf_polar = add_constr_ac_line_flow(acopf_polar, data, asym_branches)
    acopf_polar = add_constr_ac_volt_limit_polar(acopf_polar, data, asym_branches, max_ang_diff)
    acopf_polar = add_obj_quad_prod_cost(acopf_polar,data)
    
    return acopf_polar
end