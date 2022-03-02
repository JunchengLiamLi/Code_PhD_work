using CSV,DataFrames

struct power_system_data
    buses::Vector{Int64}
    bus_shunt_suscep::Vector{Float64}
    bus_shunt_conduc::Vector{Float64}
    bus_real_loads::Vector{Float64}
    bus_imag_loads::Vector{Float64}
    bus_volt_mag_max::Vector{Float64}
    bus_volt_mag_min::Vector{Float64}

    branches::Vector{Tuple{Int64,Int64,Int64}}
    tap_ratio::Dict{Int64,Float64}
    branch_shunt_suscep::Dict{Int64,Float64}
    branch_shunt_conduc::Dict{Int64,Float64}
    phase_angle::Dict{Int64,Float64}
    branch_resist::Dict{Int64,Float64}
    branch_react::Dict{Int64,Float64}
    branch_current_limit::Dict{Int64,Float64}

    gen::Vector{Int64}
    gen_bus::Vector{Int64}
    gen_pmax::Vector{Float64}
    gen_pmin::Vector{Float64}
    gen_qmax::Vector{Float64}
    gen_qmin::Vector{Float64}
    gen_real_cost::Vector{Float64}

    baseMVA::Float64
end

function compute_branch_admittance(data::power_system_data)
    # compute branch series admittance
    #------------------------------------------------------------------------------------------------
    y_tilde = Dict{Tuple{Int64,Int64},Complex{Float64}}()
    a = Dict{Tuple{Int64,Int64},Complex{Float64}}()
    for (idx,fnode,tnode) in data.branches
        y_tilde[fnode,tnode] = data.branch_resist[idx]/(data.branch_resist[idx]^2 + data.branch_react[idx]^2) - im*data.branch_react[idx]/(data.branch_resist[idx]^2 + data.branch_react[idx]^2)
        y_tilde[tnode,fnode] = y_tilde[fnode, tnode]
        a[fnode,tnode] = data.tap_ratio[idx]*(cos(data.phase_angle[idx]) + im * sin(data.phase_angle[idx]))
        a[tnode,fnode] = data.tap_ratio[idx]*(cos(-data.phase_angle[idx]) + im * sin(-data.phase_angle[idx]))
    end
    #------------------------------------------------------------------------------------------------
    
    num_bus = length(data.buses)
    Y = Matrix{Complex{Float64}}(undef, num_bus, num_bus)
    for (idx, fnode, tnode) in data.branches
        Y[fnode, tnode] = 0.0 + 0.0*im
        if (idx, fnode, tnode) in data.branches
            Y[fnode, tnode] += -sum(y_tilde[n1,n2]/conj(a[n1,n2]) for (i,n1,n2) in data.branches if n1 == fnode && n2 == tnode)
        end 
        if (idx, tnode, fnode) in data.branches
            Y[fnode, tnode] += -sum(y_tilde[n1,n2]/a[n1,n2] for (i,n1,n2) in data.branches if n2 == fnode && n1 == tnode)
        end
    end
    for (idx, tnode, fnode) in data.branches
        Y[fnode, tnode] = 0.0 + 0.0*im
        if (idx, fnode, tnode) in data.branches
            Y[fnode, tnode] += -sum(y_tilde[n1,n2]/conj(a[n1,n2]) for (i,n1,n2) in data.branches if n1 == fnode && n2 == tnode)
        end 
        if (idx, tnode, fnode) in data.branches
            Y[fnode, tnode] += -sum(y_tilde[n1,n2]/a[n1,n2] for (i,n1, n2) in data.branches if n2 == fnode && n1 == tnode)
        end
    end
    num_branch = length(data.branches)
    from_buses = Vector{Int64}(undef, num_branch)
    to_buses = Vector{Int64}(undef, num_branch)
    for i in 1:num_branch
        from_buses[i] = data.branches[i][2]
        to_buses[i] = data.branches[i][3]
    end

    for b in data.buses
        Y[b,b] = data.bus_shunt_conduc[b] + im * data.bus_shunt_suscep[b]
        if b in from_buses
            Y[b,b] += sum((y_tilde[n1,n2] + data.branch_shunt_conduc[i]/2 + im*data.branch_shunt_suscep[i]/2)/abs2(a[n1,n2]) for (i,n1,n2) in data.branches if n1 == b)
        end
        if b in to_buses
            Y[b,b] += sum(y_tilde[n1,n2] + data.branch_shunt_conduc[i]/2 + im*data.branch_shunt_suscep[i]/2 for (i,n1,n2) in data.branches if n2 == b)
        end
    end
     return Y
end

function input_matpower_data(bus_data, branch_data, gen_data, gencost_data, baseMVA)
    num_branch = length(branch_data.fbus)
    branches = Vector{Tuple{Int64,Int64,Int64}}(undef,num_branch)
    tap_ratio = Dict{Int64,Float64}()
    branch_shunt_suscep = Dict{Int64,Float64}()
    branch_shunt_conduc = Dict{Int64,Float64}()
    phase_angle = Dict{Int64,Float64}()
    branch_resist = Dict{Int64,Float64}()
    branch_react = Dict{Int64,Float64}()
    branch_current_limit = Dict{Int64,Float64}()

    for i in 1:num_branch
        branches[i] = (i,branch_data.fbus[i], branch_data.tbus[i])
        tap_ratio[i] = 1
        branch_shunt_suscep[i] = branch_data.b[i]
        branch_shunt_conduc[i] = 0.0
        phase_angle[i] = branch_data.angle[i]
        branch_resist[i] = branch_data.r[i]
        branch_react[i] = branch_data.x[i]
        branch_current_limit[i] = branch_data.rateA[i]
    end
    num_gen = length(gen_data.bus)
    gen = 1:num_gen
    network_data = power_system_data(bus_data.bus_i,bus_data.Bs,bus_data.Gs,bus_data.Pd,bus_data.Qd,bus_data.Vmax, bus_data.Vmin,
                            branches,tap_ratio,branch_shunt_suscep,branch_shunt_conduc,phase_angle,branch_resist,branch_react,branch_current_limit,
                            gen,gen_data.bus,gen_data.Pmax,gen_data.Pmin, gen_data.Qmax, gen_data.Qmin, gencost_data.c1,baseMVA)
    return network_data
end