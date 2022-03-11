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
    gen_c1::Vector{Float64}
    gen_c2::Vector{Float64}

    baseMVA::Float64
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
                            gen,gen_data.bus,gen_data.Pmax,gen_data.Pmin, gen_data.Qmax, gen_data.Qmin, gencost_data.c1, gencost_data.c2, baseMVA)
    return network_data
end