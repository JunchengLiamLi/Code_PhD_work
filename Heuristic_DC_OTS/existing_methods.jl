using PyCall, Gurobi, JuMP

nx = pyimport("networkx")
np = pyimport("numpy")

py"""
def gen_rand_span_tree(G):
    for i,j in G.edges():
        G[i][j]['weight'] = np.random.rand()
    return nx.tree.minimum_spanning_tree(G)
"""
py"""
def add_weights(T, distmx):
    for i,j in T.edges():
        T[i][j]['weight'] = distmx[i,j]
    return T
"""

function OTS(busData, branchData, generatorData, bigM, idx_unswitchable_lines, cardinality_limit, timeLimit, gap, relaxation)
    myMod = Model( optimizer_with_attributes(Gurobi.Optimizer, "Seed" => 1, "TimeLimit" => timeLimit, "MIPGap" => gap) );
    
    power_demand = Dict{Int64, Float64}()
    for i in 1:length(busData.index)
        power_demand[busData.index[i]] = busData.demand[i]
    end

    # Define variables
    #--------------------------------------------------------------------------------
    @variable(myMod, production[generatorData.index], lower_bound = 0)

    branches = Array{Tuple{Int64,Int64,Int64},1}()
    switchable_branches = Array{Tuple{Int64,Int64,Int64},1}()
    unswitchable_branches = Array{Tuple{Int64,Int64,Int64},1}()
    for i in 1:length(branchData.index)
        push!(branches, (branchData.index[i], branchData.node1[i], branchData.node2[i]) )
        if i in idx_unswitchable_lines
            push!(unswitchable_branches, (branchData.index[i], branchData.node1[i], branchData.node2[i]) )
        else
            push!(switchable_branches, (branchData.index[i], branchData.node1[i], branchData.node2[i]) )
        end
    end

    powerFlow = Dict{Tuple{Int64, Int64, Int64}, VariableRef}()
    for branch in branches
        global powerFlow[branch[1], branch[2], branch[3]] = @variable(myMod)
    end

    branchOpen = Dict{Tuple{Int64, Int64, Int64}, VariableRef}()
    for branch in switchable_branches
        if relaxation
            global branchOpen[branch[1], branch[2], branch[3]] = @variable(myMod,lower_bound = 0,upper_bound = 1)
        else
            global branchOpen[branch[1], branch[2], branch[3]] = @variable(myMod, binary = true)
        end
    end

    @variable(myMod, volAng[busData.index])
    #------------------------------------------------------------------------------------

    # Define constraints
    #-----------------------------------------------------------------------------------
    # power flow balance for generator bus
    for g in generatorData.index
        bus = generatorData.bus[g]
        @constraint(myMod,
            sum(powerFlow[index_, node1_, node2_] for (index_, node1_, node2_) in branches if node2_ == bus)
            - sum(powerFlow[index_, node1_, node2_] for (index_, node1_, node2_) in branches if node1_ == bus)
            + production[g]
            == power_demand[bus]
            )
    end

    # power flow balance for non-generator bus
    nonGenBuses = setdiff(busData.index, generatorData.bus)
    for bus in nonGenBuses
        @constraint(myMod,
            sum(powerFlow[index_, node1_, node2_] for (index_, node1_, node2_) in branches if node2_ == bus)
            - sum(powerFlow[index_, node1_, node2_] for (index_, node1_, node2_) in branches if node1_ == bus)
            == power_demand[bus]
            )
    end

    # Kirchhoff's Laws
    @constraint(myMod, [(index, node1, node2) in switchable_branches],
                -bigM[index]*(1-branchOpen[index, node1, node2])
                <= powerFlow[index, node1, node2]
                - (volAng[node1] - volAng[node2])/branchData.reactance[index]
                )

    @constraint(myMod, [(index, node1, node2) in switchable_branches],
                bigM[index]*(1-branchOpen[index, node1, node2])
                >= powerFlow[index, node1, node2]
                - (volAng[node1] - volAng[node2])/branchData.reactance[index]
                )
    @constraint(myMod, [(index, node1, node2) in unswitchable_branches],
                powerFlow[index, node1, node2]
                == (volAng[node1] - volAng[node2])/branchData.reactance[index]
                )
    # Power transmission limits
    @constraint(myMod, [(index, node1, node2) in switchable_branches],
                powerFlow[index, node1, node2]
                + branchOpen[index, node1, node2] * branchData.transLimit[index]
                >= 0
                )

    @constraint(myMod, [(index, node1, node2) in switchable_branches],
                powerFlow[index, node1, node2]
                - branchOpen[index, node1, node2] * branchData.transLimit[index]
                <= 0
                )
    @constraint(myMod, [(index, node1, node2) in unswitchable_branches],
                powerFlow[index, node1, node2]
                <=  branchData.transLimit[index]
                )
    @constraint(myMod, [(index, node1, node2) in unswitchable_branches],
                powerFlow[index, node1, node2]
                >=  -branchData.transLimit[index]
                )

    # Power production limits
    @constraint(myMod, [g in generatorData.index], production[g] <= generatorData.pmax[g])

    # Cardinality constraint
    @constraint(myMod, sum(1-branchOpen[edgeIdx, node1, node2] for (edgeIdx, node1, node2) in switchable_branches) <= cardinality_limit)
    #-----------------------------------------------------------------------------------

    @objective(myMod, Min, sum(production[g]*generatorData.cost[g] for g in generatorData.index))
    # set_silent(myMod)
    optimize!(myMod)
    if termination_status(myMod) != MOI.INFEASIBLE
        num_lines_open = length(switchable_branches) - sum(value(branchOpen[edgeIdx, node1, node2])  for (edgeIdx, node1, node2) in switchable_branches )
        return termination_status(myMod), objective_value(myMod), objective_bound(myMod), relative_gap(myMod), solve_time(myMod), num_lines_open
    else
        return termination_status(myMod), missing, missing, missing, solve_time(myMod), missing
    end
end

function OTS_with_random_fixed_ST(busData, branchData, generatorData, cardinality_limit, timeLimit, gap, relaxation)
    # Construct the graph
    graph = nx.Graph()
    distmx = Dict()
    nodes_to_idx = Dict()
    for idx in Branch_Data.index
        fbus = Branch_Data.node1[idx]
        tbus = Branch_Data.node2[idx]
        nodes_to_idx[fbus,tbus] = idx
        nodes_to_idx[tbus,fbus] = idx
        graph.add_edge(fbus, tbus, weight=Branch_Data.transLimit[idx]*Branch_Data.reactance[idx])
        distmx[fbus, tbus] = Branch_Data.transLimit[idx] * Branch_Data.reactance[idx]
        distmx[tbus, fbus] = Branch_Data.transLimit[idx] * Branch_Data.reactance[idx]
    end
    
    # generate a random spanning tree
    T = gen_rand_span_tree(graph)
    T = add_weights(T, distmx)
    
    span_tree_branch_idx = Vector{Int64}()
    for (i,j) in T.edges()
        push!(span_tree_branch_idx, nodes_to_idx[i,j])
    end
    
    sp_weight = Vector{Float64}()
    for (idx, fbus, tbus) in branches
        shortest_path = nx.dijkstra_path(T, fbus, tbus)
        sp_edges = Vector{Tuple{Int64,Int64}}()
        for i in 1:length(shortest_path)-1
            push!(sp_edges, (shortest_path[i],shortest_path[i+1]))
        end
        push!(sp_weight, sum(distmx[f,t] for (f,t) in sp_edges)) 
    end 
    bigM = sp_weight ./ branchData.reactance
    OTS(busData, branchData, generatorData, bigM, span_tree_branch_idx, cardinality_limit, timeLimit, gap, relaxation)
end




