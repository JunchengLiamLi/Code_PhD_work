using JuMP, Gurobi

function MTZ_longest_path(network_data::power_system_data, startNode, endNode, time_limit, gap)
    lwpp = Model(optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => time_limit, "MIPGap" => gap))
    nodes = network_data.buses
    edges = Vector{Tuple{Int64,Int64}}()
    for (i, n1, n2) in network_data.branches
        if !((n1, n2) in edges || (n2, n1) in edges)
            push!(edges, (n1,n2))
        end
    end
    
    # x[i,j] = 1 if the flow goes from i to j, x[i,j] = 0 otherwise
    @variable(lwpp, x[nodes, nodes], binary=true)

    # u represent the position of the node in the tour
    remaining_nodes = filter(x->x≠startNode, nodes)
    @variable(lwpp, u[remaining_nodes])

    # flow balance: what goes in must goes out
    @constraint(lwpp, [node in nodes; node != endNode && node != startNode],
        sum(x[f_node, node] for f_node in nodes) == sum(x[node, f_node] for f_node in nodes) )    
    
    # flow balance for the starting node
    @constraint(lwpp, sum(x[node,startNode] for node in nodes) == 0)
    @constraint(lwpp, sum(x[startNode,node] for node in nodes) == 1)

    # flow balance for the ending node
    @constraint(lwpp, sum(x[node,endNode] for node in nodes) == 1)
    @constraint(lwpp, sum(x[endNode,node] for node in nodes) == 0)

    # if an edge (i,j) does not exist, then x[i,j] = x[j,i] = 0
    for n1 in nodes, n2 in nodes
        if !( (n1,n2) in edges || (n2,n1) in edges )
            @constraint(lwpp, x[n1,n2] == 0)
            @constraint(lwpp, x[n2,n1] == 0)
        end
    end

    # no cancellation
    @constraint(lwpp, [(fromNode, toNode) in edges], x[fromNode, toNode] + x[toNode, fromNode] <= 1 )
    
    # at most one edge in for each node
    @constraint(lwpp, [node in nodes], sum(x[f_node,node] for f_node in nodes) <= 1)

    # at most one edge out for each node
    @constraint(lwpp, [node in nodes], sum(x[node,t_node] for t_node in nodes) <= 1)

    # the edge (startNode, endNode) is disconnected
    # @constraint(lwpp, x[startNode,endNode] == 0)

    # MTZ constraints to eliminate subtours
    num_nodes = length(nodes)
    @constraint(lwpp, [n1 in remaining_nodes, n2 in remaining_nodes; n1 ≠ n2], (num_nodes-1)*x[n1,n2] + u[n1] - u[n2] <= num_nodes-2 )
    
    @objective(lwpp, Max, sum(x[fromNode,toNode] + x[toNode,fromNode] for (fromNode, toNode) in edges) )
    set_silent(lwpp)
    optimize!(lwpp)

    return objective_value(lwpp)
end