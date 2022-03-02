using CSV, DataFrames, Random

# input data
bus_data = DataFrame(CSV.File("Data/300bus_data/IEEE300_bus.csv"))

function generate_new_bus_data(bus_data, i, l, num_buses)
    # create new data
    new_bus_data = copy(bus_data)
    random_scales = rand(num_buses)./10 .- 0.05 .+ 1
    new_bus_data.demand = new_bus_data.demand .* random_scales .* (l/100)
    
    # save new data
    CSV.write("Data/300bus_data/para_tune_loads/load_factor_$(l)_percent/Bus_Data_$i.csv", new_bus_data)
end

new_data_indices = 21:100
loads = [95,100,105]
num_buses = 300
for i in new_data_indices, l in loads
    generate_new_bus_data(bus_data, i, l, num_buses)
end
