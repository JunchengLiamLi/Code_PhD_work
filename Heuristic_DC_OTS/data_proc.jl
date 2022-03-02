
using CSV, DataFrames, XLSX
using Statistics
using Plots

function scat_plot_kSP(network, load_factor, input_data, conv_levelS, best_para, shapeS, labelS, lg, colorS, xscaleS, output_file)
    # input the data
    #--------------------------------------------------------------
    if network == 2 # analyze 300-bus data
        load_data = input_data[input_data.load_factor .== load_factor,:]
        p = plot()
        for i in 1:length(conv_levelS)
            conv_level = conv_levelS[i]
            plot_data = load_data[load_data.conv_level .== conv_level,:]
            plot!(p, plot_data.avg_time, plot_data.avg_rel_gap .* 100, seriestype = :scatter, xlims = xscale, color = :blue, legend = lg, markershape = shapeS[i], label = labelS[i])
        end
        for i in 1:length(best_para)
            (p1,p2,p3) = best_para[i]
            best_para_data = load_data[load_data.max_path .== p1, :]
            best_para_data = best_para_data[best_para_data.max_edge .== p2, :]
            best_para_data = best_para_data[best_para_data.conv_level .== p3, :]
            plot!(p, best_para_data.avg_time, best_para_data.avg_rel_gap .* 100, seriestype = :scatter, xlims = xscale, color = colorS[i], legend = lg, markershape = :star5, markersize = 8, label = "($p1,$p2,$p3)")
        end
        xlabel!(p, "average computation time (seconds)")
        ylabel!(p, "average relative gap(%)")
        if network == 1
            output_file_name = output_file*"_118bus_$(load_factor).png"
        elseif network == 2
            output_file_name = output_file*"_300bus_$(load_factor).png"
        end
        savefig(p, output_file_name)
    end
    #--------------------------------------------------------------
end

function scat_plot_kNN(network, load_factor, input_data, best_para, lg, colorS, xscale, output_file)
    # input the data
    #--------------------------------------------------------------
    load_data = input_data[input_data.load_factor .== load_factor,:]
    p = plot(load_data.avg_time, load_data.avg_rel_gap .* 100, seriestype = :scatter, xlims = xscale, color = :blue, markershape = :circle, label = "")
    for i in 1:length(best_para)
        (p1,p2,p3) = best_para[i]
        best_para_data = load_data[load_data.neigh .== p1, :]        
        best_para_data = best_para_data[best_para_data.per_lines .== p2, :]
        best_para_data = best_para_data[best_para_data.scal .== p3, :]
        plot!(p, best_para_data.avg_time, best_para_data.avg_rel_gap .* 100, seriestype = :scatter, xlims = xscale, color = colorS[i], legend = lg, markershape = :star5, markersize = 8, label = "($p1,$(Int32(p2*100)),$p3)")
    end
    xlabel!(p, "average computation time (seconds)")
    ylabel!(p, "average relative gap(%)")
    if network == 1
        output_file_name = output_file*"_118bus_$(load_factor).png"
    elseif network == 2
        output_file_name = output_file*"_300bus_$(load_factor).png"
    end
    savefig(p, output_file_name)
    #--------------------------------------------------------------
end

# analyze OTS results with the same parameters for each class of instances
# if method == 1, analyze k-shortest-path method
# if method == 2, analyze kNN simulation method
# if network == 1, use 118 bus data
# if network == 2, use 300 bus data
function analyze_OTS_results(network, method, load_scale_factors, bm_load_scale_factors, instances, parameters, input_df, bm_df, output_file)
    if method == 1 # analyze k-shortest-path method
        output_df = DataFrame(load_factor=[], bigM=[], max_path=[], max_edge=[], conv_level=[], avg_rel_gap=[], avg_time=[], num_not_opt=[], not_opt_gap=[])
    elseif  method == 2 # analyze kNN simulation method
        output_df = DataFrame(load_factor=[], bigM=[], neigh=[], per_lines=[], scal=[], avg_rel_gap=[], avg_time=[], num_not_opt=[], not_opt_gap=[])
    end
    for l in 1:length(load_scale_factors)
        if method == 1 # analyze k-shortest-path method
            for (max_path,max_edge,conv) in parameters
                sum_rel_gap = 0
                num_not_opt = 0
                not_opt_gap = 0
                sum_time = 0
                for i in instances
                    # build the dataframe with the information of relative gaps and computational time
                    if network == 1
                        instance = "IEEE118_P_$(load_scale_factors[l])_$i"
                        bm_instance = "IEEE118_P_$(bm_load_scale_factors[l])_$i"
                    elseif network == 2
                        instance = "IEEE300_P_$(load_scale_factors[l])_$i"
                        bm_instance = "IEEE300_P_$(bm_load_scale_factors[l])_$i"
                    end
                    inst_df = input_df[input_df.data_instance .== instance, :]      
                    para_df1 = inst_df[inst_df.max_path .== max_path,:]
                    para_df2 = para_df1[para_df1.max_edge .== max_edge,:]
                    para_df = para_df2[para_df2.conv_level .== conv, :]
                    inst_bm_df = bm_df[bm_df.data_instance .== bm_instance, :] # the benchmark
                    # summarize the information
                    sum_time += para_df.time[1]
                    sum_rel_gap += ((para_df.obj_val - inst_bm_df.obj_val) / inst_bm_df.obj_val)[1]
                    if para_df.sol_status[1] == "TIME_LIMIT"
                        num_not_opt += 1
                        not_opt_gap += para_df.gap[1]
                    end
                end
                
                output_row = []
                push!(output_row, load_scale_factors[l])
                push!(output_row, "ksp")
                push!(output_row,max_path)
                push!(output_row,max_edge)
                push!(output_row,conv)
                avg_rel_gap = sum_rel_gap / length(instances)
                push!(output_row, avg_rel_gap)
                avg_time = sum_time / length(instances)
                push!(output_row, avg_time)
                push!(output_row, num_not_opt)
                avg_not_opt_gap = not_opt_gap/num_not_opt
                push!(output_row, avg_not_opt_gap)
                push!(output_df, output_row)
            end
        elseif method == 2 # analyze kNN method
            for (neigh, per, scal) in parameters
                sum_rel_gap = 0
                num_not_opt = 0
                not_opt_gap = 0
                sum_time = 0
                for i in instances
                    # build the dataframe with the information of relative gaps and computational time
                    if network == 1
                        instance = "IEEE118_P_$(load_scale_factors[l])_$i"
                        bm_instance = "IEEE118_P_$(bm_load_scale_factors[l])_$i"
                    elseif network == 2
                        instance = "IEEE300_P_$(load_scale_factors[l])_$i"
                        bm_instance = "IEEE300_P_$(bm_load_scale_factors[l])_$i"
                    end
                    inst_df = input_df[input_df.data_instance .== instance, :]
                    para_df1 = inst_df[inst_df.size_neigh .== neigh,:]
                    para_df2 = para_df1[para_df1.per_lines_rm .== per,:]
                    para_df = para_df2[para_df2.scal_factor .== scal, :]
                    inst_bm_df = bm_df[bm_df.data_instance .== bm_instance, :] # the benchmark
                    # summarize the information
                    sum_time += para_df.time[1]
                    sum_rel_gap += ((para_df.obj_val - inst_bm_df.obj_val) / inst_bm_df.obj_val)[1,1]
                    if para_df.sol_status[1] == "TIME_LIMIT"
                        num_not_opt += 1
                        not_opt_gap += para_df.gap[1]
                    end
                end
                
                output_row = []
                push!(output_row, load_scale_factors[l])
                push!(output_row, "kNN")
                push!(output_row, neigh)
                push!(output_row, per)
                push!(output_row, scal)
                avg_rel_gap = sum_rel_gap / length(instances)
                push!(output_row, avg_rel_gap)
                avg_time = sum_time / length(instances)
                push!(output_row, avg_time)
                push!(output_row, num_not_opt)
                avg_not_opt_gap = not_opt_gap/num_not_opt
                push!(output_row, avg_not_opt_gap)
                push!(output_df, output_row)
            end
        end
    end
    CSV.write(output_file, output_df)
end

function analyze_opt_gap(cardinality_limits, num_instances, input_df, output_file)
    output_mx = zeros(Float64, length(cardinality_limits), num_instances)
    for (idx, card) in enumerate(cardinality_limits)
        input_card = input_df[input_df.card .== card, :]
        for i in 1:nrow(input_card)
            if input_card.sol_status[i] == "TIME_LIMIT"
                output_mx[idx, i] = input_card.gap[i]
            end
        end
    end
    output_df  = DataFrame(output_mx,:auto)
    CSV.write(output_file, output_df)
end

function analyze_fuller(network, fuller_input, lwp_input, load_factors, lwp_load_factors, num_instances)
    output_df = DataFrame(instance = String[], gap = Float64[])
    for l in 1:length(load_factors), i in 1:num_instances
        if network == 1
            instance = "IEEE118_P_$(load_factors[l])_$i"
            lwp_instance = "IEEE118_P_$(lwp_load_factors[l])_$i"
        elseif network == 2
            instance = "IEEE300_P_$(load_factors[l])_$i"
            lwp_instance = "IEEE300_P_$(lwp_load_factors[l])_$i"
        end
        fuller_df = fuller_input[fuller_input.data_instance .== instance,:]
        if fuller_df.obj_val[1] != "none"
            instance_df = lwp_input[lwp_input.data_instance .== lwp_instance,:]
            rel_gap = (parse(Float64,fuller_df.obj_val[1]) - instance_df.obj_val[1]) / instance_df.obj_val[1]
            output_row = [instance, rel_gap]
            push!(output_df, output_row)
        end
    end
    return output_df
end

# kSP scatter plots
#=
# input data
#--------------------------------------------------------------------------------------
input_data = DataFrame(CSV.File("results_300/analyze_ksp_para_tuning_300.csv"))
network = 2
load_factor = 1.05
#---------------------------------------------------------------------------------------

#conservative levels
#----------------------------------------------------------------------------------
conv_levelS = [2,3,4,5]
shapeS = [:diamond, :ltriangle, :circle, :xcross]
labelS = ["l=2","l=3","l=4","l=5"]
#----------------------------------------------------------------------------------

# best parameters
#-------------------------------------------------------------------------------
best_para = [(15,7,5), (15,5,4)]
colorS = [:pink, :yellow]
#-------------------------------------------------------------------------------

xscale = (30,80)
lg = :left
output_file = "plots/scat_OTS_kSP_para_tune_n"

scat_plot_kSP(network, load_factor, input_data, conv_levelS, best_para, shapeS, labelS, lg, colorS, xscale, output_file)
=#

# kNN scatter plots
#= input data
#--------------------------------------------------------------------------------------
input_data = DataFrame!(CSV.File("merged_branch_118/new_analyze_kNN_para_tuning_118.csv"))
network = 1
load_factor = 105
#---------------------------------------------------------------------------------------


# best parameters
#-------------------------------------------------------------------------------
best_para = [(5,0.1,15),(4,0.1,15),(4,0.05,10)]
colorS = [:red, :pink, :yellow]
#-------------------------------------------------------------------------------
=#
#=
xscale = (0,200)
lg = :topleft
output_file = "codes_and_test_instances/plots/scat_OTS_kNN_para_tune"

scat_plot_kNN(network, load_factor, input_data, best_para, lg, colorS, xscale, output_file)
=#

# analyze k values
#=
#-------------------------------------------------------------------------------------------------------------------
e_max = [3,5,7]
k_max = [11,14,17,20]
num_edges, num_instances = 180, 20
for idx_e_max in 1:length(e_max), idx_k_max in 1:length(k_max)
    worksheet_mx = Matrix{Int64}(undef, num_edges, num_instances)
    for instance in 1:num_instances
        k_val_df = DataFrame(CSV.File("results_118/new_ksp_k_values/load_100_percent_instance_$(instance).csv"))
        k_val_col = k_val_df[:,idx_e_max+1]
        for i in 1:length(k_val_col)
            k_val_col[i] = min(k_val_col[i], k_max[idx_k_max])
        end
        worksheet_mx[:,instance] = k_val_col
    end
    worksheet_df = DataFrame(worksheet_mx, :auto)
    XLSX.openxlsx("results_118/new_118N_k_val.xlsx",mode="rw") do xf
        XLSX.addsheet!(xf,"k_max_$(k_max[idx_k_max])_e_max_$(e_max[idx_e_max])")
        sheet = xf["k_max_$(k_max[idx_k_max])_e_max_$(e_max[idx_e_max])"]
        for r in 1:size(worksheet_df,1), c in 1:size(worksheet_df,2)
             sheet[XLSX.CellRef(r , c )] = worksheet_df[r,c]
        end
   end
end
#------------------------------------------------------------------------------------------------------------------------
=#

# analyze OTS results
#=
#----------------------------------------------------------------------------------------------------------------------------------
network = 2
method = 1
load_factorS = [1.0]
bm_load_factorS = [100]
instances = 1:20
parameter_list = []
for p1 in [15,20,25,30], p2 in [5,7,9], p3 in [5]
    push!(parameter_list, (p1,p2,p3))
end
input_results = DataFrame(CSV.File("results_300/OTS_ksp_para_tuning_new.csv"))
bm_df = DataFrame(CSV.File("results_300/OTS_lwp_para_tuning.csv"))
output_file = "results_300/analyze_OTS_para_tune_new.csv"
analyze_OTS_results(network, method, load_factorS, bm_load_factorS, instances, parameter_list, input_results, bm_df, output_file)
#----------------------------------------------------------------------------------------------------------------------------------
=#
