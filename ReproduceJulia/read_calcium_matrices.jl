using Plots
using MAT
using StatsBase
"""
Of course the absolute paths below will need to be wrangled to match your directory tree.
"""
read_spike_dense = matread("../datasets/M1_d1A_S.mat")["GC06_M1963_20191204_S1"]["Transients"]["Raster"]
FPS = matread("../datasets/M1_d1A_S.mat")["GC06_M1963_20191204_S1"]["Movie"]["FPS"]
frame_width = 1.0/FPS #0.08099986230023408 #second, sample_rate =  12.3457#Hz

"""
A method to re-represent dense boolean vectors as a two dense vectors of spikes, and times.
spikes is a matrix with regularly sampled windows, populated by spikes, with calcium spikes.
"""
function convert_bool_matrice_to_ts(read_spike_dense::Matrix{Bool}, frame_width::Real)
    nodes = UInt32[]
    times = Float32[]
    for (indy,row) in enumerate(eachrow(read_spike_dense))
        for (indx,x) in enumerate(row)
            if x
                push!(nodes,indy)
                push!(times,indx*frame_width)                
            end
        end
    end
    whole_duration = length(read_spike_dense[1,:])*frame_width
    (nodes::Vector{UInt32},times::Vector{Float32},whole_duration::Real)
end

"""
A method to get collect the Inter Spike Intervals (ISIs) per neuron, and then to collect them together to get the ISI distribution for the whole cell population
Also output a ragged array (Array of unequal length array) of spike trains. 
"""
function create_ISI_histogram(nodes::Vector{UInt32},times::Vector{Float32})
    spikes_ragged = []
    global_isis =Float32[] # the total lumped population ISI distribution.
    isi_s = []
    numb_neurons=Int(maximum(nodes))+1 # Julia doesn't index at 0.
    @inbounds for n in 1:numb_neurons
        push!(spikes_ragged,[])
    end
    @inbounds for i in 1:numb_neurons
        for (n,t) in zip(nodes,times)
            if i==n
                push!(spikes_ragged[i],t)
            end
        end
    end
    @inbounds for (i, times) in enumerate(spikes_ragged)
        push!(isi_s,[])
        for (ind,x) in enumerate(times)
            if ind>1
                isi_current = x-times[ind-1]
                push!(isi_s[i],isi_current)
            end
        end
        append!(global_isis,isi_s[i])
    end
    (global_isis,spikes_ragged)
end
(nodes,times,whole_duration) = convert_bool_matrice_to_ts(read_spike_dense,frame_width)
global_isis,spikes_ragged = create_ISI_histogram(nodes,times)

"""
Visualize one epoch, as a spike train raster and then an ISI histogram.
"""
p1 = Plots.plot()
Plots.scatter!(p1,times,nodes,legend = false,markersize = 0.8,markerstrokewidth=0,alpha=0.8, bgcolor=:snow2, fontcolor=:blue,xlabel="time (Seconds)",ylabel="Cell Id")
#savefig("scatter_plot.png")
b_range = range(minimum(global_isis), mean(global_isis)+std(global_isis), length=21)
p2 = Plots.plot()
Plots.histogram!(p2,global_isis, bins=b_range, normalize=:pdf, color=:gray,xlim=[0.0,mean(global_isis)+std(global_isis)])
Plots.plot(p1,p2)
savefig("Spike_raster_and_ISI_bar_plot.png")
