### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 2437222a-26f1-11ee-30b4-a96159720d68
using Revise
using StatsPlots, UnROOT, StatsBase, Polynomials, LinearAlgebra, Printf
using FHist, MPThemes, DataFramesMeta, Distributions, LaTeXStrings 

# ╔═╡ 243c0132-26f1-11ee-2333-7ba24e313a18
ENV["COLUMNS"] = 2000
ENV["LINES"] = 20

# ╔═╡ 243ce74e-26f1-11ee-3178-1bf57ca82af4
md"Custom functions are placed in the MiscFuncs.jl file."

# ╔═╡ 243ce824-26f1-11ee-10d6-d1044295b35e
include("MiscFuncs.jl")
using .MiscFuncs

# ╔═╡ 243ce85e-26f1-11ee-3a2c-6f32f20a1cbe
gr()
default(fmt = :jpg)
theme(
    :dao;
    size           = (1200, 800),
    legend         = :topleft,
    guidefontsize  = 20,
    tickfontsize   = 16,
    titlefontsize  = 20,
    legendfontsize = 20,
    left_margin    = 8Plots.mm,
    right_margin   = 8Plots.mm,
    top_margin     = 8Plots.mm,
    bottom_margin  = 6Plots.mm,
    dpi            = 200,
    :colorbar_titlefontsize => 16,
    widen = :false
);

# ╔═╡ 243ce872-26f1-11ee-37fe-d1e78aba7c33
baseDir = "/home/shoram/Work/PhD_Thesis/Job15/AngularCorrelations/"

# ╔═╡ 243ce884-26f1-11ee-32c6-532ac1a54626
figDir = joinpath("/media/shoram/Extra SSD/CernBox/Work/Presentations/20221006_SuperNEMO_Analysis_Meeting/Figs")

# ╔═╡ 243ce89a-26f1-11ee-108d-4b137230b0a0
f = ROOTFile(
    baseDir*"AngularCorrelationAllEnergies96MilEvents.root",
);
tree = DataFrame(LazyTree(f, "tree", keys(f["tree"])));
f= nothing; # free file from memory with gc


# ╔═╡ 243ce8c2-26f1-11ee-236d-e9bd0a43478c
md"### ``@transform`` adds a column ``:ESum`` to the ``tree`` which contains the sum of the electron energies"

# ╔═╡ 243ce8d6-26f1-11ee-26e1-b5225e40b167
@transform! tree :ESum = :reconstructedEnergy2 + :reconstructedEnergy1;

# ╔═╡ 243ce8fe-26f1-11ee-24f4-bfcdc4137274
md"### Initializing constants."

# ╔═╡ 243ce930-26f1-11ee-31cd-f17ad3702d2a

dEmitted = 1 # dθdif in degrees
nBins    = Int(180 / dEmitted)
minAngle = 0
maxAngle = 180
binWidth = maxAngle / nBins

minEnergy = 500
maxEnergy = 3500
dEnergy   = 500

xPts = minAngle:dEmitted:maxAngle-dEmitted

dϕ = dEmitted               # step in ϕ, if same as bin width
sign = "p"                  # sign in get_cut_edges function
maxSteps = Int(180 / dϕ)    # max number of steps (slices)

# ╔═╡ 243ce93a-26f1-11ee-3957-976266e37bae
colors = [palette(:seaborn_bright)[i] for i in 1:length(palette(:seaborn_bright))];

# ╔═╡ 243ce95a-26f1-11ee-0468-c92a72f82437
md"The 2d Histogram of $\phi$ vs $\theta$ is defined to be $f(\theta,\phi)$. For each combination of $\phi$ and $\theta$, the bin number is obtained as corresponding value of $f(\theta,\phi)$."

# ╔═╡ 243ce980-26f1-11ee-2d0c-851fdaae8542
rho(_cosdTheta) = 0.5 - 0.5 * _cosdTheta

# @transform! tree :weights = Weights(rho.(cosd.(:thetaEmitted))) # weights by pdf
@transform! tree :weights = 1                                     # unweighted



# ╔═╡ 243ce9a8-26f1-11ee-2b76-6bba0d7bf2ab
md"""
Quantitative analysis - *the g(k) method*
===


Following the qualitative analysis presented above, a more quantitative approach will now be introduced. The goal of the so-called *g(k) method* is to quantitatively describe correlation between $\theta$ and $\phi$. This method is used for evaluation of proposed data-cuts. 

To derive the method, two figures are presented:
1. $f(\theta, \phi)$ - a 2D histogram with $\theta$ distribution on the x-axis and $\phi$ distribution on the y-axis
2. $g(k)$ - a 1D histogram with *so-called* $k$-line on the x-axis, where bin heights correspond to the integral over individual $k_i$ lines. $k$-lines are defined as diagonal lines in $f(\theta, \phi)$, fulfilling condition $k_i=\phi - \theta; (\phi - \theta) \in b_i$; binned as necessary, in general binning of $\Delta\phi = \Delta\theta = 1^{\circ}$ is used. 
"""

# ╔═╡ 243ce9c6-26f1-11ee-35e9-852f3fdeab95
md"## $f(\theta, \phi)$"

# ╔═╡ 243ce9e4-26f1-11ee-0f12-111fbfb0e83c
nBins = 180

h2d1 = histogram2d(
    tree.thetaEmitted,
    tree.thetaEscaped;
    nbins          = (nBins, nBins),
    xlabel         = "ϑ",
    ylabel         = "φ",
    legend         = :topright,
    title          = string("\nϑ vs φ distribution, bin width Δ = 1°"),
    lims           = (0, 180),
#     aspect_ratio   = 1,
    c              = :coolwarm,
    colorbar_title = "\ncounts [#/°]",
    thickness_scaling = 1.3,
    size= (1200,800)
)

plot!(y, c = :black, lw= 3, ls =:solid, label = "ϑ = φ", legend = :topleft)

# ╔═╡ 243ce9f8-26f1-11ee-30ef-75c9c32de503
savefig("Figs/2DThetaVsPhiAll.pdf")

# ╔═╡ 243cea0c-26f1-11ee-14b4-b3ead850f465
md"The correlation of $\theta$ and $\phi$ in $f(\theta, \phi)$ shows an S-shaped structure. There are two *hotspots* visible in the 2D histogram. First, a *smaller* hotspot is visible in the lower left quadrant of $f(\theta, \phi)$, here in-general $\phi > \theta$ - escape angles are overestimated over decay angles. The second *larger* hotspot is visible in the upper right quadrant of $f(\theta, \phi)$ where $\phi < \theta$ - escape angles are underestimed over decay angles. "

# ╔═╡ 243cea1e-26f1-11ee-1a53-4fddd3aeafaf
md"We can look more closely at the individual horizontal slices to view the $\theta$ distribution in each $\Delta\phi$. A few sample figures are shown. One can see that the $\theta$ distribution is very wide."

# ╔═╡ 243cea7a-26f1-11ee-1a5d-d7476036e22c
ts = 1.6
h1 = stephist(
    tree[ tree.thetaEscaped .< 5.0, :thetaEmitted ], 
    weights=tree[tree.thetaEscaped .< 5.0, :weights], label ="", c =:orange,
    ylabel = "counts [#/°]",
    legend = :topright, xlims = (0,180), ylims = (0, 200), 
    xlabel ="θ", lw = 4, nbins = 180, dpi = 200,
    thickness_scaling = ts
)
vspan!([0,5], alpha = 0.2, c =:orange, label = "φ ∈ (0,5)°" )

h2 = stephist( 
    tree[ 30 .< tree.thetaEscaped .< 35, :thetaEmitted ],
    ylabel = "counts [#/°]",
    weights=tree[30 .< tree.thetaEscaped .< 35, :weights],label ="", c =:red,
    legend = :topright, xlims = (0,180), ylims = (0, 6000), 
    xlabel ="θ", lw = 4, nbins = 180, dpi = 200,
    thickness_scaling = ts
    
)
vspan!([30,35], alpha = 0.2, c =:red, label = "φ ∈ (30,35)°" )

h3 = stephist( 
    tree[ 85 .< tree.thetaEscaped .< 90, :thetaEmitted ],
    weights=tree[85 .< tree.thetaEscaped .< 90, :weights], label ="", c =:blue,
    ylabel = "counts [#/°]",
    thickness_scaling = ts,
    legend = :topright, xlims = (0,180), ylims = (0, 6000), 
    xlabel ="θ", lw = 4, nbins = 180, dpi = 200 
)
vspan!([85,90], alpha = 0.2, c =:blue, label = "φ ∈ (85,90)°" )

h4 = stephist( 
    tree[ 120 .< tree.thetaEscaped .< 125, :thetaEmitted ],
    weights=tree[120 .< tree.thetaEscaped .< 125, :weights],label ="", c =:green, legend = :topleft,
    ylabel = "counts [#/°]",
    thickness_scaling = ts,
    xlims = (0,180), ylims = (0, 9000), 
    xlabel ="θ", lw = 4, nbins = 180, dpi = 200 
)
vspan!([120,125], alpha = 0.2, c =:green, label = "φ ∈ (120,125)°" )

h5 = stephist( 
    tree[ 150 .< tree.thetaEscaped .< 155, :thetaEmitted ],
    weights=tree[150 .< tree.thetaEscaped .< 155, :weights],label ="", c =5, legend = :topleft,
    ylabel = "counts [#/°]",
    thickness_scaling = ts,
    xlims = (0,180), 
    xlabel ="θ", lw = 4, nbins = 180, dpi = 200 
)
vspan!([150,155], alpha = 0.2, c =5, label = "φ ∈ (150,155)°" )


h6 = stephist( 
    tree[ 175 .< tree.thetaEscaped .< 180, :thetaEmitted ],
    weights=tree[175 .< tree.thetaEscaped .< 180, :weights],label ="", c =6, legend = :topleft,
    ylabel = "counts [#/°]",
    thickness_scaling = ts,
    xlims = (0,180), ylims = (0, 1500), 
    xlabel ="θ", lw = 4, nbins = 180, dpi = 200 
)
vspan!([175,180], alpha = 0.2, c =6, label = "φ ∈ (175,180)°" )

plot(h1, h2, h3, h4, h5, h6, layout = @layout [a b ; c d ; e f])



# ╔═╡ 243cea90-26f1-11ee-3bc4-c51afe379f0a
let hs = [h1, h2, h3, h4, h5, h6]
    for (i,h) in enumerate(hs)
        savefig(h, "Figs/Slice$i.pdf")
    end
    
end

# ╔═╡ 243ceaac-26f1-11ee-0a5a-b36800325949
md"Furthermore we can look even closer at each slice and calculate the statistical estimators. These will be used later in the analysis."

# ╔═╡ 243ceac2-26f1-11ee-2725-9128a5488720
using Measures

# ╔═╡ 243ceafc-26f1-11ee-07d8-cffedeccbc13
function shaded_interval(thetas, lowerbound, upperbound, Δ = 5, binrange =(0:Δ:180); kwargs...)
    h = Hist1D(thetas, binrange)
    lowerQuantile = quantile(h, lowerbound)    # find the corresponding quantiles
    upperQuantile = quantile(h, upperbound)  

    minBinEdge = floor(lowerQuantile/Δ)*Δ+Δ  # find the upper edge of the bin with lower quantile
    maxBinEdge = floor(upperQuantile/Δ)*Δ  # find the lower edge of the bin with upper quantile
    binEdges = h.hist.edges[1][Int(binrange[1] + minBinEdge/Δ+1):Int(binrange[1] + maxBinEdge/Δ+1)] # get the binedges in between

    xvals = vcat(lowerQuantile, binEdges, upperQuantile)
    yvals = lookup.( h,  xvals) # find the binheights for the corresponding interval

    return plot( xvals, yvals, st =:step; kwargs... )
end

function shaded_interval!(p, thetas, lowerbound, upperbound, Δ = 5, binrange = (0:Δ:180); kwargs...)
    h = Hist1D(thetas, binrange)
    lowerQuantile = quantile(h, lowerbound)    # find the corresponding quantiles
    upperQuantile = quantile(h, upperbound)  

    minBinEdge = floor(lowerQuantile/Δ)*Δ+Δ  # find the upper edge of the bin with lower quantile
    maxBinEdge = floor(upperQuantile/Δ)*Δ  # find the lower edge of the bin with upper quantile
    binEdges = h.hist.edges[1][Int(binrange[1] + minBinEdge/Δ+1):Int(binrange[1] + maxBinEdge/Δ+1)] # get the binedges in between

    xvals = vcat(lowerQuantile, binEdges, upperQuantile)
    yvals = lookup.( h,  xvals) # find the binheights for the corresponding interval
    plot!(p, xvals, yvals, st =:step; kwargs... )
    
    return p
end

# ╔═╡ 243ceb38-26f1-11ee-1973-052f01ea98b7
hists = []
for (i, phiMin, phiMax) in zip(1:6, [10., 25., 45., 90., 135., 175.], [15., 30., 50., 95., 140., 180.])
    thetas = tree[phiMin .< tree.thetaEscaped .< phiMax, :thetaEmitted]
    stats =  get_slice_stats(phiMin, phiMax,  0.0, 3500., thetas, float(5))
        
    lb, ub = (1-0.68)/2,(1+0.68)/2 #lower and upper 1σ about median 

    h = Hist1D(thetas, (0:5:180))
#     mediansFHist = StatsBase.median(h)

    h1 = stephist( 
        thetas;
        nbins = Int(180/5), 
        xlims = (0,180), 
        ylims = (0, 35000),
        yticks = (1e4:1e4:3e4),
        label ="", 
        ylabel = ifelse(i == 1 || i ==3 || i == 5, "counts [#/°]", ""),
        c =i, 
        lw = 4,
        xlabel = ifelse(i<5, "", "ϑ [°]"),
        thickness_scaling = 1.2,
        legend = ifelse( i<4, :topright, :topleft ),
        bottom_margin = ifelse(i<5, -7mm, 2mm),
        left_margin = ifelse(i%2==0, 14mm, 2mm),
        legendfontsize = 18
        )
    
    @show phiMin, phiMax
    @show "median, mode, mead"
    @show median(h), (stats[6] - 2.5), mean(h)
    
    shaded_interval!(h1, thetas, lb, ub; fillrange = 0, fa = 0.4, fillcolor= :grey, label = "", lw=0, fillstyle =:x)
    vspan!([phiMin,phiMax], alpha = 0.8, c =i, label = "φ ∈ ($(Int(phiMin)),$(Int(phiMax)))°" )
    vline!([median(h)],c=3,lw =2, label ="")
    vline!([stats[6]-2.5],c=4,lw =2, label ="")
    vline!([mean(h)],lw =2,c=5, label ="")

    
    push!(hists, h1)
end

# ╔═╡ 243ceb6a-26f1-11ee-2883-69c9f62c651d
legendPlot = plot(
    [1 2 3], xlim = (4,5), legend = :top, framestyle = :none, label = ["median" "mode" "mean"], c = [3 4 5]);
plot!([1], [1],fillrange = 0, c= :grey,fillstyle =:x, lw = 0, label = "central 68%", legend_column =-1, size =(400,400), legendfontsize = 18, bottom_margin = -3mm)

# ╔═╡ 243ceb7e-26f1-11ee-0789-9530a207fa5f
layout = @layout [a{0.05h};
                  b c;
                  e f ;
                  g h]

# ╔═╡ 243ceba6-26f1-11ee-1260-65be19b68358
plot( legendPlot, hists[1], hists[2], hists[3], hists[4],  hists[5], hists[6], layout = layout, size =(1200, 1400),
        thickness_scaling = 1.2, guidefont=12, xtickfont=12, ytickfont=12, legendfont=14, gridalpha = 0.1, plot_title = "ϑ distribution for given partitions", top_margin = 2mm)

# ╔═╡ 243cebba-26f1-11ee-26e7-3b4cb7fdc77c
savefig("Figs/Slices_stats.pdf")

# ╔═╡ 243cebec-26f1-11ee-345d-950cb82da77d
let dϕ = 1., len = Int(180/dϕ), means= fill(-1., len), 
    medians = fill(-1., len), modes = fill(-1., len), 
    lowers = fill(-1., len),uppers = fill(-1.,len)
    
    xrange = 1:dϕ:180
    lb, ub = (1 - 0.68)/2, (1 + 0.68)/2
    for (i, ϕ) in enumerate(1:dϕ:180)
        cutEdges1 = get_cut_edges(ϕ - 1, 1, dϕ, "p")                    # provides the lower and upper cut 
        sdf       = @chain tree begin                                   # filter out the dataframe
            @select(:thetaEscaped, :thetaEmitted, :weights)                       # keeps only the two angles columns
            @subset((cutEdges1[1] .<= :thetaEscaped .<= cutEdges1[2]))  # keeps only rows where ϕ is within the cut edges
        end

        stats = get_slice_stats(
                                cutEdges1[1],  # ϕmin
                                cutEdges1[2],  # ϕmax
                                0,             # Emin
                                3500,          # Emax
                                sdf.thetaEmitted, 
                                dϕ,
                                Weights(sdf.weights)
                                )

        means[i] = stats[5] 
        modes[i] = stats[6] - dϕ/2
        medians[i] =  quantile(sdf.thetaEmitted, 0.5) 
        lowers[i] = quantile(sdf.thetaEmitted, lb)
        uppers[i] = quantile(sdf.thetaEmitted, ub)
    end
    plot(xrange, uppers, fillrange=lowers, label = "", c = :grey, fa = 0.3, lw = 0)
    plot!(xrange, xrange, c =:black, lw =3, label = "")
    scatter!(xrange, means, label = "mean", c= 1, ms =5, lims = (0,180), ticks = 0:30:180, xlabel = L"\mathrm{partition}~\Delta\varphi_i", ylabel = "estimator")
    scatter!(xrange, modes, label = "mode", c=4, ms =5)
    scatter!(xrange, medians, label = "median", c=3, ms =5, title ="Statistical estimators for each partition" )
    plot!([1], [1],fillrange = 0, c=:grey, lw = 0, fa=0.3, label = "central 68%", legend_column =2, legend =:bottomright) 

end

# ╔═╡ 243cec00-26f1-11ee-24e7-bf8fe02ade17
savefig("Figs/Slices_stats_all.pdf")

# ╔═╡ 243cec14-26f1-11ee-0d28-033b7e625dd2
md"## k-lines"

# ╔═╡ 243cec2a-26f1-11ee-16e8-25f63d7cb2dc
md"In the figure above of $f(\theta, \phi)$, a reference k = 0 line is shown by black dashed line. This line represents the bins where $\phi = \theta$, perfect correlation. In the figure below, two more k-lines are depicted. $k = -20$ and $k = +20$ lines are show in red and blue, respectively. These lines in turn represent $\phi - \theta = -20$ and $\phi - \theta = 20$. "

# ╔═╡ 243cec50-26f1-11ee-2271-0596a4eec19c
pts = 0:0.1:180

# plot!(h2d1, y, c = :black, lw= 3, ls =:dash, label = "k = 0", legend = :topleft)
plot!(h2d1, pts, pts .+ 20,  lw= 4, c  =:red , label = "k = +20")
plot!(h2d1, pts, pts .- 20,  lw= 4, c  =:blue, label = "k = -20")

# ╔═╡ 243cec64-26f1-11ee-28e4-798cb7b12f20
savefig("Figs/Theta_vs_Phi_with_k_lines.pdf")

# ╔═╡ 243cec78-26f1-11ee-196b-cbd9a31cdf4c
md"## $g(k)$"

# ╔═╡ 243cec8a-26f1-11ee-387e-6f6064f3b56f
md"""
As stated earlier, the $g(k)$ is a 1D histgoram of the integrals over the k-lines. The total number of k-lines to integrate $f(\theta, \phi)$ over is equal to $(180/\Delta\phi) * 2 + 1$, dependent on the binwidth $\Delta\phi$. To avoid double binning, $g(k)$ is not calculated from $f(\theta, \phi)$ itself, but rather from the definition of $k$-lines. Thus for each event, $\phi - \theta$ is calculated and **then** binned in the 1D histogram $g(k)$. The figure below shows $g(k)$ of the original $f(\theta, \phi)$ presented above. In the text following, some data-cuts will be introduced which produce different $f_i(\theta, \phi)$ distributions and corresponding $g_i(k)$ histograms.  

#### To quantitatively describe the correlation using $g(k)$-method, we calculate the RMS of the distribution as the standard error and the area representing the total amount of events that pass the cuts. 
"""

# ╔═╡ 243cecbc-26f1-11ee-0d3c-d76879f6822f
difs0 = tree.thetaEscaped .- tree.thetaEmitted    # array of ϕ - θ
rms0  = round(get_rms(difs0), digits = 2)
h0    = StatsBase.fit(Histogram, difs0, -180:1:180)
nEvents = @sprintf "%.2E" nrow(tree)

hk = histogram(difs0, 
    xlims  = (-180,180), 
    xticks = (-180:60:180),
    ylims  = (0, 1.1*maximum(h0.weights)), 
    nbins  = length(h0.edges[1]),
    xlabel = "k [°]",
    ylabel = "counts",
    label  = "k = φ - ϑ",
    lw     = 4,
    fillrange = 0,
    fillalpha = 0.2,
)

plot!(1, label = "RMS = $rms0 °", lw = 0, c =:white)
plot!(1, label = "nEvents = $nEvents", lw = 0, c =:white)


# ╔═╡ 243cecd2-26f1-11ee-27f0-d3a8869360f6
savefig("Figs/k_full_2Dist.pdf")

# ╔═╡ 243cece6-26f1-11ee-0b45-19ccd54c8e04
md"The figure above depicts the $g(k)$ histogram calculated from $f(\theta, \phi)$. The $g(k)$ distribution is centered around $k = 0$ line, which represents perfect correlation. However, the distribution is quite wide as the $RMS = 46.99^{\circ}$. All events from the original distribution are represented. "

# ╔═╡ 243ced0e-26f1-11ee-1597-d54ff9c0f831
md"""
$g(k)$ - analysis of energy cuts
===

The goal of calculating $g(k)$ and its $RMS$ is to evaluate different data-cuts. In the analysis below, various energy cuts are tested. Six data-cuts are presented - the sum of electron energies $E_{sum}$ must fulfill $E_{sum} \in (500*i, 500*i + \Delta E); \Delta E = 500 keV; i = (1, 2, 3, 4, 5, 6) keV$. (Events with $E_{sum} \in (0, 500) keV$ are omitted as there is not enough statistics). 

First, a set of $f\_i(\theta, \phi)$ is presented.  
"""

# ╔═╡ 243ced22-26f1-11ee-36cf-138b1cab41af
hs = []
for e in 500:500:3000
    e2 = e+500
    gdf = @chain tree begin
        @subset((e .<= :ESum .<= e2))
        @select(:thetaEscaped, :thetaEmitted, :ESum)
    end
    h= histogram2d(gdf.thetaEmitted, gdf.thetaEscaped, 
                 nbins = nBins, lims=(0, 180), xlabel="ϑ [°]", ylabel="φ [°]",
                 c = :coolwarm,  legend =:none, aspect_ratio=0.8,
                 ticks = (0:60:180), margin= 2mm
        )
    plot!(h, 1, label =L"E_{sum}\in (%$e, %$e2) \mathrm{keV}", lw = 0, c =:white, legend = :topleft)

    push!(hs,h)
end

# ╔═╡ 243ced36-26f1-11ee-00a6-8bf814c97fe5
plot(hs[1], hs[2], hs[3], hs[4], hs[5], hs[6], size = (1400, 1800), plot_titlefontsize = 20,plot_title = "Effects of Energy cuts on f(ϑ, φ)", layout = @layout [a b; c d; e f ])

# ╔═╡ 243ced4a-26f1-11ee-1180-1f267dff3cb6
savefig("Figs/Ecuts.pdf")

# ╔═╡ 243ced60-26f1-11ee-2874-c162cd20607f
md"Now the corresponding $g(k)$"

# ╔═╡ 243ced86-26f1-11ee-00b4-3994b80dd1fc
p1 = plot(size = (800, 800), legend=:topright, xlims=(-180, 180), xlabel="k [°]", title="unnormalized", lw = 4)
p3 = plot()
for (i,e) in enumerate(500:500:3000)
    e2 = e+500
    gdf = @chain tree begin
        @subset((e .<= :ESum .<= e2))
        @select(:thetaEscaped, :thetaEmitted, :ESum)
    end

    difs   = gdf.thetaEscaped .- gdf.thetaEmitted
    rms    = round(get_rms(difs), digits = 2)
    h1     = StatsBase.fit(Histogram, difs, -180:180)
    nEvents = @sprintf "%.2E" nrow(gdf)
    
    
    stephist!(p1, difs, 
        nbins  = Int(180/dϕ*2+1), 
        lw     = 4, 
        label  = "",#"\nE ∈ ($e, $e2)keV \nrms = $rms °\n",
        legend =:outertopright, 
        xlims  = (-180,180), 
        xlabel = "k [°]",
        xticks = -180:60:180,
        ylims  = (0, 8e4),
    )
    plot!(p3,1, label = "\nE ∈ ($e, $e2)keV \nrms = $rms °\nevents = $nEvents",legend_column =3, legendfontsize = 14, framestyle =:none, legend =:topright)
    
end
# plot!(p, 1, legend_column =3, label ="", size = (1200, 1000), legendfontsize = 14)
p1

# ╔═╡ 243cedae-26f1-11ee-39a9-35106f739dd6
md"""
It is visible from both figures $f_i(\theta, \phi)$ that the number of events increases with increasing energy. 

The $RMS$ is represented in the legend. Again, with increasing energy data-cut $RMS$ decreases. To view the width of each distribution, the histograms are normalized to area of 1.
"""

# ╔═╡ 243cedcc-26f1-11ee-1c61-5d15e60f7582
p2 = plot(size = (800, 800), legend=:topright, xlims=(-180, 180), xlabel="k [°]", title="normalized", lw = 4)
for (i,e) in enumerate(500:500:3000)
    e2 = e+500
    gdf = @chain tree begin
        @subset((e .<= :ESum .<= e2))
        @select(:thetaEscaped, :thetaEmitted, :ESum)
    end

    difs   = gdf.thetaEscaped .- gdf.thetaEmitted
    rms    = round(get_rms(difs), digits = 2)
    h1     = StatsBase.fit(Histogram, difs, -180:180)
    
    stephist!(
            p2, difs, 
            nbins  = Int(180/dϕ*2+1), 
            lw     = 4, 
            label  = "", #"\nE ∈ ($e, $e2)keV \nrms = $rms °\n",
            legend =:outertop, 
            xlims  = (-180,180), 
            xticks = -180:60:180,
            xlabel = "k [°]",
            ylims  = (0, 0.015),
            norm =:pdf
    )

end
# plot!(p2, 1, legend_column =3, label ="", size = (1200, 1000), legendfontsize = 14)
p2

# ╔═╡ 243cedea-26f1-11ee-260d-f5c15b96e0c9
plot( p1, p2, p3, size = (1400, 1000), bottom_margin = -4mm,plot_title = "k-lines distribution for various energy cuts",layout =@layout [b c;_{0.9w} a{0.1h} _] )

# ╔═╡ 243cedf2-26f1-11ee-3209-19415bb8f4e7
savefig("Figs/klines_Ecuts.pdf")

# ╔═╡ 243cee08-26f1-11ee-1b38-e3f1a1ad55e4
md"### From the figure, it is visible that increasing the energy cut improves the corellation. The $g(k)$-method can be used for quantification of correlation."

# ╔═╡ 243cee1c-26f1-11ee-3386-9ff4c615b4ec
md"""
Slicing horizontally 
===
"""

# ╔═╡ 243cee3a-26f1-11ee-3860-97fb6ce779d9
md"We can see that while applying an energy cut on the data results in decreased statistics, it did provide for a better reconstruction precision. We thus have a tool for comparing the effects of data cuts on the data."

# ╔═╡ 243cee4e-26f1-11ee-1fdf-31a66321e8ba
md"Next we look more in detail at individual $\Delta\phi$ slices. We will slice up $f(\theta,\phi)$ horizontally in slices of $\Delta\phi$ = $1^{\circ}$. However, for better visualisation of what is happening we first show a horizontal slice with $\phi \in (10, 15)\deg$ and its $g(k)$ as:"

# ╔═╡ 243cee80-26f1-11ee-18ed-f9104a3c3951
gdf = @chain tree begin
    @subset((10 .<= :thetaEscaped .<= 15))
    @select(:thetaEscaped, :thetaEmitted, :weights)
end

fh2d = Hist2D(                                          
(gdf[!,2], gdf[!,1]),
    Weights(gdf.weights),
(minAngle:dEmitted:maxAngle, minAngle:dEmitted:maxAngle), 
) 

gs4 = get_diagonal_sums(fh2d)
ks4 = get_k_factors(fh2d);

h = histogram2d(gdf[!,2], gdf[!,1];
    nbins        = (nBins, nBins),
    weights      = gdf.weights,
    xlabel       = "θemitted -> θ",
    ylabel       = "θescaped -> ϕ",
    legend       = :topright,
    title        = string("f(ϕ, θ): dϕ ∈ (10, 15)°, ", nrow(gdf), " entries"),
    lims         = (0, 180),
    c = :coolwarm
    # aspect_ratio = 1,
    )
h1d = stephist(gdf.thetaEmitted, nbins = nBins, weights = gdf.weights, 
                label ="ϕ ∈ (10,15)",xlabel ="θ", ylabel ="counts",
                title = "1d Histogram", legend =:topright)
vspan!([10,15], label ="(10-15)° marker", c =:black, alpha= 0.2)
p = plot(ks4, gs4, label = "", xlabel = "k-factor", ylabel ="g(k)")

l = @layout [a{0.5w} b; c{0.5w} _]
plot(h,p,h1d, layout = l, plot_title = "Horizontal slice of f(ϕ,θ)", size = (1000,800))

# ╔═╡ 243cee94-26f1-11ee-04e2-67bad5b72bc0
md"This procedure is repeated for each slice with $\Delta\phi = 1^{\circ}$. Slicing $f(\phi, \theta)$ horizontally to cover the whole 0 - 180 degree range yields $g(k)$s:"

# ╔═╡ 243ceebc-26f1-11ee-02b9-91fe10bf7bf5
p = plot(size = (800, 800), legend=:none, xlims=(-180, 180), xlabel="k [°]", ylabel="g(k) [#/°]", lw = 4, dpi =3)
for (i,p) in enumerate(0:dϕ:180)
    p2 = p+dϕ
    gdf = @chain tree begin
        @subset((p .<= :thetaEscaped .<= p2))
        @select(:thetaEscaped, :thetaEmitted)
    end

    difs   = gdf.thetaEscaped .- gdf.thetaEmitted
    rms    = round(get_rms(difs), digits = 2)
    h1     = StatsBase.fit(Histogram, difs, -180:180)
    
    stephist!(p, difs, 
        nbins  = Int(180/dϕ*2+1), 
        lw     = 4, 
        label  ="",
        legend =:topright, 
        xlims  = (-180,180), 
        xlabel = "k [°]",
        ylims  = (0, 1e3)
    )

end
p;

# ╔═╡ 243ceed0-26f1-11ee-23a2-f93c2573a36b
md"Not much can be deduced in this example. So many lines are difficult to decipher. However, if one were to look at the individual $\Delta\phi$ cuts as a new dimension, we can look at the graph in the plane of $(k, \Delta\phi)$ with z-direction being the value of $g(k, \Delta\phi)$. "

# ╔═╡ 243ceeee-26f1-11ee-0167-690c3953dbd9
dfGOrig = get_gs_df(tree, dϕ, sign)
matGOrig = df_to_mat(dfGOrig);

# ╔═╡ 243cef16-26f1-11ee-2f87-bfd5a7df2204
xRange = dϕ-180:dϕ:180-dϕ
yRange = 0:dϕ:180-dϕ
sf1 = surface(xRange, yRange,  matGOrig, legend =:none, xlabel ="k", ylabel ="Δφ", zlabel="counts", tickfontsize = 12, c= :turbo, yticks = (0:60:180), xticks = (-180:60:180))
hm1 = heatmap(xRange, yRange .+ dϕ/2,  matGOrig, ylabel ="Δφ", xlabel ="k" , tickfontsize = 14, bottom_margin  = 6Plots.mm, c= :turbo, yticks = (0:60:180), xticks = (-180:60:180), ylims = (0,180), xlims = (-180,180), title = "Partition vs k distribution", xrotation = 45)
vline!([0], label ="", c = :black, lw = 2, s=:dash)
plot(sf1,hm1, size =(1000,400), layout = @layout [a{0.4w} b])

# ╔═╡ 243cef2c-26f1-11ee-22e7-b954d192ba52
savefig("Figs/3D_klines.pdf")

# ╔═╡ 243cef3e-26f1-11ee-3d02-adee0624b712
md"Now we can see a few important features. First of all, there are two peaks visible in the left figure, with the higher peak (more statistics) being in the region of $130^{\circ} < \phi < 17^{\circ}0$ . Secondly,  we can see the deviation of the peaks from the `` k = 0`` line in the right figure. There are two hotspots visible. First hotspot (corresponding to the lower peak in  figure) is centered around $\Delta\phi \approx 30^{\circ}$ and is shifted slightly to the right of the ``k = 0`` line. The escaped angle overestimates the emitted angle. Second hotspot (corresponding to the higher peak in figure) is centered around $\Delta\phi \approx 150^{\circ}$ and is shifted visibly to the left of the ``k = 0`` line. The escaped angle underestimates the emitted angle. Lastly, we can see that the regions $\phi \approx 0^{\circ}$ and $\phi \approx 180^{\circ}$ are squeezed toward higher, lower angles, respectively. "

# ╔═╡ 243cef52-26f1-11ee-2467-119f6f3a61c5
md"Furthermore, we can also look at how the individual $\Delta\phi$ slices look in terms of statistical variables (mean, mode, median). For each variable, obtained from the $\phi(\theta)$ distributions. "

# ╔═╡ 243cef7a-26f1-11ee-0b31-7744f0e87bc3
means   = Vector{Float64}(undef, Int(180/dϕ))
modes   = Vector{Float64}(undef, Int(180/dϕ))
medians = Vector{Float64}(undef, Int(180/dϕ))

for (i, ϕ) in enumerate(1:dϕ:180)
    cutEdges1 = get_cut_edges(ϕ - 1, 1, dϕ, "p")                    # provides the lower and upper cut 
    sdf       = @chain tree begin                                   # filter out the dataframe
        @select(:thetaEscaped, :thetaEmitted, :weights)                       # keeps only the two angles columns
        @subset((cutEdges1[1] .<= :thetaEscaped .<= cutEdges1[2]))  # keeps only rows where ϕ is within the cut edges
    end
    
    stats = get_slice_stats(
                            cutEdges1[1],  # ϕmin
                            cutEdges1[2],  # ϕmax
                            0,             # Emin
                            3500,          # Emax
                            sdf.thetaEmitted, 
                            dϕ,
                            Weights(sdf.weights)
                            )
    
    means[i] = stats[5] 
    modes[i] = stats[6] .+ dϕ/2
    medians[i] =  stats[7] 
end

# ╔═╡ 243cef8c-26f1-11ee-121d-7f87ac0992d7
scatter( means, xPts, ms = 4, label = "means" , xlabel = "mean g_i(k)", ylabel ="Δϕ_i(k)", aspect_ratio =1, legend = :topleft, lims = (0,180))
scatter!( modes, xPts, ms = 3, label ="modes" )
scatter!( medians, xPts, ms = 3, label ="medians" )

# ╔═╡ 243cefac-26f1-11ee-16ad-1336fe5e865d
md"""
Shift to reduce RMS
===

"""

# ╔═╡ 243cefca-26f1-11ee-250e-2d34a6a17ca3
md"""
The goal of the analysis which is presented in the pages below is to alter the $\phi$ data so that the overall RMS is reduced - in other words, we want to find such representation of $\phi$ which leads to least error. 

First, we define (or describe) three variables:
1. $\phi$ - the escape angle is an experimentally measurable variable,
2. $\theta$ - the decay angle is inaccesible through experiment,
3. $\phi'$ - a new angle which we define as the *representation* of measured $\phi$ (we **want** $\phi' \approx \theta$).

Since $\theta$ is no accessible, we have introduced a new variable $\phi'$ which is supposed to represent the **most likely $\theta$ which the measured $\phi$ originated from**. We could see in the individual slice histograms that each event in $\Delta\phi$ slice can have originated from any $\theta$, with varying probabilites. We want $\phi'$ to represent the most likely $\theta$. 

To obtain $\phi'(\phi)$ we define it as:
$\phi'(\phi) = \phi + s$
Where $s$ is a constant *shift* which minimizes RMS for each $\phi$ obtained from $\Delta\phi$ slices. 

Here we introduce three various possibilities for what $s$ could be. The most obvious shift to use would be to take advantage of statistical estimators presented earlier. We present *shift* for each $\phi \in \Delta\phi$ so that when we apply the shift, the given statistical estimator will lie within the $\Delta\phi$ slice. 

For example, for mean, to obtain $s$ we use:
$s = \bar{\theta} - \Delta\phi_{center}$. Where $\bar{\theta}$ is the mean of the $\theta$ distribution of the given $\Delta\phi$ slice and $\Delta\phi_{center}$ is the bincenter of the slice. Analogous for mode and median.
"""

# ╔═╡ 243cefe8-26f1-11ee-3b15-05c62a9e56ef
res_means = [ means[i] - y(i*dϕ) for i in 1:length(means) ]
res_modes = [ modes[i] - y(i*dϕ) for i in 1:length(modes) ]
res_medians = [ medians[i] - y(i*dϕ) for i in 1:length(medians) ]

scatter( res_means, xPts, ms = 3, label = "mean" , xlabel = "shift", ylabel ="Δϕ_i(k)", ylims = (0,180), aspect_ratio =1, legend = :topright)
scatter!( res_modes, xPts, ms = 3, label ="mode" )
scatter!( res_medians, xPts, ms = 3, label ="median" )
vline!([0], label ="", c = :black, lw = 3, s=:dash)

# ╔═╡ 243ceffc-26f1-11ee-08a9-a92001550db8
md"We can see that for each $\Delta\phi_i$ slice the three estimators provide different values to shift the angles by. The most drastic shift (ie. farthest away from ``k=0``) is given by mean, the least on the other hand by mode. "

# ╔═╡ 243cf010-26f1-11ee-1d15-8fc9416287e8
md"To avoid undesirable discretization of our data, we fit the shifts. We also flipped the axes so that we get $s(\Delta\phi_i)$."

# ╔═╡ 243cf042-26f1-11ee-07d9-0d5a5145d62f
scatter(  xPts,res_means, ms = 3, label = "shift by" , ylabel = "shift", xlabel ="Δϕ_i(k)", xlims = (0,180), aspect_ratio =1, legend = :topright)
f1 = Polynomials.fit(xPts,  res_means, 6 )

plot!(f1, extrema(xPts)..., label="Fit")

# ╔═╡ 243cf062-26f1-11ee-043f-8d2e1198c20a
md"""

Now we shift each $\phi$ in the original data set to obtain a new set of $\phi'$, we do so by $\phi' = \phi +s$. 
"""

# ╔═╡ 243cf074-26f1-11ee-388d-c9642d413e88
md"""
*Shift by mean*
====
"""

# ╔═╡ 243cf094-26f1-11ee-3c06-db94de90db7b
modTree2 = @chain tree begin
    @select(:thetaEmitted, :thetaEscaped, :weights)
    @rtransform :bin =  get_bin_center(:thetaEscaped, Int(180/dϕ))  # create a vector of bin centers (which bin ϕ falls inside)
    @transform :thetaEscapedOld = :thetaEscaped                     # create a copy of the old ϕ (for comparison only)

    @rtransform :thetaEscapedDisc = :thetaEscapedOld + res_means[Int(ceil(:bin/dϕ))] # shift ϕ by s: ϕ' = ϕ + s 
    @rtransform :thetaEscaped = :thetaEscapedOld + f1(:thetaEscapedOld) # shift ϕ by s: ϕ' = ϕ + s 
    @subset( 0 .< :thetaEscaped .< 180) # keep only physical angles
end

# ╔═╡ 243cf09c-26f1-11ee-396e-f76f88359683
dfG2 = get_gs_df(modTree2, dϕ, sign)
matG2 = df_to_mat(dfG2);

# ╔═╡ 243cf0b0-26f1-11ee-2cba-75a34b6d43e5
md"We look at the $f(\theta, \phi')$ figure."

# ╔═╡ 243cf0ce-26f1-11ee-18f1-d756d8d021ec
h2d2 = histogram2d(modTree2.thetaEmitted, modTree2.thetaEscaped,
    nbins        = (nBins, nBins),
    weights      = modTree2.weights,
    xlabel       = "θ",
    ylabel       = "ϕ'",
    legend       = :topright,
    title        = string("f(θ, ϕ'), ", nrow(modTree2), " entries"),
    lims         = (0, 180),
    aspect_ratio = 1,
    right_margin = 6Plots.mm,
    c= :turbo
)
plot!(xPts, xPts, label ="", c= :black, style= :dash, lw =3)

# ╔═╡ 243cf0e2-26f1-11ee-1295-eb9d5613b7a5
md"We can see that shifting by mean value resulted in *squeezing* the phase-space. We have reduced the range of angles which we can interpret in our measuremt. However, this should lead toward reduces RMS. We look at that in the following figures."

# ╔═╡ 243cf100-26f1-11ee-0ea4-e376ab7a28a7
md"First, for comparison the original dataset with $f(\theta, \phi)$, $g(k, \Delta\phi)$, calculated $RMS$ for each $g_i(k)$ and total $RMS$.  "

# ╔═╡ 243cf13c-26f1-11ee-30fd-ffd76c8c0b79
h2d1 = histogram2d(
    tree.thetaEmitted,
    tree.thetaEscaped;
    nbins        = (nBins, nBins),
    weights      = tree.weights,
    xlabel       = "θ",
    ylabel       = "φ",
    legend       = :topright,
    title        = string("f(θ, φ)"),
    lims         = (0, 180),
#     aspect_ratio = 1,
    c = :turbo,
)
plot!(xPts, xPts, label ="", c= :black, style= :dash, lw =3)
rms1 = [ get_rms(dfGOrig[:,i], dfGOrig[:,1]) for i in 2:ncol(dfGOrig) ]
rmsTotalUnModded = round(get_rms(gs1, ks1 .* dEmitted ), digits = 2)

sct1 = scatter( xPts, rms1, ms=4, legend=:top, xlabel ="Δφ", ylabel ="RMS", c= :red, label ="original RMS = $rmsTotalUnModded", title = "RMS per partition", xlims = (0,180), xticks = 0:60:180 )

fh2d1 = Hist2D((tree.thetaEmitted,tree.thetaEscaped), 
    Weights(tree.weights),
    (minAngle:dEmitted:maxAngle, minAngle:dEmitted:maxAngle))
gs1 = get_diagonal_sums(fh2d1)
ks1 = get_k_factors(fh2d1);

gk1 = plot(ks1 .* dEmitted, gs1, legend=:topright, xlims=(-180, 180), xlabel="k", ylabel="count", label="")

# plot(title, h2d1, hm1, sct1, layout = @layout[a{0.05h};b c; d _] , size = (1100, 800))
plot(h2d1, hm1, sct1, gk1, layout = @layout[a b; c d] , size = (1200, 800), plot_title= "Unmodified angles")

# ╔═╡ 243cf146-26f1-11ee-0da6-9dc072a83446
column_unmodded = plot(hm1, sct1, layout = grid(1,2), size = (1200, 600), bottom_margin = 12mm)

# ╔═╡ 243cf15a-26f1-11ee-2bff-e7cdb192c94a
savefig("Figs/2D_klines_and_RMS_unmodded.pdf")

# ╔═╡ 243cf16e-26f1-11ee-3af3-75f9696d1199
savefig(sct1, "Figs/RMS_unmodded.pdf")

# ╔═╡ 243cf182-26f1-11ee-2707-8bbc39b2651c
md"And the **modified** dataset."

# ╔═╡ 243cf1ca-26f1-11ee-2839-037bdc7f5d72
hm2 = Plots.heatmap(xRange, yRange,  matG2, ylabel ="Δϕ'", xlabel ="k", c =:turbo )
vline!([0], label ="", c = :black, lw = 3, s=:dash)

rms2 = [ get_rms(dfG2[:,i], dfG2[:,1]) for i in 2:ncol(dfG2) ]
yMin = 0.9*minimum(filter(x -> x .> 0, rms2)) # get the minimum rms value, excluding 0
yMax = 1.1*maximum(filter(x -> x .> 0, rms2))

sct2 = scatter( xPts, rms2, label ="rms(Δϕ')", ms=2, 
                legend=:top, xlabel ="Δϕ' slice", c =:blue,
                ylabel ="rms", ylims = (yMin, yMax) )

fh2d2 = Hist2D((modTree2.thetaEmitted, modTree2.thetaEscaped),
    Weights(modTree2.weights),
    (minAngle:dEmitted:maxAngle, minAngle:dEmitted:maxAngle))
gs2 = get_diagonal_sums(fh2d2)
ks2 = get_k_factors(fh2d2);

gk2 = plot(ks2 .* dEmitted, gs2, legend=:topright, xlims=(-179, 179), 
            xlabel="k-factor", ylabel="g(k)", label="g_2(k)")



plot(h2d2, hm2, sct2, gk2, layout = @layout[a b; c d] , size = (1200, 800), 
    plot_title= "Modified by mean angles")

# ╔═╡ 243cf1dc-26f1-11ee-20c4-9b8f378fab60
rmsTotalUnModded = round(get_rms(gs1, ks1 .* dEmitted ), digits = 2)

# ╔═╡ 243cf1f0-26f1-11ee-11a2-6768bc542c93
rmsTotalModded = round(get_rms(gs2, ks2 .* dEmitted ), digits = 2)

# ╔═╡ 243cf204-26f1-11ee-17a1-d33c9b43bb7e
md"Now we compare the four $RMS$ figures together. "

# ╔═╡ 243cf222-26f1-11ee-0498-3f3d35575cb7
scComp2 = scatter( [xPts xPts], [rms1 rms2],  label =["rms_1 = $rmsTotalUnModded °" "rms_2 = $rmsTotalModded °"], 
        ms=3, c = [:red :blue], bottom_margin = 8Plots.mm,
        legend=:topleft, xlabel ="Δϕ_i slice", ylabel ="rms", ylims = (yMin, 1.1*maximum(rms1)))

gkComp2 = plot( [ks1 .* dEmitted, ks2 .* dEmitted], [gs1, gs2], label = ["g_1(k)" "g_2(k)"], legend = :topleft,
                xlabel = "k", ylabel = "g(k)", seriestype = :stepmid, lw = 4, c = [:red :blue])

plot(scComp2, gkComp2, size = (1200, 600))

# ╔═╡ 243cf22a-26f1-11ee-144c-8326ffd944ff
md"We can see that shifting by mean value results in reduced $RMS$, the goal is achieved. "

# ╔═╡ 243cf24a-26f1-11ee-1db7-2f0f06fb55fd
md"## Finally, we can now provide a function, which as an input takes the measured angle $\phi$ and as an output provides $\phi'$ (the most likely $\theta$): **$\phi'(\phi)$**."

# ╔═╡ 243cf268-26f1-11ee-20c0-a744732e3706
plot(
    xPts, xPts .+ f1.(xPts), lims = (0,180), lw = 4, c= :blue, 
    label ="", xlabel = "ϕ", ylabel ="ϕ'", aspect_ratio = 1, 
    title = "ϕ'(ϕ); measured (ϕ) vs reported (ϕ') angle"
)


# ╔═╡ 243cf27c-26f1-11ee-1b4b-39eef8360fa2


# ╔═╡ Cell order:
# ╠═2437222a-26f1-11ee-30b4-a96159720d68
# ╠═243c0132-26f1-11ee-2333-7ba24e313a18
# ╟─243ce74e-26f1-11ee-3178-1bf57ca82af4
# ╠═243ce824-26f1-11ee-10d6-d1044295b35e
# ╠═243ce85e-26f1-11ee-3a2c-6f32f20a1cbe
# ╠═243ce872-26f1-11ee-37fe-d1e78aba7c33
# ╠═243ce884-26f1-11ee-32c6-532ac1a54626
# ╠═243ce89a-26f1-11ee-108d-4b137230b0a0
# ╟─243ce8c2-26f1-11ee-236d-e9bd0a43478c
# ╠═243ce8d6-26f1-11ee-26e1-b5225e40b167
# ╟─243ce8fe-26f1-11ee-24f4-bfcdc4137274
# ╠═243ce930-26f1-11ee-31cd-f17ad3702d2a
# ╠═243ce93a-26f1-11ee-3957-976266e37bae
# ╟─243ce95a-26f1-11ee-0468-c92a72f82437
# ╠═243ce980-26f1-11ee-2d0c-851fdaae8542
# ╟─243ce9a8-26f1-11ee-2b76-6bba0d7bf2ab
# ╟─243ce9c6-26f1-11ee-35e9-852f3fdeab95
# ╠═243ce9e4-26f1-11ee-0f12-111fbfb0e83c
# ╠═243ce9f8-26f1-11ee-30ef-75c9c32de503
# ╟─243cea0c-26f1-11ee-14b4-b3ead850f465
# ╟─243cea1e-26f1-11ee-1a53-4fddd3aeafaf
# ╠═243cea7a-26f1-11ee-1a5d-d7476036e22c
# ╠═243cea90-26f1-11ee-3bc4-c51afe379f0a
# ╟─243ceaac-26f1-11ee-0a5a-b36800325949
# ╠═243ceac2-26f1-11ee-2725-9128a5488720
# ╠═243ceafc-26f1-11ee-07d8-cffedeccbc13
# ╠═243ceb38-26f1-11ee-1973-052f01ea98b7
# ╠═243ceb6a-26f1-11ee-2883-69c9f62c651d
# ╠═243ceb7e-26f1-11ee-0789-9530a207fa5f
# ╠═243ceba6-26f1-11ee-1260-65be19b68358
# ╠═243cebba-26f1-11ee-26e7-3b4cb7fdc77c
# ╠═243cebec-26f1-11ee-345d-950cb82da77d
# ╠═243cec00-26f1-11ee-24e7-bf8fe02ade17
# ╟─243cec14-26f1-11ee-0d28-033b7e625dd2
# ╟─243cec2a-26f1-11ee-16e8-25f63d7cb2dc
# ╠═243cec50-26f1-11ee-2271-0596a4eec19c
# ╠═243cec64-26f1-11ee-28e4-798cb7b12f20
# ╟─243cec78-26f1-11ee-196b-cbd9a31cdf4c
# ╟─243cec8a-26f1-11ee-387e-6f6064f3b56f
# ╠═243cecbc-26f1-11ee-0d3c-d76879f6822f
# ╠═243cecd2-26f1-11ee-27f0-d3a8869360f6
# ╟─243cece6-26f1-11ee-0b45-19ccd54c8e04
# ╟─243ced0e-26f1-11ee-1597-d54ff9c0f831
# ╠═243ced22-26f1-11ee-36cf-138b1cab41af
# ╠═243ced36-26f1-11ee-00a6-8bf814c97fe5
# ╠═243ced4a-26f1-11ee-1180-1f267dff3cb6
# ╟─243ced60-26f1-11ee-2874-c162cd20607f
# ╠═243ced86-26f1-11ee-00b4-3994b80dd1fc
# ╟─243cedae-26f1-11ee-39a9-35106f739dd6
# ╠═243cedcc-26f1-11ee-1c61-5d15e60f7582
# ╠═243cedea-26f1-11ee-260d-f5c15b96e0c9
# ╠═243cedf2-26f1-11ee-3209-19415bb8f4e7
# ╟─243cee08-26f1-11ee-1b38-e3f1a1ad55e4
# ╟─243cee1c-26f1-11ee-3386-9ff4c615b4ec
# ╟─243cee3a-26f1-11ee-3860-97fb6ce779d9
# ╟─243cee4e-26f1-11ee-1fdf-31a66321e8ba
# ╠═243cee80-26f1-11ee-18ed-f9104a3c3951
# ╟─243cee94-26f1-11ee-04e2-67bad5b72bc0
# ╠═243ceebc-26f1-11ee-02b9-91fe10bf7bf5
# ╟─243ceed0-26f1-11ee-23a2-f93c2573a36b
# ╠═243ceeee-26f1-11ee-0167-690c3953dbd9
# ╠═243cef16-26f1-11ee-2f87-bfd5a7df2204
# ╠═243cef2c-26f1-11ee-22e7-b954d192ba52
# ╟─243cef3e-26f1-11ee-3d02-adee0624b712
# ╟─243cef52-26f1-11ee-2467-119f6f3a61c5
# ╠═243cef7a-26f1-11ee-0b31-7744f0e87bc3
# ╠═243cef8c-26f1-11ee-121d-7f87ac0992d7
# ╟─243cefac-26f1-11ee-16ad-1336fe5e865d
# ╟─243cefca-26f1-11ee-250e-2d34a6a17ca3
# ╠═243cefe8-26f1-11ee-3b15-05c62a9e56ef
# ╟─243ceffc-26f1-11ee-08a9-a92001550db8
# ╟─243cf010-26f1-11ee-1d15-8fc9416287e8
# ╠═243cf042-26f1-11ee-07d9-0d5a5145d62f
# ╟─243cf062-26f1-11ee-043f-8d2e1198c20a
# ╟─243cf074-26f1-11ee-388d-c9642d413e88
# ╠═243cf094-26f1-11ee-3c06-db94de90db7b
# ╠═243cf09c-26f1-11ee-396e-f76f88359683
# ╟─243cf0b0-26f1-11ee-2cba-75a34b6d43e5
# ╠═243cf0ce-26f1-11ee-18f1-d756d8d021ec
# ╟─243cf0e2-26f1-11ee-1295-eb9d5613b7a5
# ╟─243cf100-26f1-11ee-0ea4-e376ab7a28a7
# ╠═243cf13c-26f1-11ee-30fd-ffd76c8c0b79
# ╠═243cf146-26f1-11ee-0da6-9dc072a83446
# ╠═243cf15a-26f1-11ee-2bff-e7cdb192c94a
# ╠═243cf16e-26f1-11ee-3af3-75f9696d1199
# ╟─243cf182-26f1-11ee-2707-8bbc39b2651c
# ╠═243cf1ca-26f1-11ee-2839-037bdc7f5d72
# ╠═243cf1dc-26f1-11ee-20c4-9b8f378fab60
# ╠═243cf1f0-26f1-11ee-11a2-6768bc542c93
# ╟─243cf204-26f1-11ee-17a1-d33c9b43bb7e
# ╠═243cf222-26f1-11ee-0498-3f3d35575cb7
# ╟─243cf22a-26f1-11ee-144c-8326ffd944ff
# ╟─243cf24a-26f1-11ee-1db7-2f0f06fb55fd
# ╠═243cf268-26f1-11ee-20c0-a744732e3706
# ╠═243cf27c-26f1-11ee-1b4b-39eef8360fa2
