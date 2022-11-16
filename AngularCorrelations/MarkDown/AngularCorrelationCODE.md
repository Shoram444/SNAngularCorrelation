```julia
using Revise
using StatsPlots, UnROOT, StatsBase, Polynomials, LinearAlgebra
using FHist, MPThemes, DataFramesMeta, Distributions 
```


```julia
ENV["COLUMNS"] = 2000
ENV["LINES"] = 20
```




    20



Custom functions are placed in the MiscFuncs.jl file.


```julia
include("MiscFuncs.jl")
using .MiscFuncs
```


```julia
gr()
default(fmt = :png)
theme(:dao;
    size           = (800, 800),
    legend         = :topleft,
    guidefontsize  = 16,
    tickfontsize   = 12,
    titlefontsize  = 16,
    legendfontsize = 12,
    left_margin    = 4Plots.mm,
    right_margin   = 8Plots.mm,
    top_margin     = 4Plots.mm,
    bottom_margin  = 6Plots.mm,
    dpi            = 100,
    :colorbar_titlefontsize => 16,
);
```


```julia
baseDir = "/home/shoram/Work/PhD_Thesis/Job15/AngularCorrelations/"
```




    "/home/shoram/Work/PhD_Thesis/Job15/AngularCorrelations/"




```julia
figDir = joinpath("/media/shoram/Extra SSD/CernBox/Work/Presentations/20221006_SuperNEMO_Analysis_Meeting/Figs")
```




    "/media/shoram/Extra SSD/CernBox/Work/Presentations/20221006_SuperNEMO_Analysis_Meeting/Figs"




```julia
f = ROOTFile(
    baseDir*"AngularCorrelationAllEnergies96MilEvents.root",
);
tree = DataFrame(LazyTree(f, "tree", keys(f["tree"])));
f= nothing; # free file from memory with gc

```

### ``@transform`` adds a column ``:ESum`` to the ``tree`` which contains the sum of the electron energies


```julia
@transform! tree :ESum = :reconstructedEnergy2 + :reconstructedEnergy1;
```

### Initializing constants.


```julia

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
```




    180




```julia
colors = [palette(:seaborn_bright)[i] for i in 1:length(palette(:seaborn_bright))];
```

The 2d Histogram of $\phi$ vs $\theta$ is defined to be $f(\theta,\phi)$. For each combination of $\phi$ and $\theta$, the bin number is obtained as corresponding value of $f(\theta,\phi)$.


```julia
rho(_cosdTheta) = 0.5 - 0.5 * _cosdTheta

# @transform! tree :weights = Weights(rho.(cosd.(:thetaEmitted))) # weights by pdf
@transform! tree :weights = 1                                     # unweighted


```




<div class="data-frame"><p>15,747,968 rows × 18 columns</p><table class="data-frame"><thead><tr><th></th><th>momentumEmitted2x</th><th>reconstructedEnergy2</th><th>momentumEscaped1z</th><th>momentumEmitted1x</th><th>thetaEmitted</th><th>momentumEscaped1y</th><th>momentumEmitted1z</th><th>momentumEmitted1y</th><th>momentumEscaped2z</th><th>thetaEscaped</th><th>momentumEscaped2y</th><th>momentumEmitted2y</th><th>reconstructedEnergy1</th><th>momentumEscaped1x</th><th>momentumEmitted2z</th><th>momentumEscaped2x</th><th>ESum</th><th>weights</th></tr><tr><th></th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Int64">Int64</th></tr></thead><tbody><tr><th>1</th><td>0.723365</td><td>604.269</td><td>-0.471071</td><td>-1.41249</td><td>89.6031</td><td>2.26532</td><td>-0.946003</td><td>2.12515</td><td>-0.155021</td><td>94.6325</td><td>0.351932</td><td>0.615735</td><td>2258.02</td><td>-1.36057</td><td>0.283402</td><td>0.777892</td><td>2862.28</td><td>1</td></tr><tr><th>2</th><td>-1.87869</td><td>1651.11</td><td>0.111843</td><td>0.620354</td><td>96.965</td><td>-0.427881</td><td>0.585939</td><td>-0.857937</td><td>0.547949</td><td>139.497</td><td>-0.460427</td><td>-0.935666</td><td>802.522</td><td>1.07408</td><td>0.0929118</td><td>-1.94376</td><td>2453.63</td><td>1</td></tr><tr><th>3</th><td>-0.451873</td><td>450.772</td><td>-0.204413</td><td>0.143294</td><td>154.233</td><td>-0.344668</td><td>0.461237</td><td>-0.617821</td><td>-0.00393392</td><td>138.388</td><td>-0.0885292</td><td>0.444109</td><td>425.001</td><td>0.557768</td><td>-0.51231</td><td>-0.766414</td><td>875.773</td><td>1</td></tr><tr><th>4</th><td>0.00622859</td><td>2520.27</td><td>0.0712075</td><td>-0.27381</td><td>67.3403</td><td>0.440423</td><td>-0.0579716</td><td>-1.01549</td><td>0.518941</td><td>154.281</td><td>-2.14362</td><td>-1.34793</td><td>659.761</td><td>0.717866</td><td>2.66656</td><td>-1.67234</td><td>3180.03</td><td>1</td></tr><tr><th>5</th><td>-0.660664</td><td>834.025</td><td>-0.282475</td><td>-1.77769</td><td>11.2133</td><td>-2.1054</td><td>-0.0196651</td><td>-1.94104</td><td>0.221755</td><td>53.0724</td><td>-0.966977</td><td>-1.05076</td><td>2170.29</td><td>-1.37587</td><td>0.0859927</td><td>0.291419</td><td>3004.32</td><td>1</td></tr><tr><th>6</th><td>-0.44519</td><td>1530.67</td><td>0.917931</td><td>0.0695716</td><td>168.417</td><td>1.51856</td><td>2.09515</td><td>0.648145</td><td>-1.07627</td><td>81.0901</td><td>0.0381113</td><td>-0.472476</td><td>1741.93</td><td>-0.999117</td><td>-1.86705</td><td>-1.51857</td><td>3272.61</td><td>1</td></tr><tr><th>7</th><td>-0.445143</td><td>1784.88</td><td>0.277989</td><td>0.889556</td><td>122.632</td><td>-0.0529092</td><td>0.257764</td><td>-0.238635</td><td>0.00632242</td><td>69.6012</td><td>1.73728</td><td>1.75164</td><td>573.353</td><td>0.858823</td><td>-1.32043</td><td>0.80815</td><td>2358.23</td><td>1</td></tr><tr><th>8</th><td>-1.42443</td><td>1723.99</td><td>0.250602</td><td>-1.16971</td><td>4.96695</td><td>0.834453</td><td>0.0445409</td><td>1.26977</td><td>-0.894417</td><td>34.9839</td><td>1.23976</td><td>1.64022</td><td>1290.01</td><td>-1.4441</td><td>-0.121338</td><td>-1.4647</td><td>3014.0</td><td>1</td></tr><tr><th>9</th><td>-0.0553668</td><td>1281.87</td><td>0.457432</td><td>-0.344509</td><td>79.857</td><td>1.30726</td><td>-0.133961</td><td>1.66448</td><td>0.94093</td><td>61.1214</td><td>0.715475</td><td>0.160922</td><td>1268.96</td><td>0.500476</td><td>-1.71006</td><td>-0.744009</td><td>2550.83</td><td>1</td></tr><tr><th>10</th><td>-1.7249</td><td>2124.13</td><td>-0.525607</td><td>-0.976017</td><td>43.0942</td><td>0.765786</td><td>-0.71208</td><td>-0.308521</td><td>-2.15062</td><td>33.994</td><td>0.664738</td><td>1.23561</td><td>836.582</td><td>-0.544891</td><td>-1.47674</td><td>-1.08732</td><td>2960.71</td><td>1</td></tr><tr><th>11</th><td>-0.19787</td><td>467.234</td><td>-0.513265</td><td>-0.80288</td><td>35.1387</td><td>-0.210073</td><td>-0.13689</td><td>-0.695541</td><td>-0.230839</td><td>104.069</td><td>-0.438105</td><td>-0.806524</td><td>675.699</td><td>-0.763111</td><td>-0.0786541</td><td>0.483935</td><td>1142.93</td><td>1</td></tr><tr><th>12</th><td>-0.529524</td><td>2307.59</td><td>-0.682419</td><td>-0.307922</td><td>150.65</td><td>-0.818542</td><td>-0.238325</td><td>-1.27159</td><td>1.96584</td><td>101.923</td><td>0.305447</td><td>2.70679</td><td>913.668</td><td>0.601471</td><td>-0.27609</td><td>1.57859</td><td>3221.26</td><td>1</td></tr><tr><th>13</th><td>1.22009</td><td>976.306</td><td>-0.64402</td><td>0.895972</td><td>56.3877</td><td>0.122387</td><td>-0.727497</td><td>0.989891</td><td>0.620286</td><td>55.6631</td><td>0.386989</td><td>0.454814</td><td>1093.07</td><td>1.29631</td><td>0.505461</td><td>1.11372</td><td>2069.37</td><td>1</td></tr><tr><th>14</th><td>1.62335</td><td>1841.63</td><td>-0.772472</td><td>0.148038</td><td>128.995</td><td>-1.01428</td><td>-0.825084</td><td>-1.266</td><td>0.520843</td><td>94.8403</td><td>1.28805</td><td>1.23633</td><td>1091.04</td><td>0.781192</td><td>1.05356</td><td>1.81786</td><td>2932.68</td><td>1</td></tr><tr><th>15</th><td>-0.975718</td><td>619.331</td><td>-1.0753</td><td>-0.587855</td><td>63.2692</td><td>1.86365</td><td>-1.97465</td><td>1.62362</td><td>-0.608049</td><td>40.9504</td><td>0.0886557</td><td>0.0910889</td><td>2161.47</td><td>-1.43294</td><td>-0.23707</td><td>-0.75854</td><td>2780.8</td><td>1</td></tr><tr><th>16</th><td>-0.181153</td><td>985.07</td><td>1.34426</td><td>-1.8497</td><td>40.4733</td><td>-0.640471</td><td>1.69434</td><td>-0.0995104</td><td>0.984389</td><td>55.1476</td><td>0.757152</td><td>0.0829501</td><td>2050.88</td><td>-1.97838</td><td>1.39191</td><td>-0.54679</td><td>3035.95</td><td>1</td></tr><tr><th>17</th><td>1.62458</td><td>2530.76</td><td>0.0707918</td><td>-1.20434</td><td>124.468</td><td>0.512761</td><td>0.35071</td><td>0.34492</td><td>1.19175</td><td>142.771</td><td>-2.09416</td><td>-2.11764</td><td>886.684</td><td>-1.14899</td><td>1.36657</td><td>1.72853</td><td>3417.44</td><td>1</td></tr><tr><th>18</th><td>0.685257</td><td>639.204</td><td>0.0196268</td><td>-0.00221343</td><td>128.263</td><td>0.374426</td><td>0.395698</td><td>1.14861</td><td>-0.17065</td><td>20.8424</td><td>0.0968995</td><td>-0.742877</td><td>806.954</td><td>0.864823</td><td>0.20101</td><td>0.966251</td><td>1446.16</td><td>1</td></tr><tr><th>19</th><td>-0.211174</td><td>628.139</td><td>0.502857</td><td>-2.74734</td><td>105.706</td><td>2.36161</td><td>-0.622121</td><td>1.39787</td><td>-0.258807</td><td>172.673</td><td>-0.709862</td><td>-0.736769</td><td>2674.92</td><td>-1.86425</td><td>0.670143</td><td>0.514478</td><td>3303.06</td><td>1</td></tr><tr><th>20</th><td>0.0476415</td><td>707.158</td><td>-0.736812</td><td>-1.20521</td><td>142.67</td><td>1.44929</td><td>-0.496897</td><td>1.82196</td><td>0.813798</td><td>112.558</td><td>0.299922</td><td>-0.860809</td><td>1786.85</td><td>-1.10866</td><td>0.692477</td><td>0.550832</td><td>2494.01</td><td>1</td></tr><tr><th>&vellip;</th><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td></tr></tbody></table></div>



Quantitative analysis - *the g(k) method*
===


Following the qualitative analysis presented above, a more quantitative approach will now be introduced. The goal of the so-called *g(k) method* is to quantitatively describe correlation between $\theta$ and $\phi$. This method is used for evaluation of proposed data-cuts. 

To derive the method, two figures are presented:
1. $f(\theta, \phi)$ - a 2D histogram with $\theta$ distribution on the x-axis and $\phi$ distribution on the y-axis
2. $g(k)$ - a 1D histogram with *so-called* $k$-line on the x-axis, where bin heights correspond to the integral over individual $k_i$ lines. $k$-lines are defined as diagonal lines in $f(\theta, \phi)$, fulfilling condition $k_i=\phi - \theta; (\phi - \theta) \in b_i$; binned as necessary, in general binning of $\Delta\phi = \Delta\theta = 1^{\circ}$ is used. 

## $f(\theta, \phi)$


```julia
nBins = 180

h2d1 = histogram2d(
    tree.thetaEmitted,
    tree.thetaEscaped;
    nbins          = (nBins, nBins),
    xlabel         = "θ",
    ylabel         = "ϕ",
    legend         = :topright,
    title          = string("\nf(ϕ, θ), ", nrow(tree), " entries"),
    lims           = (0, 180),
    aspect_ratio   = 1,
    c              = :coolwarm,
    colorbar_title = "\ncounts [#/°]",

)

plot!(y, c = :black, lw= 3, ls =:dash, label = "", legend = :topleft)
```




    
![svg](output_17_0.svg)
    



The correlation of $\theta$ and $\phi$ in $f(\theta, \phi)$ shows an S-shaped structure. There are two *hotspots* visible in the 2D histogram. First, a *smaller* hotspot is visible in the lower left quadrant of $f(\theta, \phi)$, here in-general $\phi > \theta$ - escape angles are overestimated over decay angles. The second *larger* hotspot is visible in the upper right quadrant of $f(\theta, \phi)$ where $\phi < \theta$ - escape angles are underestimed over decay angles. 

We can look more closely at the individual horizontal slices to view the $\theta$ distribution in each $\Delta\phi$. A few sample figures are shown. One can see that the $\theta$ distribution is very wide.


```julia
h1 = stephist(
    tree[ tree.thetaEscaped .< 5.0, :thetaEmitted ], 
    weights=tree[tree.thetaEscaped .< 5.0, :weights], label ="", c =:orange,  
    size=(1000,1400), legend = :topright, xlims = (0,180), ylims = (0, 200), 
    xlabel ="θ", ylabel ="", lw = 4, nbins = 180, dpi = 400 
)
vspan!([0,5], alpha = 0.2, c =:orange, label = "ϕ ∈ (0,5)" )

h2 = stephist( 
    tree[ 30 .< tree.thetaEscaped .< 35, :thetaEmitted ],
    weights=tree[30 .< tree.thetaEscaped .< 35, :weights],label ="", c =:red,
    size=(1000,1400), legend = :topright, xlims = (0,180), ylims = (0, 6000), 
    xlabel ="θ", ylabel ="", lw = 4, nbins = 180, dpi = 400 
)
vspan!([30,35], alpha = 0.2, c =:red, label = "ϕ ∈ (30,35)" )

h3 = stephist( 
    tree[ 85 .< tree.thetaEscaped .< 90, :thetaEmitted ],
    weights=tree[85 .< tree.thetaEscaped .< 90, :weights], label ="", c =:blue,
    size=(1000,1400), legend = :topright, xlims = (0,180), ylims = (0, 6000), 
    xlabel ="θ", ylabel ="", lw = 4, nbins = 180, dpi = 400 
)
vspan!([85,90], alpha = 0.2, c =:blue, label = "ϕ ∈ (85,90)" )

h4 = stephist( 
    tree[ 145 .< tree.thetaEscaped .< 150, :thetaEmitted ],
    weights=tree[145 .< tree.thetaEscaped .< 150, :weights],label ="", c =:green, legend = :topleft,
    size=(1000,1400), xlims = (0,180), ylims = (0, 9000), 
    xlabel ="θ", ylabel ="", lw = 4, nbins = 180, dpi = 400 
)
vspan!([145,150], alpha = 0.2, c =:green, label = "ϕ ∈ (145,150)" )


h5 = stephist( 
    tree[ 175 .< tree.thetaEscaped .< 180, :thetaEmitted ],
    weights=tree[175 .< tree.thetaEscaped .< 180, :weights],label ="", c =:black, legend = :topleft,
    size=(1000,1400), xlims = (0,180), ylims = (0, 1500), 
    xlabel ="θ", ylabel ="", lw = 4, nbins = 180, dpi = 400 
)
vspan!([175,180], alpha = 0.2, c =:black, label = "ϕ ∈ (175,180)" )

plot(h1, h2, h3, h4, h5, layout = @layout [a b ; c d ; e _])


```




    
![svg](output_20_0.svg)
    



Furthermore we can look even closer at each slice and calculate the statistical estimators. These will be used later in the analysis.


```julia
thetas = tree[20 .< tree.thetaEscaped .< 25, :thetaEmitted ]
# w = Weights(rho.(cosd.(thetas)))

stats =  get_slice_stats(
                            20.0,  # ϕmin
                            25.0,  # ϕmax
                            0,             # Emin
                            3500,          # Emax
                            thetas, 
                            dϕ,
                            )

h = Hist1D(thetas, (0:dϕ:180))
mediansFHist = StatsBase.median(h)

h1 = stephist( 
    thetas;
    nbins = 180, 
    xlims = (0,180), 
    ylims = (0, 1.2*maximum(bincounts(h))),
    label ="", 
    c =:orange, 
    legend= :topright, 
    lw = 4,
    xlabel = "θ",
    guidefontsize = 20,
    legendfontsize = 12,
    size = (800, 800)
    )

vspan!([20,25], alpha = 0.2, c =:orange, label = "ϕ ∈ (20,25)" )
vline!([stats[7]],lw =3, label ="median")
vline!([stats[6]],lw =3, label ="mode")
vline!([stats[5]],lw =3, label ="mean")

```




    
![svg](output_22_0.svg)
    



## k-lines

In the figure above of $f(\theta, \phi)$, a reference k = 0 line is shown by black dashed line. This line represents the bins where $\phi = \theta$, perfect correlation. In the figure below, two more k-lines are depicted. $k = -20$ and $k = +20$ lines are show in red and blue, respectively. These lines in turn represent $\phi - \theta = -20$ and $\phi - \theta = 20$. 


```julia
pts = 0:0.1:180

plot!(y, c = :black, lw= 3, ls =:dash, label = "k = 0", legend = :topleft)
plot!(pts, pts .+ 20,  lw= 4, c  =:red , label = "k = +20")
plot!(pts, pts .- 20,  lw= 4, c  =:blue, label = "k = -20")
```




    
![svg](output_25_0.svg)
    



## $g(k)$

As stated earlier, the $g(k)$ is a 1D histgoram of the integrals over the k-lines. The total number of k-lines to integrate $f(\theta, \phi)$ over is equal to $(180/\Delta\phi) * 2 + 1$, dependent on the binwidth $\Delta\phi$. To avoid double binning, $g(k)$ is not calculated from $f(\theta, \phi)$ itself, but rather from the definition of $k$-lines. Thus for each event, $\phi - \theta$ is calculated and **then** binned in the 1D histogram $g(k)$. The figure below shows $g(k)$ of the original $f(\theta, \phi)$ presented above. In the text following, some data-cuts will be introduced which produce different $f_i(\theta, \phi)$ distributions and corresponding $g_i(k)$ histograms.  

#### To quantitatively describe the correlation using $g(k)$-method, we calculate the RMS of the distribution as the standard error and the area representing the total amount of events that pass the cuts. 


```julia
difs0 = tree.thetaEscaped .- tree.thetaEmitted    # array of ϕ - θ
rms0  = round(get_rms(difs0), digits = 2)
h0    = StatsBase.fit(Histogram, difs0, -180:180)
nEvents = nrow(tree)

histogram(difs0, 
    xlims  = (-180,180), 
    ylims  = (0, 1.1*maximum(h0.weights)), 
    nbins  = length(h0.edges[1]),
    xlabel = "k [°]",
    ylabel = "g0(k) [#/°]",
    label  = "g0(k)",
    lw     = 4,
    fillrange = 0,
    fillalpha = 0.2
)

annotate!((-150, 1.7e5, text("RMS₀ = $rms0 °", 14, :left)))
annotate!((-150, 1.6e5, text("nEvents = $nEvents ", 14, :left)))
```




    
![svg](output_28_0.svg)
    



The figure above depicts the $g(k)$ histogram calculated from $f(\theta, \phi)$. The $g(k)$ distribution is centered around $k = 0$ line, which represents perfect correlation. However, the distribution is quite wide as the $RMS = 46.99^{\circ}$. All events from the original distribution are represented. 

$g(k)$ - analysis of energy cuts
===

The goal of calculating $g(k)$ and its $RMS$ is to evaluate different data-cuts. In the analysis below, various energy cuts are tested. Six data-cuts are presented - the sum of electron energies $E_{sum}$ must fulfill $E_{sum} \in (500*i, 500*i + \Delta E); \Delta E = 500 keV; i = (1, 2, 3, 4, 5, 6) keV$. (Events with $E_{sum} \in (0, 500) keV$ are omitted as there is not enough statistics). 

First, a set of $f\_i(\theta, \phi)$ is presented.  


```julia
hs = []
for e in 500:500:3000
    e2 = e+500
    gdf = @chain tree begin
        @subset((e .<= :ESum .<= e2))
        @select(:thetaEscaped, :thetaEmitted, :ESum)
    end

    push!(hs, histogram2d(gdf.thetaEmitted, gdf.thetaEscaped, 
                 nbins = nBins, lims=(0, 180), xlabel="θ [°]", ylabel="ϕ [°]",
                 c = :coolwarm, title ="Esum ∈ ($e, $e2) keV"
        )
    )
end
```


```julia
plot(hs[1], hs[2], hs[3], hs[4], hs[5], hs[6], size = (1000, 1200), layout = @layout [a b; c d; e f ])
```




    
![svg](output_32_0.svg)
    



Now the corresponding $g(k)$


```julia
p = plot(size = (800, 800), legend=:topright, xlims=(-180, 180), xlabel="k [°]", ylabel="g(k) [#/°]", lw = 4)
for (i,e) in enumerate(500:500:3000)
    e2 = e+500
    gdf = @chain tree begin
        @subset((e .<= :ESum .<= e2))
        @select(:thetaEscaped, :thetaEmitted, :ESum)
    end

    difs   = gdf.thetaEscaped .- gdf.thetaEmitted
    rms    = round(get_rms(difs), digits = 2)
    h1     = StatsBase.fit(Histogram, difs, -180:180)
    
    stephist!(p, difs, 
        nbins  = Int(180/dϕ*2+1), 
        lw     = 4, 
        label  ="g$i; E ∈ ($e, $e2), \n rms = $rms",
        legend =:topright, 
        xlims  = (-180,180), 
        xlabel = "k [°]",
        ylims  = (0, 8e4)
    )

end
p

```




    
![svg](output_34_0.svg)
    



It is visible from both figures $f_i(\theta, \phi)$ that the number of events increases with increasing energy. 

The $RMS$ is represented in the legend. Again, with increasing energy data-cut $RMS$ decreases. To view the width of each distribution, the histograms are normalized to area of 1.


```julia
p = plot(size = (800, 800), legend=:topright, xlims=(-180, 180), xlabel="k [°]", ylabel="g(k) [#/°]", lw = 4)
for (i,e) in enumerate(500:500:3000)
    e2 = e+500
    gdf = @chain tree begin
        @subset((e .<= :ESum .<= e2))
        @select(:thetaEscaped, :thetaEmitted, :ESum)
    end

    difs   = gdf.thetaEscaped .- gdf.thetaEmitted
    rms    = round(get_rms(difs), digits = 2)
    h1     = StatsBase.fit(Histogram, difs, -180:180)
    
    stephist!(p, difs, 
              nbins  = Int(180/dϕ*2+1), 
              lw     = 4, 
              label  = "g$i; E ∈ ($e, $e2), \n rms = $rms",
              legend = :topright, 
              xlims  = (-180,180), 
              xlabel = "k [°]",
              norm   = :pdf,
              ylims  = (0, 0.016)
    )

end
p
```




    
![svg](output_36_0.svg)
    



### From the figure, it is visible that increasing the energy cut improves the corellation. The $g(k)$-method can be used for quantification of correlation.

Slicing horizontally 
===

We can see that while applying an energy cut on the data results in decreased statistics, it did provide for a better reconstruction precision. We thus have a tool for comparing the effects of data cuts on the data.

Next we look more in detail at individual $\Delta\phi$ slices. We will slice up $f(\theta,\phi)$ horizontally in slices of $\Delta\phi$ = $1^{\circ}$. However, for better visualisation of what is happening we first show a horizontal slice with $\phi \in (10, 15)\deg$ and its $g(k)$ as:


```julia
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
```




    
![svg](output_41_0.svg)
    



This procedure is repeated for each slice with $\Delta\phi = 1^{\circ}$. Slicing $f(\phi, \theta)$ horizontally to cover the whole 0 - 180 degree range yields $g(k)$s:


```julia
p = plot(size = (800, 800), legend=:none, xlims=(-180, 180), xlabel="k [°]", ylabel="g(k) [#/°]", lw = 4)
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
        ylims  = (0, 2e3)
    )

end
p
```




    
![svg](output_43_0.svg)
    



Not much can be deduced in this example. So many lines are difficult to decipher. However, if one were to look at the individual $\Delta\phi$ cuts as a new dimension, we can look at the graph in the plane of $(k, \Delta\phi)$ with z-direction being the value of $g(k, \Delta\phi)$. 


```julia
dfGOrig = get_gs_df(tree, dϕ, sign)
matGOrig = df_to_mat(dfGOrig);
```


```julia
xRange = dϕ-180:dϕ:180-dϕ
yRange = 0:dϕ:180-dϕ
sf1 = surface(xRange, yRange,  matGOrig, legend =:none, xlabel ="k", ylabel ="Δϕ", zlabel="g(k, Δϕ)", tickfontsize = 12, c= :turbo)
hm1 = heatmap(xRange, yRange,  matGOrig, ylabel ="Δϕ", xlabel ="k" , tickfontsize = 10, bottom_margin  = 6Plots.mm, c= :turbo)
vline!([0], label ="", c = :black, lw = 2, s=:dash)
plot(sf1,hm1, size =(1000,400), layout = @layout [a{0.4w} b])
```




    
![svg](output_46_0.svg)
    



Now we can see a few important features. First of all, there are two peaks visible in the left figure, with the higher peak (more statistics) being in the region of $130^{\circ} < \phi < 17^{\circ}0$ . Secondly,  we can see the deviation of the peaks from the `` k = 0`` line in the right figure. There are two hotspots visible. First hotspot (corresponding to the lower peak in  figure) is centered around $\Delta\phi \approx 30^{\circ}$ and is shifted slightly to the right of the ``k = 0`` line. The escaped angle overestimates the emitted angle. Second hotspot (corresponding to the higher peak in figure) is centered around $\Delta\phi \approx 150^{\circ}$ and is shifted visibly to the left of the ``k = 0`` line. The escaped angle underestimates the emitted angle. Lastly, we can see that the regions $\phi \approx 0^{\circ}$ and $\phi \approx 180^{\circ}$ are squeezed toward higher, lower angles, respectively. 

Furthermore, we can also look at how the individual $\Delta\phi$ slices look in terms of statistical variables (mean, mode, median). For each variable, obtained from the $\phi(\theta)$ distributions. 


```julia
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
```


```julia
scatter( means, xPts, ms = 4, label = "means" , xlabel = "mean g_i(k)", ylabel ="Δϕ_i(k)", aspect_ratio =1, legend = :topleft, lims = (0,180))
scatter!( modes, xPts, ms = 3, label ="modes" )
scatter!( medians, xPts, ms = 3, label ="medians" )
```




    
![svg](output_50_0.svg)
    



Shift to reduce RMS
===


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


```julia
res_means = [ means[i] - y(i*dϕ) for i in 1:length(means) ]
res_modes = [ modes[i] - y(i*dϕ) for i in 1:length(modes) ]
res_medians = [ medians[i] - y(i*dϕ) for i in 1:length(medians) ]

scatter( res_means, xPts, ms = 3, label = "mean" , xlabel = "shift", ylabel ="Δϕ_i(k)", ylims = (0,180), aspect_ratio =1, legend = :topright)
scatter!( res_modes, xPts, ms = 3, label ="mode" )
scatter!( res_medians, xPts, ms = 3, label ="median" )
vline!([0], label ="", c = :black, lw = 3, s=:dash)
```




    
![svg](output_53_0.svg)
    



We can see that for each $\Delta\phi_i$ slice the three estimators provide different values to shift the angles by. The most drastic shift (ie. farthest away from ``k=0``) is given by mean, the least on the other hand by mode. 

To avoid undesirable discretization of our data, we fit the shifts. We also flipped the axes so that we get $s(\Delta\phi_i)$.


```julia
scatter(  xPts,res_means, ms = 3, label = "shift by" , ylabel = "shift", xlabel ="Δϕ_i(k)", xlims = (0,180), aspect_ratio =1, legend = :topright)
f1 = Polynomials.fit(xPts,  res_means, 6 )

plot!(f1, extrema(xPts)..., label="Fit")
```




    
![svg](output_56_0.svg)
    




Now we shift each $\phi$ in the original data set to obtain a new set of $\phi'$, we do so by $\phi' = \phi +s$. 

*Shift by mean*
====


```julia
modTree2 = @chain tree begin
    @select(:thetaEmitted, :thetaEscaped, :weights)
    @rtransform :bin =  get_bin_center(:thetaEscaped, Int(180/dϕ))  # create a vector of bin centers (which bin ϕ falls inside)
    @transform :thetaEscapedOld = :thetaEscaped                     # create a copy of the old ϕ (for comparison only)

    @rtransform :thetaEscapedDisc = :thetaEscapedOld + res_means[Int(ceil(:bin/dϕ))] # shift ϕ by s: ϕ' = ϕ + s 
    @rtransform :thetaEscaped = :thetaEscapedOld + f1(:thetaEscapedOld) # shift ϕ by s: ϕ' = ϕ + s 
    @subset( 0 .< :thetaEscaped .< 180) # keep only physical angles
end
```




<div class="data-frame"><p>15,747,968 rows × 6 columns</p><table class="data-frame"><thead><tr><th></th><th>thetaEmitted</th><th>thetaEscaped</th><th>weights</th><th>bin</th><th>thetaEscapedOld</th><th>thetaEscapedDisc</th></tr><tr><th></th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Int64">Int64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th></tr></thead><tbody><tr><th>1</th><td>89.6031</td><td>92.2988</td><td>1</td><td>94.5</td><td>94.6325</td><td>92.5524</td></tr><tr><th>2</th><td>96.965</td><td>120.357</td><td>1</td><td>139.5</td><td>139.497</td><td>120.722</td></tr><tr><th>3</th><td>154.233</td><td>119.633</td><td>1</td><td>138.5</td><td>138.388</td><td>120.144</td></tr><tr><th>4</th><td>67.3403</td><td>129.17</td><td>1</td><td>154.5</td><td>154.281</td><td>129.132</td></tr><tr><th>5</th><td>11.2133</td><td>67.4548</td><td>1</td><td>53.5</td><td>53.0724</td><td>67.4527</td></tr><tr><th>6</th><td>168.417</td><td>84.8354</td><td>1</td><td>81.5</td><td>81.0901</td><td>85.1413</td></tr><tr><th>7</th><td>122.632</td><td>78.2451</td><td>1</td><td>69.5</td><td>69.6012</td><td>79.0214</td></tr><tr><th>8</th><td>4.96695</td><td>53.6993</td><td>1</td><td>34.5</td><td>34.9839</td><td>53.8336</td></tr><tr><th>9</th><td>79.857</td><td>72.9471</td><td>1</td><td>61.5</td><td>61.1214</td><td>73.0331</td></tr><tr><th>10</th><td>43.0942</td><td>52.9266</td><td>1</td><td>33.5</td><td>33.994</td><td>53.2077</td></tr><tr><th>11</th><td>35.1387</td><td>97.6664</td><td>1</td><td>104.5</td><td>104.069</td><td>97.3047</td></tr><tr><th>12</th><td>150.65</td><td>96.4213</td><td>1</td><td>101.5</td><td>101.923</td><td>96.4782</td></tr><tr><th>13</th><td>56.3877</td><td>69.2744</td><td>1</td><td>55.5</td><td>55.6631</td><td>69.4157</td></tr><tr><th>14</th><td>128.995</td><td>92.4145</td><td>1</td><td>94.5</td><td>94.8403</td><td>92.7602</td></tr><tr><th>15</th><td>63.2692</td><td>58.3689</td><td>1</td><td>40.5</td><td>40.9504</td><td>58.4208</td></tr><tr><th>16</th><td>40.4733</td><td>68.9162</td><td>1</td><td>55.5</td><td>55.1476</td><td>68.9002</td></tr><tr><th>17</th><td>124.468</td><td>122.461</td><td>1</td><td>142.5</td><td>142.771</td><td>123.096</td></tr><tr><th>18</th><td>128.263</td><td>43.44</td><td>1</td><td>20.5</td><td>20.8424</td><td>43.8765</td></tr><tr><th>19</th><td>105.706</td><td>135.303</td><td>1</td><td>172.5</td><td>172.673</td><td>135.385</td></tr><tr><th>20</th><td>142.67</td><td>102.768</td><td>1</td><td>112.5</td><td>112.558</td><td>102.595</td></tr><tr><th>&vellip;</th><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td></tr></tbody></table></div>




```julia
dfG2 = get_gs_df(modTree2, dϕ, sign)
matG2 = df_to_mat(dfG2);
```

We look at the $f(\theta, \phi')$ figure.


```julia
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
```




    
![svg](output_62_0.svg)
    



We can see that shifting by mean value resulted in *squeezing* the phase-space. We have reduced the range of angles which we can interpret in our measuremt. However, this should lead toward reduces RMS. We look at that in the following figures.

First, for comparison the original dataset with $f(\theta, \phi)$, $g(k, \Delta\phi)$, calculated $RMS$ for each $g_i(k)$ and total $RMS$.  


```julia
h2d1 = histogram2d(
    tree.thetaEmitted,
    tree.thetaEscaped;
    nbins        = (nBins, nBins),
    weights      = tree.weights,
    xlabel       = "θ",
    ylabel       = "ϕ",
    legend       = :topright,
    title        = string("f(θ, ϕ), ", nrow(tree), " entries"),
    lims         = (0, 180),
    aspect_ratio = 1,
    c = :turbo
)
plot!(xPts, xPts, label ="", c= :black, style= :dash, lw =3)
rms1 = [ get_rms(dfGOrig[:,i], dfGOrig[:,1]) for i in 2:ncol(dfGOrig) ]

sct1 = scatter( xPts, rms1, label ="rms(Δϕ)", ms=2, legend=:top, xlabel ="Δϕ slice", ylabel ="rms", c= :red )

fh2d1 = Hist2D((tree.thetaEmitted,tree.thetaEscaped), 
    Weights(tree.weights),
    (minAngle:dEmitted:maxAngle, minAngle:dEmitted:maxAngle))
gs1 = get_diagonal_sums(fh2d1)
ks1 = get_k_factors(fh2d1);

gk1 = plot(ks1 .* dEmitted, gs1, legend=:topright, xlims=(-179, 179), xlabel="k-factor", ylabel="g(k)", label="g_1(k)")

# plot(title, h2d1, hm1, sct1, layout = @layout[a{0.05h};b c; d _] , size = (1100, 800))
plot(h2d1, hm1, sct1, gk1, layout = @layout[a b; c d] , size = (1200, 800), plot_title= "Unmodified angles")
```




    
![svg](output_65_0.svg)
    



And the **modified** dataset.


```julia
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
```




    
![svg](output_67_0.svg)
    




```julia
rmsTotalUnModded = round(get_rms(gs1, ks1 .* dEmitted ), digits = 2)
```




    46.99




```julia
rmsTotalModded = round(get_rms(gs2, ks2 .* dEmitted ), digits = 2)
```




    43.85



Now we compare the four $RMS$ figures together. 


```julia
scComp2 = scatter( [xPts xPts], [rms1 rms2],  label =["rms_1 = $rmsTotalUnModded °" "rms_2 = $rmsTotalModded °"], 
        ms=3, c = [:red :blue], bottom_margin = 8Plots.mm,
        legend=:topleft, xlabel ="Δϕ_i slice", ylabel ="rms", ylims = (yMin, 1.1*maximum(rms1)))

gkComp2 = plot( [ks1 .* dEmitted, ks2 .* dEmitted], [gs1, gs2], label = ["g_1(k)" "g_2(k)"], legend = :topleft,
                xlabel = "k", ylabel = "g(k)", seriestype = :stepmid, lw = 4, c = [:red :blue])

plot(scComp2, gkComp2, size = (1200, 600))
```




    
![svg](output_71_0.svg)
    



We can see that shifting by mean value results in reduced $RMS$, the goal is achieved. 

## Finally, we can now provide a function, which as an input takes the measured angle $\phi$ and as an output provides $\phi'$ (the most likely $\theta$): **$\phi'(\phi)$**.


```julia
plot(
    xPts, xPts .+ f1.(xPts), lims = (0,180), lw = 4, c= :blue, 
    label ="", xlabel = "ϕ", ylabel ="ϕ'", aspect_ratio = 1, 
    title = "ϕ'(ϕ); measured (ϕ) vs reported (ϕ') angle"
)

```




    
![svg](output_74_0.svg)
    




```julia

```
