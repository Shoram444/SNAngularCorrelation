



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



The 2d Histogram of $\phi$ vs $\theta$ is defined to be $f(\theta,\phi)$. For each combination of $\phi$ and $\theta$, the bin number is obtained as corresponding value of $f(\theta,\phi)$.

<div class="data-frame"><p>15,747,968 rows × 18 columns</p><table class="data-frame"><thead><tr><th></th><th>momentumEmitted2x</th><th>reconstructedEnergy2</th><th>momentumEscaped1z</th><th>momentumEmitted1x</th><th>thetaEmitted</th><th>momentumEscaped1y</th><th>momentumEmitted1z</th><th>momentumEmitted1y</th><th>momentumEscaped2z</th><th>thetaEscaped</th><th>momentumEscaped2y</th><th>momentumEmitted2y</th><th>reconstructedEnergy1</th><th>momentumEscaped1x</th><th>momentumEmitted2z</th><th>momentumEscaped2x</th><th>ESum</th><th>weights</th></tr><tr><th></th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Float64">Float64</th><th title="Int64">Int64</th></tr></thead><tbody><tr><th>1</th><td>0.723365</td><td>604.269</td><td>-0.471071</td><td>-1.41249</td><td>89.6031</td><td>2.26532</td><td>-0.946003</td><td>2.12515</td><td>-0.155021</td><td>94.6325</td><td>0.351932</td><td>0.615735</td><td>2258.02</td><td>-1.36057</td><td>0.283402</td><td>0.777892</td><td>2862.28</td><td>1</td></tr><tr><th>2</th><td>-1.87869</td><td>1651.11</td><td>0.111843</td><td>0.620354</td><td>96.965</td><td>-0.427881</td><td>0.585939</td><td>-0.857937</td><td>0.547949</td><td>139.497</td><td>-0.460427</td><td>-0.935666</td><td>802.522</td><td>1.07408</td><td>0.0929118</td><td>-1.94376</td><td>2453.63</td><td>1</td></tr><tr><th>3</th><td>-0.451873</td><td>450.772</td><td>-0.204413</td><td>0.143294</td><td>154.233</td><td>-0.344668</td><td>0.461237</td><td>-0.617821</td><td>-0.00393392</td><td>138.388</td><td>-0.0885292</td><td>0.444109</td><td>425.001</td><td>0.557768</td><td>-0.51231</td><td>-0.766414</td><td>875.773</td><td>1</td></tr><tr><th>4</th><td>0.00622859</td><td>2520.27</td><td>0.0712075</td><td>-0.27381</td><td>67.3403</td><td>0.440423</td><td>-0.0579716</td><td>-1.01549</td><td>0.518941</td><td>154.281</td><td>-2.14362</td><td>-1.34793</td><td>659.761</td><td>0.717866</td><td>2.66656</td><td>-1.67234</td><td>3180.03</td><td>1</td></tr><tr><th>5</th><td>-0.660664</td><td>834.025</td><td>-0.282475</td><td>-1.77769</td><td>11.2133</td><td>-2.1054</td><td>-0.0196651</td><td>-1.94104</td><td>0.221755</td><td>53.0724</td><td>-0.966977</td><td>-1.05076</td><td>2170.29</td><td>-1.37587</td><td>0.0859927</td><td>0.291419</td><td>3004.32</td><td>1</td></tr><tr><th>6</th><td>-0.44519</td><td>1530.67</td><td>0.917931</td><td>0.0695716</td><td>168.417</td><td>1.51856</td><td>2.09515</td><td>0.648145</td><td>-1.07627</td><td>81.0901</td><td>0.0381113</td><td>-0.472476</td><td>1741.93</td><td>-0.999117</td><td>-1.86705</td><td>-1.51857</td><td>3272.61</td><td>1</td></tr><tr><th>7</th><td>-0.445143</td><td>1784.88</td><td>0.277989</td><td>0.889556</td><td>122.632</td><td>-0.0529092</td><td>0.257764</td><td>-0.238635</td><td>0.00632242</td><td>69.6012</td><td>1.73728</td><td>1.75164</td><td>573.353</td><td>0.858823</td><td>-1.32043</td><td>0.80815</td><td>2358.23</td><td>1</td></tr><tr><th>8</th><td>-1.42443</td><td>1723.99</td><td>0.250602</td><td>-1.16971</td><td>4.96695</td><td>0.834453</td><td>0.0445409</td><td>1.26977</td><td>-0.894417</td><td>34.9839</td><td>1.23976</td><td>1.64022</td><td>1290.01</td><td>-1.4441</td><td>-0.121338</td><td>-1.4647</td><td>3014.0</td><td>1</td></tr><tr><th>9</th><td>-0.0553668</td><td>1281.87</td><td>0.457432</td><td>-0.344509</td><td>79.857</td><td>1.30726</td><td>-0.133961</td><td>1.66448</td><td>0.94093</td><td>61.1214</td><td>0.715475</td><td>0.160922</td><td>1268.96</td><td>0.500476</td><td>-1.71006</td><td>-0.744009</td><td>2550.83</td><td>1</td></tr><tr><th>10</th><td>-1.7249</td><td>2124.13</td><td>-0.525607</td><td>-0.976017</td><td>43.0942</td><td>0.765786</td><td>-0.71208</td><td>-0.308521</td><td>-2.15062</td><td>33.994</td><td>0.664738</td><td>1.23561</td><td>836.582</td><td>-0.544891</td><td>-1.47674</td><td>-1.08732</td><td>2960.71</td><td>1</td></tr><tr><th>11</th><td>-0.19787</td><td>467.234</td><td>-0.513265</td><td>-0.80288</td><td>35.1387</td><td>-0.210073</td><td>-0.13689</td><td>-0.695541</td><td>-0.230839</td><td>104.069</td><td>-0.438105</td><td>-0.806524</td><td>675.699</td><td>-0.763111</td><td>-0.0786541</td><td>0.483935</td><td>1142.93</td><td>1</td></tr><tr><th>12</th><td>-0.529524</td><td>2307.59</td><td>-0.682419</td><td>-0.307922</td><td>150.65</td><td>-0.818542</td><td>-0.238325</td><td>-1.27159</td><td>1.96584</td><td>101.923</td><td>0.305447</td><td>2.70679</td><td>913.668</td><td>0.601471</td><td>-0.27609</td><td>1.57859</td><td>3221.26</td><td>1</td></tr><tr><th>13</th><td>1.22009</td><td>976.306</td><td>-0.64402</td><td>0.895972</td><td>56.3877</td><td>0.122387</td><td>-0.727497</td><td>0.989891</td><td>0.620286</td><td>55.6631</td><td>0.386989</td><td>0.454814</td><td>1093.07</td><td>1.29631</td><td>0.505461</td><td>1.11372</td><td>2069.37</td><td>1</td></tr><tr><th>14</th><td>1.62335</td><td>1841.63</td><td>-0.772472</td><td>0.148038</td><td>128.995</td><td>-1.01428</td><td>-0.825084</td><td>-1.266</td><td>0.520843</td><td>94.8403</td><td>1.28805</td><td>1.23633</td><td>1091.04</td><td>0.781192</td><td>1.05356</td><td>1.81786</td><td>2932.68</td><td>1</td></tr><tr><th>15</th><td>-0.975718</td><td>619.331</td><td>-1.0753</td><td>-0.587855</td><td>63.2692</td><td>1.86365</td><td>-1.97465</td><td>1.62362</td><td>-0.608049</td><td>40.9504</td><td>0.0886557</td><td>0.0910889</td><td>2161.47</td><td>-1.43294</td><td>-0.23707</td><td>-0.75854</td><td>2780.8</td><td>1</td></tr><tr><th>16</th><td>-0.181153</td><td>985.07</td><td>1.34426</td><td>-1.8497</td><td>40.4733</td><td>-0.640471</td><td>1.69434</td><td>-0.0995104</td><td>0.984389</td><td>55.1476</td><td>0.757152</td><td>0.0829501</td><td>2050.88</td><td>-1.97838</td><td>1.39191</td><td>-0.54679</td><td>3035.95</td><td>1</td></tr><tr><th>17</th><td>1.62458</td><td>2530.76</td><td>0.0707918</td><td>-1.20434</td><td>124.468</td><td>0.512761</td><td>0.35071</td><td>0.34492</td><td>1.19175</td><td>142.771</td><td>-2.09416</td><td>-2.11764</td><td>886.684</td><td>-1.14899</td><td>1.36657</td><td>1.72853</td><td>3417.44</td><td>1</td></tr><tr><th>18</th><td>0.685257</td><td>639.204</td><td>0.0196268</td><td>-0.00221343</td><td>128.263</td><td>0.374426</td><td>0.395698</td><td>1.14861</td><td>-0.17065</td><td>20.8424</td><td>0.0968995</td><td>-0.742877</td><td>806.954</td><td>0.864823</td><td>0.20101</td><td>0.966251</td><td>1446.16</td><td>1</td></tr><tr><th>19</th><td>-0.211174</td><td>628.139</td><td>0.502857</td><td>-2.74734</td><td>105.706</td><td>2.36161</td><td>-0.622121</td><td>1.39787</td><td>-0.258807</td><td>172.673</td><td>-0.709862</td><td>-0.736769</td><td>2674.92</td><td>-1.86425</td><td>0.670143</td><td>0.514478</td><td>3303.06</td><td>1</td></tr><tr><th>20</th><td>0.0476415</td><td>707.158</td><td>-0.736812</td><td>-1.20521</td><td>142.67</td><td>1.44929</td><td>-0.496897</td><td>1.82196</td><td>0.813798</td><td>112.558</td><td>0.299922</td><td>-0.860809</td><td>1786.85</td><td>-1.10866</td><td>0.692477</td><td>0.550832</td><td>2494.01</td><td>1</td></tr><tr><th>&vellip;</th><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td><td>&vellip;</td></tr></tbody></table></div>



Quantitative analysis - *the g(k) method*
===


Following the qualitative analysis presented above, a more quantitative approach will now be introduced. The goal of the so-called *g(k) method* is to quantitatively describe correlation between $\theta$ and $\phi$. This method is used for evaluation of proposed data-cuts. 

To derive the method, two figures are presented:
1. $f(\theta, \phi)$ - a 2D histogram with $\theta$ distribution on the x-axis and $\phi$ distribution on the y-axis
2. $g(k)$ - a 1D histogram with *so-called* $k$-line on the x-axis, where bin heights correspond to the integral over individual $k_i$ lines. $k$-lines are defined as diagonal lines in $f(\theta, \phi)$, fulfilling condition $k_i=\phi - \theta; (\phi - \theta) \in b_i$; binned as necessary, in general binning of $\Delta\phi = \Delta\theta = 1^{\circ}$ is used. 

## $f(\theta, \phi)$

    
![svg](output_17_0.svg)
    



The correlation of $\theta$ and $\phi$ in $f(\theta, \phi)$ shows an S-shaped structure. There are two *hotspots* visible in the 2D histogram. First, a *smaller* hotspot is visible in the lower left quadrant of $f(\theta, \phi)$, here in-general $\phi > \theta$ - escape angles are overestimated over decay angles. The second *larger* hotspot is visible in the upper right quadrant of $f(\theta, \phi)$ where $\phi < \theta$ - escape angles are underestimed over decay angles. 

We can look more closely at the individual horizontal slices to view the $\theta$ distribution in each $\Delta\phi$. A few sample figures are shown. One can see that the $\theta$ distribution is very wide.

    
![svg](output_20_0.svg)
    



Furthermore we can look even closer at each slice and calculate the statistical estimators. These will be used later in the analysis.



    
![svg](output_22_0.svg)
    



## k-lines

In the figure above of $f(\theta, \phi)$, a reference k = 0 line is shown by black dashed line. This line represents the bins where $\phi = \theta$, perfect correlation. In the figure below, two more k-lines are depicted. $k = -20$ and $k = +20$ lines are show in red and blue, respectively. These lines in turn represent $\phi - \theta = -20$ and $\phi - \theta = 20$. 




    
![svg](output_25_0.svg)
    



## $g(k)$

As stated earlier, the $g(k)$ is a 1D histgoram of the integrals over the k-lines. The total number of k-lines to integrate $f(\theta, \phi)$ over is equal to $(180/\Delta\phi) * 2 + 1$, dependent on the binwidth $\Delta\phi$. To avoid double binning, $g(k)$ is not calculated from $f(\theta, \phi)$ itself, but rather from the definition of $k$-lines. Thus for each event, $\phi - \theta$ is calculated and **then** binned in the 1D histogram $g(k)$. The figure below shows $g(k)$ of the original $f(\theta, \phi)$ presented above. In the text following, some data-cuts will be introduced which produce different $f_i(\theta, \phi)$ distributions and corresponding $g_i(k)$ histograms.  

#### To quantitatively describe the correlation using $g(k)$-method, we calculate the RMS of the distribution as the standard error and the area representing the total amount of events that pass the cuts. 

    
![svg](output_28_0.svg)
    



The figure above depicts the $g(k)$ histogram calculated from $f(\theta, \phi)$. The $g(k)$ distribution is centered around $k = 0$ line, which represents perfect correlation. However, the distribution is quite wide as the $RMS = 46.99^{\circ}$. All events from the original distribution are represented. 

$g(k)$ - analysis of energy cuts
===

The goal of calculating $g(k)$ and its $RMS$ is to evaluate different data-cuts. In the analysis below, various energy cuts are tested. Six data-cuts are presented - the sum of electron energies $E_{sum}$ must fulfill $E_{sum} \in (500*i, 500*i + \Delta E); \Delta E = 500 keV; i = (1, 2, 3, 4, 5, 6) keV$. (Events with $E_{sum} \in (0, 500) keV$ are omitted as there is not enough statistics). 

First, a set of $f\_i(\theta, \phi)$ is presented.  


    
![svg](output_32_0.svg)
    



Now the corresponding $g(k)$


    
![svg](output_34_0.svg)
    



It is visible from both figures $f_i(\theta, \phi)$ that the number of events increases with increasing energy. 

The $RMS$ is represented in the legend. Again, with increasing energy data-cut $RMS$ decreases. To view the width of each distribution, the histograms are normalized to area of 1.



    
![svg](output_36_0.svg)
    



### From the figure, it is visible that increasing the energy cut improves the corellation. The $g(k)$-method can be used for quantification of correlation.

Slicing horizontally 
===

We can see that while applying an energy cut on the data results in decreased statistics, it did provide for a better reconstruction precision. We thus have a tool for comparing the effects of data cuts on the data.

Next we look more in detail at individual $\Delta\phi$ slices. We will slice up $f(\theta,\phi)$ horizontally in slices of $\Delta\phi$ = $1^{\circ}$. However, for better visualisation of what is happening we first show a horizontal slice with $\phi \in (10, 15)\deg$ and its $g(k)$ as:



    
![svg](output_41_0.svg)
    



This procedure is repeated for each slice with $\Delta\phi = 1^{\circ}$. Slicing $f(\phi, \theta)$ horizontally to cover the whole 0 - 180 degree range yields $g(k)$s:



    
![svg](output_43_0.svg)
    



Not much can be deduced in this example. So many lines are difficult to decipher. However, if one were to look at the individual $\Delta\phi$ cuts as a new dimension, we can look at the graph in the plane of $(k, \Delta\phi)$ with z-direction being the value of $g(k, \Delta\phi)$. 


    
![svg](output_46_0.svg)
    



Now we can see a few important features. First of all, there are two peaks visible in the left figure, with the higher peak (more statistics) being in the region of $130^{\circ} < \phi < 17^{\circ}0$ . Secondly,  we can see the deviation of the peaks from the `` k = 0`` line in the right figure. There are two hotspots visible. First hotspot (corresponding to the lower peak in  figure) is centered around $\Delta\phi \approx 30^{\circ}$ and is shifted slightly to the right of the ``k = 0`` line. The escaped angle overestimates the emitted angle. Second hotspot (corresponding to the higher peak in figure) is centered around $\Delta\phi \approx 150^{\circ}$ and is shifted visibly to the left of the ``k = 0`` line. The escaped angle underestimates the emitted angle. Lastly, we can see that the regions $\phi \approx 0^{\circ}$ and $\phi \approx 180^{\circ}$ are squeezed toward higher, lower angles, respectively. 

Furthermore, we can also look at how the individual $\Delta\phi$ slices look in terms of statistical variables (mean, mode, median). For each variable, obtained from the $\phi(\theta)$ distributions. 




    
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






    
![svg](output_53_0.svg)
    



We can see that for each $\Delta\phi_i$ slice the three estimators provide different values to shift the angles by. The most drastic shift (ie. farthest away from ``k=0``) is given by mean, the least on the other hand by mode. 

To avoid undesirable discretization of our data, we fit the shifts. We also flipped the axes so that we get $s(\Delta\phi_i)$.


    
![svg](output_56_0.svg)
    




Now we shift each $\phi$ in the original data set to obtain a new set of $\phi'$, we do so by $\phi' = \phi +s$. 

*Shift by mean*
====



We look at the $f(\theta, \phi')$ figure.




    
![svg](output_62_0.svg)
    



We can see that shifting by mean value resulted in *squeezing* the phase-space. We have reduced the range of angles which we can interpret in our measuremt. However, this should lead toward reduces RMS. We look at that in the following figures.

First, for comparison the original dataset with $f(\theta, \phi)$, $g(k, \Delta\phi)$, calculated $RMS$ for each $g_i(k)$ and total $RMS$.  




    
![svg](output_65_0.svg)
    



And the **modified** dataset.





    
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




    
![svg](output_71_0.svg)
    



We can see that shifting by mean value results in reduced $RMS$, the goal is achieved. 

## Finally, we can now provide a function, which as an input takes the measured angle $\phi$ and as an output provides $\phi'$ (the most likely $\theta$): **$\phi'(\phi)$**.




    
![svg](output_74_0.svg)
    

