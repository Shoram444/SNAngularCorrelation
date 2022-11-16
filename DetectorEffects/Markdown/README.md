Detector Effects 
===============	

Study of Detector Effects on the angular correlation between decay ($\theta$) and escape ($\phi$) angles for the SuperNEMO experiment. $\theta$ is defined as the angle between the two emitted electrons as they are emitted in the $\beta\beta$ decay - the initial angle. $\phi$ is defined as the angle between the two electrons as they escape the source foil and enter the tracker - the final angle. The angles are obtained from the scalar products of the momentum vectors as found in the relevant G4Step in the simulation process. 


$cos(\theta) = \frac{p_1 \cdot p_2}{|p_1||p_2|}$


In order to get study only the detector effects separated from the physics of $\beta\beta$ decay, we have simulated $10^8$ events with uniform angular distribution as well as uniform single electron spectra. This is to make sure that each angle is sampled with the same probability and that the input events are uniformly distributed, thus we do not favour certain events as would be the case for standard $\beta\beta$ decay which has a preference for back-to-back electrons. 



Once the events have been simulated in Falaise, a number of data-cuts has been applied in order to filter only the $\beta\beta$ candidate events. A possible $\beta\beta$ candidate must have:

 1. Two negatively charged tracks reconstructed,
 2. Two vertices on the source foil, within given distance from each other,
 3. Sum of electron energies within the range: $0 ~keV \leq E_{sum} \leq 3500 ~keV$,
 4. Two individual Optical Module hits,
 5. Two associated Optical Module hits. 

 The result of applying the data-cuts is a change in $\theta$ distribution from uniform to "*increasing toward higher angles*". The comparison is shown in the following figure:


    
![svg](output_8_0.svg)
    



Two important features can be noticed in the figure. First, the overall number of events has been reduced from $10^8$ to $\sim 1.6 \cdot 10^7$ events. The efficiency of the applied data-cuts is around $16 \%$. Second, the spectrum shape has been changed - from uniform INPUT spectrum to spectrum in which higher angles are prefered.

 Once the path of the particles through the detector has been simulated using Falaise, $\phi$ can be extracted and studied. The comparison of $\theta$ and $\phi$ distributions is shown below. From now on, the OUTPUT spectrum of $\theta$ will be marked only by $\theta$. 


    
![svg](output_11_0.svg)
    



 The shape of the $\phi$ distribution is vastly different from $\theta$. Five features (comparison between $\theta$ and $\phi$) are noted and will be looked into more in detail:

 1. Dip in $\phi$ spectrum at low angles, roughly $0^{\circ} \leq \phi \leq 30^{\circ} $, we will call this **region 1**;
 2. Rise in $\phi$ spectrum at angles roughly $30^{\circ} \leq \phi \leq 70^{\circ} $, we will call this **region 2**;
 3. Dip in $\phi$ spectrum at angles roughly $70^{\circ} \leq \phi \leq 110^{\circ} $, we will call this **region 3**;
 4. Rise in $\phi$ spectrum at angles roughly $110^{\circ} \leq \phi \leq 150^{\circ} $, we will call this **region 4**;
 5. Dip in $\phi$ spectrum at angles roughly $150^{\circ} \leq \phi \leq 180^{\circ} $, we will call this **region 5**. 

The regions are marked in the following figure using separate color shadings. 

sdf2.thetaEscaped, 
nbins  = 40, 




    
![svg](output_13_0.svg)
    



 First of, it must be emphasized that the underlying principle between the change in shape from $\theta$ to $\phi$ overall is due to the fact that an event with $\theta$ as the initial decay angle between the electrons can (after the electrons traverse the source foil) change into **any physical** $\phi$. The probability of this change is not uniform, however. 

 For example, a single $\beta\beta$ decay event has $\theta = 10.5^{\circ}$. The electrons travel through the foil, scattering multiple times, and escape with $\phi = 110.5^{\circ}$. In such a case, in the $\theta$ distribution with bin-width $\Delta \phi = 1^{\circ}$, such an event would fall within the bin with edges $b_i \in (10, 11)^{\circ}$, whereas in the $\phi$ distribution it would be shifted to bin with edges $b_j \in (110, 111)^{\circ}$. 

 Of course, the initial $\theta = 10.5^{\circ}$ could have changed into any other physical angle (in range $(0,180)^{\circ}$), however with various probabilities. Quantitatively this effect will be studied later. Here we only focus on the qualitative effects of the detector. 

 With this in mind, there is another important factor to keep in mind. When $\theta$ is very small (the angle is *closed*) there is more *room* for the shift to be toward higher angles. When angle is *closed* it can only *open*. The opposite is true for the situation with large $\theta$ events. 

Qualitative analysis of $\phi$ distribution **by Region**
=========

It is noted, from the studies made in (DocDB #4816) that the track reconstruction is most successfull with electrons that escape the source foil perpendicular to the y-z plane of the source foil. On the other hand, the vertex reconstruction is decreased with events being parallel to the foil. 


## Region 1 and Region 5

The two *dips* in regions 1 (very small $\phi$) and 5 (very large $\phi$) can be explained with similar reasoning. First of all, the previosly mentioned *openning* (*closing*) of *closed* (*open*) angles results in events with $\theta$ within said regions being more likely to be shifted toward the other regions, thus these regions are more likely to be *drained* of their statistics.
Secondly, we can look at these regions in terms of how well events generated within them pass through the aforementioned data-cuts: 

 1. <font color="red">An event must have two negatively charged tracks reconstructed</font>, 
 2. <font color="red">Two vertices on the source foil, within given distance from each other</font>,
 3. <font color="green">Sum of electron energies within the range: $0~keV \leq E_{sum} \leq 3500~keV$</font>,
 4. <font color="red">Two individual Optical Module hits</font>,
 5. <font color="red">Two associated Optical Module hits</font>. 

Since $\phi$ is very small/large, both electrons will follow roughly the same path along a single line. Thus if one of the electrons escapes in the path **along** the source foil (which is difficult for track reconstruction due to the tracker cell design), *both* will travel along the foil. The first and second data-cut are not fulfilled. Furthermore, for events from region 1, it is possible that both electrons hit the same OM, thus fourth and fifth data-cuts are not fulfilled. Lastly, electrons with very small $\theta$ will fly along a very similar trajectory, which may be evaluated as a single particle in the reconstruction, rather than two separate electrons. The more apart the trajectories are from each other, the easier it is to reconstruct them as two individual particles. 

Adding these *difficulties* in the fulfilment of of data-cut condition, along with the shift between regions, results in the dips in the $\phi$ distribution shape. 

## Region 2 and Region 4

The $\phi$ distribution shows two rising peaks in the region 2 and region 4. These are regions with events of moderately opened/closed angles. In terms of the *shifting* of angles from regions 1 and 5, here most events will be shifted into, as it is the next closest (see quantitative analysis). Furthermore, from the geometry of the detector design, when one of the electrons escapes the foil perpendicular, the other electron will be within moderate angle off of it, however it will not be parallel to the foil, yet. These two regions are composed of electrons with *easiest* to reconstruct trajectories. 

Lastly, there is an apparent discrepancy between the height of the *peak* of region 2 and region 4 (with region 4 being the *taller* one). This discrepancy is likely due to the fact that - as noted in DocDB #4816 - with the smaller angles it is still probable that the reconstruction will not be successfull in separating the two individual tracks. 

 1. <font color="green">An event must have two negatively charged tracks reconstructed</font>, 
 2. <font color="green">Two vertices on the source foil, within given distance from each other</font>,
 3. <font color="green">Sum of electron energies within the range: $0~keV \leq E_{sum} \leq 3500~keV$</font>,
 4. <font color="green">Two individual Optical Module hits</font>,
 5. <font color="green">Two associated Optical Module hits</font>. 

## Region 3 

The *dip* in the region 3 is explained in DocDB #4816. The events with $\phi$ close to $90^{\circ}$ suffer from the detectors preference for events perpendicular to the source foil. When an electron escapes the source foil perpendicular to it, the other electron will escape parallel. Such events are then difficult for vertex reconstruction. 

 1. <font color="red">An event must have two negatively charged tracks reconstructed</font>, 
 2. <font color="red">Two vertices on the source foil, within given distance from each other</font>,
 3. <font color="green">Sum of electron energies within the range: $0~keV \leq E_{sum} \leq 3500~keV$</font>,
 4. <font color="green">Two individual Optical Module hits</font>,
 5. <font color="green">Two associated Optical Module hits</font>. 
 


Further investigation into region 2 and 4 discrepancy
=======

From analysis in DocDB #4816, it is easier to reconstruct events when the electrons escape on opposite sides of the foil, rather than same side. 



    
![svg](output_18_0.svg)
    



It can be seen that the vast majority of events in the region 4 originate from opposite side escape electrons, thus the better efficiency of reconstruction is apparent in the taller peak. 


### From the figure, it is visible that increasing the energy cut improves the corellation. 

Pending Tasks
===

* Check the frequency of how often electrons from region 1 hit the same OM.
* Check for each region, how often which data-cut is failed
* Check for absolute angle (w.r.t. the foil)

