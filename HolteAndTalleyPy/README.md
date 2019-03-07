# Holte and Talley Mixed Layer Depth algorithm

The paper outlining this algorithm is [here](http://http://mixedlayer.ucsd.edu/data/HolteTalley2009.pdf)


## The Temperature Algorithm
### Summer Flow Graph
![summer flow graph](readmeimages/summer.png)
### Winter Flow Graph
![winter flow graph](readmeimages/winter.png)

## The Density and Salinity Algorithms
The original paper does not include flowcharts for the density and salinity algorithms so I have included those below transcribed from the matlab supplementary materials provided.

![](readmeimages/salinitysummer.jpg)
![](readmeimages/salinitywinter.jpg)
![](readmeimages/densitysummer.jpg)
![](readmeimages/densitywinter.jpg)

## Glossary of Terms

## MLTFit

The MLTFit is defined as the interection of the line of best fit of the thermocline and the mixed layer. The best fit line of the thermocline is the line tangent to the profile at the point where d(T)/d(P) is the greatest. To determine the best fit line of the mixed layer, the algorithm calculates the best fit over progressively more points in the profile, calculates the error over the points included in the fit and then finds the deepest fit line which satisfies an error tolerance E<sub>T</sub> < 10<sup>-10</sup>.

## TTMLD

The deepest point in the profile which is less than 0.2C compared to the SST

## ΔT

T(i<sub>MLTFIT</sub>)-T(i<sub>MLTFIT+2</sub>)

## range

25 dbar

## TDTM

The minimum of the gradient maximum and the temperature maximum if that is less than 100 dbar. Otherwise it is set to 0

## DTM

Deepest point such that the gradient is less than 0.005 C/dbar or max if that criteria is not met.


## Some results

![](readmeimages/summerprofile.png)
![](readmeimages/winterprofile.png)
![](readmeimages/h&Texamplerun.png)
