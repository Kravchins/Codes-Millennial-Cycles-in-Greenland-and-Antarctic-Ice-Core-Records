# Codes-Millennial-Cycles-in-Greenland-and-Antarctic-Ice-Core-Records
MATLAB Codes for Millennial-Scale Spectral Analysis of Paleoclimate Proxy Data and LODE algorithm for calculating the ice cover and temperature

(1) Synchrosqueezing Transform and Ridge Extraction for paleoclimate records; (2) Logistic delay-differential equation (LODE) algorithm for calculating the ice cover and temperature.

Any use of this code MUST refer to the original publication: 

Kravchinsky,V.A., Zhang,R., Borowiecki,R., Goguitchaichvili,A., Czarnecki,J., Czarnecki,A., Boers,N., Berger,A., & van der Baan,M. (2025). Millennial Cycles in Greenland and Antarctic Ice Core Records: Evidence of Astronomical Influence on Global Climate. Journal of Geophysical Research: Atmospheres. (the publication doi will be added here) https://doi.org/

and

Borowiecki, R., Kravchinsky, V.A., van der Baan, M., & Herrera, R.H. (2023). The synchrosqueezing transform to evaluate paleoclimate cyclicity. Computers & Geosciences 175, 105336. https://doi.org/10.1016/j.cageo.2023.105336

Questions/suggestions to vadim@ualberta.ca

This script was tested to run in Matlab R2024a with dependencies on the wavelet toolbox and SST toolbox. Download the folder "synchrosqueezing" from https://github.com/ebrevdo/synchrosqueezing/tree/master and place this folder in your folder with the codes and input data.

The main script loads data, performs the necessary computations, and plots out the results of both the ridge extraction and the windowed reconstruction.

Please ensure the main script, along with both functions (coi and red noise), are added to your working path, along with the XLS file containing sample synthetic data.
