DCE-MRI:
Dual input single compartment model for liver perfusion mapping 

By setting the input parameters below, the algorithm generates tissue
enhancement curve, adds complex noise to dynamic curves, solves for
perfusion parameters using whole-field linear least-squares method, and 
generates maps of relative error.  

Please refer below for more info:
Jafari, R., Chhabra, S., Prince, M. R., Wang, Y. and Spincemaille, P. (2017),
Vastly accelerated linear least-squares fitting with numerical optimization
for dual-input delay-compensated quantitative liver perfusion mapping. Magn.
Reson. Med. doi:10.1002/mrm.26888

Contact: rj259@cornell.edu

Set Input Parameters:
 
af_thr = 0.9; %arterial fraction

ecv_thr = 0.8; % extracellular volume

mtt_thr = 10; %mean transit time

m = 100; % generate a grid with m by m pixels

cnr = 1e9; %contrast-to-noise ratio

