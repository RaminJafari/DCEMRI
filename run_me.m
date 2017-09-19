close all;
clear;

%By setting the input parameters below, the algorithm generates tissue
%enhancement curve, adds complex noise to dynamic curves, solves for
%perfusion parameters using whole-field linear least-squares method, and 
%generates maps of relative error.  

%Please refer below for more info:
%Jafari, R., Chhabra, S., Prince, M. R., Wang, Y. and Spincemaille, P. (2017),
%Vastly accelerated linear least-squares fitting with numerical optimization
%for dual-input delay-compensated quantitative liver perfusion mapping. Magn.
%Reson. Med. doi:10.1002/mrm.26888
%Contact: rj259@cornell.edu

%%%%%%%%%%%% Set Input Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
af_thr = 0.9; %arterial fraction
ecv_thr = 0.8; % extracellular volume
mtt_thr = 10; %mean transit time
m = 100; % generate a grid with m by m pixels
cnr = 1e9; %contrast-to-noise ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('inputs.mat') %load experimental arterial and portal venous enahancemnet curves
t_s = 0.2298; %sampling time
time_sim_thr = (1:size(aif,1))'*t_s; % generates time axis for enhancement curves

k2 = 1/mtt_thr;
ka = af_thr*(ecv_thr*k2);
kp = ecv_thr*k2 - ka;
coeff = [ka, kp, k2];
curve = [aif, pvif ,time_sim_thr];     
enhancement_sim_thr = dualinput_nlls(coeff,curve);

time_sim_mat = repmat(time_sim_thr,1,m^2);
aif_mat_thr = repmat(aif,1,m^2);
pvif_mat_thr = repmat(pvif,1,m^2);
enhancement_mat_thr = repmat(enhancement_sim_thr,1,m^2);
%%
%generate dynamic curves and add complex noise
pre = 500;
dynamic_mat_thr = enhancement_mat_thr*pre + pre;
dynamic_aif_mat_thr = aif_mat_thr*pre + pre;
dynamic_pvif_mat_thr = pvif_mat_thr*pre + pre;

sigma = max(squeeze(dynamic_mat_thr(1,1,:)))/cnr;
mu = 0;

dynamic_mat_thr_noise = abs(dynamic_mat_thr+complex(normrnd(mu,sigma,size(dynamic_mat_thr)), normrnd(mu,sigma,size(dynamic_mat_thr))));
enhancement_mat_thr_noise = (dynamic_mat_thr_noise-pre)./pre;
dynamic_aif_mat_thr_noise = abs(dynamic_aif_mat_thr+complex(normrnd(mu,sigma,size(dynamic_mat_thr)), normrnd(mu,sigma,size(dynamic_mat_thr))));
aif_mat_thr_noise = (dynamic_aif_mat_thr_noise-pre)./pre;
dynamic_pvif_mat_thr_noise = abs(dynamic_pvif_mat_thr+complex(normrnd(mu,sigma,size(dynamic_mat_thr)), normrnd(mu,sigma,size(dynamic_mat_thr))));
pvif_mat_thr_noise = (dynamic_pvif_mat_thr_noise-pre)./pre;

enhancement_mat_thr_noise(enhancement_mat_thr_noise<0) = 1e-6;
aif_mat_thr_noise(aif_mat_thr_noise < 0) = 1e-6;
pvif_mat_thr_noise(pvif_mat_thr_noise < 0) = 1e-6;
%%
%calculate perfusion parameters from enhancement curves using whole-field linear least squares fit
enhancement_mat = reshape(enhancement_mat_thr_noise(:),[size(enhancement_mat_thr_noise,1)*size(enhancement_mat_thr_noise,2),size(enhancement_mat_thr_noise,3)]);
aif_mat = reshape(aif_mat_thr(:),[size(aif_mat_thr_noise,1)*size(aif_mat_thr_noise,2),size(aif_mat_thr_noise,3)]);
pvif_mat = reshape(pvif_mat_thr_noise(:),[size(pvif_mat_thr_noise,1)*size(pvif_mat_thr_noise,2),size(pvif_mat_thr_noise,3)]);

[A,B] = cg_A_B_mats(aif_mat_thr_noise,pvif_mat_thr_noise,enhancement_mat_thr_noise,time_sim_mat);
[g,h] = size(A);
reg = 0.00004*(speye(g,h));
[perf,flag_per,relres,iter,resvec] = cgs(A+reg,B,1e-12,10000);
%%
%allocate perfusion paramters to corresponsing voxel in the image
add = -1;
for i = 1:m
 for j = 1:m      
  add =add+1;
  cnt = ((add)*3)+1;
  model = [perf(cnt);perf(cnt+1);perf(cnt+2)];
  model = model./t_s;
  
  mtt (i,j) = -1./(model(3));
  af(i,j) =  model(1)/(model(1)+model(2));
  ecv(i,j) = (model(1)+model(2))/( -1*model(3));
 end
end
%%
%calculate relative error and plot results

ecv_err = 100*(abs(ecv-ecv_thr)./ecv_thr);
af_err = 100*(abs(af-af_thr)./af_thr);
mtt_err = 100*(abs(mtt-mtt_thr)./mtt_thr);

figure
subplot(221)
imagesc(af_err)
title('arterial fraction relative error (%)')
colorbar
pbaspect([1 1 1])

subplot(222)
imagesc(ecv_err)
title('extracellular volume relative error (%)')
colorbar
pbaspect([1 1 1])

subplot(223)
imagesc(mtt_err)
title('mean transit time relative error(%)')
colorbar
pbaspect([1 1 1])

subplot(224)
plot(time_sim_thr,aif)
hold on
plot(time_sim_thr,pvif)
hold on
plot(time_sim_thr,enhancement_sim_thr)
legend('aif','pvif','tissue')
xlabel('time (sec)')
ylabel('noiseless relative enhancement')
pbaspect([1 1 1])
