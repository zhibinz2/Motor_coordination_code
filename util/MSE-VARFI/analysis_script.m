%% EXEMPLARY SCRIPT FOR COMPUTING MULTISCALE ENTROPY (CE_x) BASED ON MULTIVARIATE ARFI MODELS
clear; close all; clc;

%% 1) open data and set analysis parameters
load('example_multivariate_series.mat')
Yo=[series_HP; series_SAP; series_RESP]; % all series in a matrix (series in row)

tau_vect = 1:5; % range of time scales to be explored
q_trunc = 50;    % lag at which truncate
pmax=12; % maximum  order for model order selection
ncoeff_FIR=48; % number of coeffs of FIR lowpass filter for linear MSE
d_min=-0.5; d_max=1; %raneg for d_estimation()

f_taus = 1./(2*tau_vect); %cutoff frequencies of rescaling filters
nscales=length(f_taus);  %number of time scales spanned

%% 2) estimate differencing parameter and filter each series to remove long-range correlations
M=size(Yo,1);
for m=1:M
    %d estimation (Whittle semiparametric estimator)
    est_d=d_estimation(Yo(m,:),d_min,d_max);
    ed(m)=est_d;
    %d filter
    [Ytmp,~]=remove_d(Yo(m,:)',est_d); %filtered data
    Y(m,:)=Ytmp';
end

% model identification
[~,p_ARd,~,~] = mos_idMVAR(Y,pmax,0); % model order selection (Bayesian Information criterion)
[eAm_ARd,eSu_ARd,~,~]=idMVAR(Y,p_ARd,0); % identification (ordinary least squares)


%% 3) Approximation of the VARFI process with a finite-order VAR process
eCpol=[eye(M) -eAm_ARd]; 
[eAm_ARFI,eApol] = arfi2ar2(eCpol,ed,q_trunc); 

%% 4) Computation of MSE for VAR process
for is=1:nscales
    clc; disp(['scale ' int2str(is) ' of ' int2str(nscales) '...'])
    tau=round(tau_vect(is));

    % LINEAR MSE on high-order VAR (computes both information storage IS and complexity IC)
    % core function: FIR filter for rescaling, VAR to VARMA to SS models,
    % downsampled SS process, estimation of partial variances as in [Faes et al Entropy 2017, ref. 19 of main paper]
    % and finally computation of global complexity (Eq. 3) and individual complexity (Eqs. 4a,b,c) for all target series
    [IS,IC]=ar2mse2(eAm_ARFI,eSu_ARd,tau,ncoeff_FIR);
    Complexity{is}=IC;
    
end



%% visualization
% Plot the time series
serie = {'HP','SBP','RESP'};  %labels
for m = 1:size(Yo,1)
    figure(1); subplot(3,1,m); 
    plot(Yo(m,:)); title([serie{m} ', original']); xlim([1 size(Yo,2)])
    figure(2); subplot(3,1,m); 
    plot(Y(m,:),'r'); title([serie{m} ', after removing d and normalization']); xlim([1 size(Yo,2)])
end

% retrieve and plot global complexity and complexity for heart rate (j=1)
j=1; i=2; k=3;
for is=1:nscales
    MSEall(is)=Complexity{is}.glob;
    MSE_j_j(is)=Complexity{is}.ind(j);
    MSE_j_ji(is)=Complexity{is}.ind_ij(j,i);
    MSE_j_jk(is)=Complexity{is}.ind_ij(j,k);
    MSE_j_jki(is)=Complexity{is}.ind_all(j);
end

figure;
subplot(1,2,1)
plot(MSEall,'-ok'); 
title('multivariate MSE')
subplot(1,2,2)
plot(MSE_j_j,'-ok'); hold on;
plot(MSE_j_ji,'-ob');
plot(MSE_j_jk,'-or');
plot(MSE_j_jki,'-og');
legend('C(j|{j})','C(j|{ji})','C(j|{jk})','C(j|{ijk})')
title(['MSE of series j=' int2str(j) '(i=' int2str(i) ',k=' int2str(k) ')'])
xlabel('\tau');
ylabel('CE');