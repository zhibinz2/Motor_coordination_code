%% Computation of multiscale entropy and information storage from AR parameters
%  INPUT:  AR parameters Am, Su
%              scale factor tau
%              upsampling factor s, such that cutoff f_tau=s/(2*tau)

function [S_dws_n, C_dws_n, MSEdws,  MSEavg, Edws,Eavg, B0, Bm] = ar2mse2(Am,Su,tau,ncoeff)

M = size(Am,1); 
N = size(Am,2)/M;
s=1; %upsampling parameter of [Faes et al  2017]; not tested here

%%% upsampling
Am_cell = mat2cell(Am,M,repmat(M,[1 N])); %make a cell vector with the coefficients separated
Ams=[];
for iv=1:N
    Ams=[Ams zeros(M,M*(floor(s)-1)) Am_cell{iv}];
end
if isempty(Ams), Ams=0; end % to handle the case of white noise (AR process of order 0)      

% %%% stability check
M = size(Ams,1); eps=size(Ams,2)/M; 
Es=eye(M*eps);AAs=[Ams;Es(1:end-M,:)];lambda=eig(AAs);lambdamax=max(abs(lambda));
if lambdamax>=1% stability check
    error('The simulated VAR process is not stable');
end

%%% FILTERING
if tau==1
    q=0; b=1;
else
    q=ncoeff; % number of filter coeffs
    ft=1/(2*tau); %cutoff frequency
    Wn=2*ft; %normalized cutoff frequency (fNyquist=1)
    b=fir1(q,Wn,'noscale'); %Hamming window, linear phase (symmetry of b coeffs)
end
B0=b(1)*eye(M);
Bm=[];
for i=1:length(b)-1; Bm = [Bm, eye(M) * b(i+1)]; end

%ISS parameters
[A,C,K,V] = varma2iss(Ams,Bm,Su,B0); % max(abs(eig(A-K*C)))
%ret_flt = iss_PV(A,C,K,V);
%Eavg=0.5*log(2*pi*exp(1)*ret_flt.Sigmaj);
%MSEavg=0.5*log(2*pi*exp(1)*ret_flt.Sigmaj_j);
ret_flt = iss_PV2(A,C,K,V);
Eavg=0.5*log((2*pi*exp(1))^M*det(ret_flt.Sigmaj));
MSEavg=0.5*log((2*pi*exp(1))^M*det(ret_flt.Sigmaj_j));

%%% DOWNSAMPLING
[Ad,Kd,Vd] = iss_ds(A,C,K,V,tau);
Cd=C; %Rtau=R;

%% MSE and IS based on the determinant of Vd (Sigma_e) and Sigma_x
ret_dws_global = iss_PV2(Ad,Cd,Kd,Vd); 
% Edws=0.5*log(2*pi*exp(1)*ret_dws.Sigmaj);
% MSEdws=0.5*log(2*pi*exp(1)*ret_dws.Sigmaj_j);
% 
% Edws_n=0;%0.5*log(2*pi*exp(1)); 
% MSEdws_n=MSEdws-Edws+Edws_n; % MSE for normalized process

Edws = 0.5 * log(((2*pi*exp(1))^M)*det(ret_dws_global.Sigmaj)); %H_x
MSEdws = 0.5 * log(((2*pi*exp(1))^M)*det(ret_dws_global.Sigmaj_j)); %H_xx

Edws_n = 0.5*log((2*pi*exp(1))^M); 
S_dws_n.glob = Edws - MSEdws;
C_dws_n.glob = MSEdws - Edws + Edws_n; % equivalent: 0.5*log(((2*pi*exp(1))^M)*det(ret_dws_global.Sigmaj_j)/det(ret_dws_global.Sigmaj))


%% MSE and IS based on the partial covariances considering only the processes specified by jj - msID
%%%% complexity analysis of j given j
for jj = 1:M
    [VR, lambda0] = iss_PCOV_msID(Ad,Cd,Kd,Vd,jj);
    ret_dws_ind.Sigmaj(jj)=lambda0(jj,jj); %variance of series j
    ret_dws_ind.Sigmaj_j(jj)=VR; %partial variance of j given j
    
    Edws(jj) = 0.5 * log(2*pi*exp(1)*det(ret_dws_ind.Sigmaj(jj))); %entropy of j
    MSEdws(jj) = 0.5 * log(2*pi*exp(1)*det(ret_dws_ind.Sigmaj_j(jj))); %conditional entropy of j given j
    
    Edws_n = 0.5*log(2*pi*exp(1));
    S_dws_n.ind(jj) = Edws(jj) - MSEdws(jj); % information storage of j
    C_dws_n.ind(jj) = MSEdws(jj) - Edws(jj) + Edws_n; %complexity of j
end
    
%%%% complexity analysis of j given all (i,j,k)
VR = iss_PCOV_msID(Ad,Cd,Kd,Vd,1:M);
for jj=1:M
    ret_dws_ind.Sigmaj_ijk(jj)=VR(jj,jj); %partial variance of j given all (i,j,k)
    tmp = 0.5 * log(2*pi*exp(1)*det(ret_dws_ind.Sigmaj_ijk(jj))); %conditional entropy of j given all (i,j,k)
    
    Edws_n = 0.5*log(2*pi*exp(1));
    S_dws_n.ind_all(jj) = Edws(jj) - tmp; % information storage of j given all (i,j,k)
    C_dws_n.ind_all(jj) = tmp - Edws(jj) + Edws_n; %complexity of j given all (i,j,k)
end

%%%% complexity analysis of j given i,j (target and one driver)
for jj = 1:M
    for ii=1:M
        if ii~=jj
            VR = iss_PCOV_msID(Ad,Cd,Kd,Vd,[jj ii]);
            ret_dws_ind.Sigmaj_ij(jj,ii)=VR(1,1); %partial variance of j given i,j
            tmp = 0.5 * log(2*pi*exp(1)*det(ret_dws_ind.Sigmaj_ij(jj,ii))); %conditional entropy of j given i,j
            
            Edws_n = 0.5*log(2*pi*exp(1));
            S_dws_n.ind_ij(jj,ii) = Edws(jj) - tmp; % information storage of j given i,j
            C_dws_n.ind_ij(jj,ii) = tmp - Edws(jj) + Edws_n; %complexity of j given i,j
            
        end
    end
end





end


