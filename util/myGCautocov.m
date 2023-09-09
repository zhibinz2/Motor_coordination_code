function [F,pval,sig] = myGCautocov(X);
% this function runs Granger Causality based on MVGC's autocov method
% Parameters
ntrials   = size(X,3);     % number of trials
nobs      = size(X,2);  % number of observations per trial
regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default) 
morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 10;     % maximum model order for model order estimation
% recomended acmaxlags minimum = 105
acmaxlags = 10;   % maximum autocovariance lags (empty for automatic calculation)
tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
% for plotting only
fs        = 1.3;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)
nvars = size(X,1); % number of variables
% model order estimation
[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
% Select model order.
if     strcmpi(morder,'actual')
    morder = amo;
    fprintf('\nusing actual model order = %d\n',morder);
elseif strcmpi(morder,'AIC')
    morder = moAIC;
    fprintf('\nusing AIC best model order = %d\n',morder);
elseif strcmpi(morder,'BIC')
    morder = moBIC;
    fprintf('\nusing BIC best model order = %d\n',morder);
else
    fprintf('\nusing specified model order = %d\n',morder);
end
% VAR model estimation 
[A,SIG] = tsdata_to_var(X,morder,regmode);
% Check for failed regression
assert(~isbad(A),'VAR estimation failed');
% Autocovariance calculation
[G,info] = var_to_autocov(A,SIG,acmaxlags);
var_acinfo(info,true); % report results (and bail out on error)
% Granger causality calculation: time domain 
F = autocov_to_pwcgc(G);
assert(~isbad(F,false),'GC calculation failed');
pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);
end

