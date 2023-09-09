% Compute partial variances for a state space model from innovations form parameters 
% assumes one target and one source, and iterates
% A,C,K,V - innovations form state space parameters

function ret = iss_PV2(A,C,K,V)

% determine the variance of the processes lambda0=E[Yn Yn'] for each target Y
O=dlyap(A,K*V*K'); %EQ.16A LMSE
lambda0=C*O*C'+V;  %EQ.16B LMSE
Sj = lambda0;%diag(lambda0);
Sj_j = V; %Sj_j- = S_j - Sj,j- * Sj-^(-1) * Sj,j-^T  = V ? eq.1 e 3 Seth and Barret (2009)

ret.Sigmaj=Sj;
ret.Sigmaj_j=Sj_j;

end

