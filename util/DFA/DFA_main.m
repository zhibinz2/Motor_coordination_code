function [D,Alpha1,n,F_n,FitValues]=DFA_main(DATA)

n=5:5:length(DATA)/2;
Nw=length(n); % number of different win_lengths
F_n=zeros(Nw,1); % initialize RMS values in the FDA time series
 for i=1:Nw
     F_n(i)=DFA(DATA,n(i),1);
 end
n=n';
 
A=polyfit(log10(n),log10(F_n),1);
Alpha1=A(1); 
D=3-A(1); 

% plot the fit
FitValues=polyval(A,log10(n(1:end)));

return
