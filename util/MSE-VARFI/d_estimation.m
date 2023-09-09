function d = d_estimation( x , d_min, d_max)
% I use d_min=-0.5; d_max=1;
% A. Leite, PhD Thesis page 80 

d = fminbnd(@(d) myfun(x,d),d_min,d_max);

end

function r=myfun(x,d)

% A. Leite, PhD Thesis page 80 

xl=reshape(x,[],1);
N=length(xl);
t=(1:N)';
m=round(N.^(0.5)); %justification in A. Leite, PhD Thesis page 80 
j=1:m;
wj=2*pi.*j/N;
r=log(1/m*sum(wj.^(2*d).*1/N*abs(xl.'*exp(-1i.*t*wj)).'.^2))-2*d/m*sum(log(wj));
end

