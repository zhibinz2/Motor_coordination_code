%% generation of ARFI polynomial from AR polynomial and differencing parameter

function [A,Apol,G,C] = arfi2ar2(Cpol,d,q)

% clear; close all; clc;
% % NB: inserire un check per avere q>p!
% q=10; % truncating lag of differencing polynomial (high)
% % differencing parameter
% d=0.3;
% % AR polynomial
% ro0=0.65; %r0=0.65 (fr0=0)
% ro1=0.8; fr1=0.1; % ro1=0.8 fr1=0.1
% ro2=0.9; fr2=0.25; %ro2=0.92; fr2=0.25;
% Su=1; 
% z0=ro0;
% z1a=ro1*(cos(2*pi*fr1)+1i*sin(2*pi*fr1)); z1b=ro1*(cos(2*pi*fr1)-1i*sin(2*pi*fr1));
% z2a=ro2*(cos(2*pi*fr2)+1i*sin(2*pi*fr2)); z2b=ro2*(cos(2*pi*fr2)-1i*sin(2*pi*fr2));
% z=[z0; z1a; z1b; z2a; z2b];
% Cpol=poly(z);

M = size(Cpol,1);

%EACH CELL OF C CONTAINS ONE COEFFICIENT - MATRIX M*M
C = cell(0);
for i = 1:size(Cpol,2)/M-1
    C{i} = -Cpol(:,M*i+1:M*i+M);
end


if d==0 %CHANGED SIZES
    Apol = Cpol;
    A = -Apol(:,M+1:end);
    G = [eye(M) zeros(M,M*q)];
else
    p=length(C);
    if p>=q, error('this code works for differencing polynomial longer than AR polynomial'); end
    m=q+p;

    % coeffs of differencing part
    G = cell(0);
    for i=0:q
        %CHANGED TO DIAG
        %MATRIX M*M: Tsey, Maximum likelihood estimation of stationary
        %multivariate ARFIMA processes, 2010
        G{i+1}=diag(gamma(i-d)./(gamma(i+1).*gamma(-d)));
    end

    %CHANGED C WITH G TO HAVE C*G INSTEAD OF G*C
    
    Apol{1}=G{1}; %G(0)
    for k=1:p
        soma = zeros(M);
        for j=1:k
            soma = soma + C{j}*G{k+1-j};
        end
        Apol{k+1}=G{k+1}-soma;
    end

    for k=p+1:q
        soma = zeros(M);
        for j=1:p
            soma = soma + C{j}*G{k+1-j};
        end
        Apol{k+1}=G{k+1}-soma;
    end

    for k=q+1:m
        soma = zeros(M);
        for j=k-q:p
            soma = soma + C{j}*G{k+1-j};
        end
        Apol{k+1}=-soma;
    end

    A = -cat(2,Apol{2:end});
    Apol = cat(2,Apol{:});
    G = cat(2,G);

end

end

