function output1=DFA(DATA,win_length,order)
% output1 is the RMS (root mean squared value)
N=length(DATA);   % length of the data
n=floor(N/win_length); % the number of windows
N1=n*win_length; % length ot the data again

% initialize y
y=zeros(N1,1); % sum of deviations from the grand mean in all data length

% initialize fitcoef
fitcoef=zeros(n,order+1); % initialize fitcoef: the polynomial coeficiences for the fits

 % grand mean of the data
mean1=mean(DATA(1:N1)); % grand mean of the data

% sum of deviations from the grand mean in all data length
for i=1:N1      
    y(i)=sum(DATA(1:i)-mean1); 
end
y=y'; % transpose for polyfit

% polynomial coefficients in each window
for j=1:n % the number of windows
    fitcoef(j,:)=polyfit(1:win_length,y(((j-1)*win_length+1):j*win_length),order);
end

% initialize Yn (the new Y values in each window)
Yn=zeros(N1,1);
% the new Y values after fitting in each window
for j=1:n
    Yn(((j-1)*win_length+1):j*win_length)=polyval(fitcoef(j,:),1:win_length);
end

sum1=sum((y'-Yn).^2)/N1; % take the mean square
sum1=sqrt(sum1); % take the root
output1=sum1; % the final RMS (root mean squared value) 
% the bigger the output, the bigger the fluctuation is in this win_length scale