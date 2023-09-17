clear
addpath(genpath('./util'))
run color_scheme.m

data_path= '../Cleaned_data/'; % Cleaned EEG data can be found at https://osf.io/rstpu/. 
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

% Compute and extract H
H_all=[];
for s=1:numSes
    clear intervals
    runid=num2str(seeds(s,:));
    load([data_path 'clean_' runid '.mat'],'intervals')   
    for b=1:12
        [~,H_all(1,s,b)]=DFA_main(intervals{b}(:,1));
        [~,H_all(2,s,b)]=DFA_main(intervals{b}(:,2));
    end
end

% organize H_all 
H_all; % 2x12numSesx12trials
H_Lall=squeeze(H_all(1,:,:)); % 12session by 12trials
H_Lall=reshape(H_Lall',[],1); % in time sequence
H_Rall=squeeze(H_all(2,:,:)); % 12session by 12trials
H_Rall=reshape(H_Rall',[],1); % in time sequence
H_all_LR=[H_Lall; H_Rall];

% organize conditions
condition_all=[];
for s=1:numSes
    clear conditions
    runid=num2str(seeds(s,:));
    load([data_path 'clean_' runid '.mat'],'conditions');
    condition_all(s,:)=conditions; % 12 session x 12 trials
end
% reshape into a vector in time sequence
condition_all=reshape(condition_all',[],1); % 144 x 1 

% indices for 4 states
states4names={'Independent','Leading','Following','Bidrectional'};
% find the indices for each condition in L R conbined sequence (4 states)
uncoupleInd_LR=[find(condition_all==1);12*numSes+find(condition_all==1)];
leadingInd_LR=[find(condition_all==2);12*numSes+find(condition_all==3)];
followingInd_LR=[find(condition_all==3);12*numSes+find(condition_all==2)];
mutualInd_LR=[find(condition_all==4);12*numSes+find(condition_all==4)];
Inds4_LR=[uncoupleInd_LR leadingInd_LR followingInd_LR mutualInd_LR];

% indicies in the synch and synco time sequence
synchind=[1:3 7:9 13:15 19:21 25:27 31:33]; % 3 trials x 6 sessions
syncoind=[4:6 10:12 16:18 22:24 28:30 34:36]; % 3 trials x 6 sessions
synind=[synchind;syncoind];

%% Plot Fig 4
canvas(0.15,1);

subplot(3,1,1); 
plot(H_Lall(uncoupleInd_LR(synchind(1:3))),H_Rall(uncoupleInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen);hold on;
plot(H_Lall(uncoupleInd_LR(synchind(4:6))),H_Rall(uncoupleInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(7:9))),H_Rall(uncoupleInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(10:12))),H_Rall(uncoupleInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(13:15))),H_Rall(uncoupleInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(synchind(16:18))),H_Rall(uncoupleInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(uncoupleInd_LR(syncoind(1:3))),H_Rall(uncoupleInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(4:6))),H_Rall(uncoupleInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(7:9))),H_Rall(uncoupleInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(10:12))),H_Rall(uncoupleInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(13:15))),H_Rall(uncoupleInd_LR(syncoind(13:15))),'*','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(uncoupleInd_LR(syncoind(16:18))),H_Rall(uncoupleInd_LR(syncoind(16:18))),'diamond','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
ax = gca;
ax.FontSize = 13; 
xlabel('DFA Hurst Exponent, Participant A','FontSize',12);
ylabel('DFA Husrt Exponent, Participant B','FontSize',12);
title('Independent','FontSize',13);

% fit the combined
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(uncoupleInd_LR(1:36)),H_Rall(uncoupleInd_LR(1:36)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); 
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'k-','LineWidth',3);
text(max(xxx)+0.01,min(FitValues),sprintf('\\rho=%.2f',RHO),'Color',[0 0 0],'FontSize', 13)

% fit the synch data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(uncoupleInd_LR(1:2:35)),H_Rall(uncoupleInd_LR(1:2:35)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1);
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',darkgreen,'LineWidth',3);
text(max(xxx),min(FitValues),sprintf('\\rho=%.2f',RHO),'Color',darkgreen,'FontSize', 13)

% fit the synco data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(uncoupleInd_LR(2:2:36)),H_Rall(uncoupleInd_LR(2:2:36)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1);
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',pink,'LineWidth',3);
xlim([0.2 1.4]);ylim([0.2 1.4]);
plot([0.5 1], [0.5 1],'m--');
hold off;
text(max(xxx),max(FitValues)+0.01,sprintf('\\rho=%.2f',RHO),'Color',pink,'FontSize', 13)
grid on;

lg=legend('synch - subj pair 1','synch - subj pair 2','synch - subj pair 3','synch - subj pair 4', ...
    'synch - subj pair 5','synch - subj pair 6', ...
    'synco - subj pair 1','synco - subj pair 2','synco - subj pair 3','synco - subj pair 4',...
    'synco - subj pair 5','synco - subj pair 6','location','eastoutside');


subplot(3,1,2); 
plot(H_Lall(leadingInd_LR(synchind(1:3))),H_Rall(leadingInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen); hold on;
plot(H_Lall(leadingInd_LR(synchind(4:6))),H_Rall(leadingInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(7:9))),H_Rall(leadingInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(10:12))),H_Rall(leadingInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(13:15))),H_Rall(leadingInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(synchind(16:18))),H_Rall(leadingInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(leadingInd_LR(syncoind(1:3))),H_Rall(leadingInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink);
plot(H_Lall(leadingInd_LR(syncoind(4:6))),H_Rall(leadingInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(7:9))),H_Rall(leadingInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(10:12))),H_Rall(leadingInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(13:15))),H_Rall(leadingInd_LR(syncoind(13:15))),'*','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(leadingInd_LR(syncoind(16:18))),H_Rall(leadingInd_LR(syncoind(16:18))),'diamond','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(synchind(1:3))),H_Lall(followingInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen);
plot(H_Rall(followingInd_LR(synchind(4:6))),H_Lall(followingInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(7:9))),H_Lall(followingInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(10:12))),H_Lall(followingInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(13:15))),H_Lall(followingInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(synchind(16:18))),H_Lall(followingInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Rall(followingInd_LR(syncoind(1:3))),H_Lall(followingInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink);
plot(H_Rall(followingInd_LR(syncoind(4:6))),H_Lall(followingInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(7:9))),H_Lall(followingInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(10:12))),H_Lall(followingInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(13:15))),H_Lall(followingInd_LR(syncoind(13:15))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Rall(followingInd_LR(syncoind(16:18))),H_Lall(followingInd_LR(syncoind(16:18))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
ax = gca;
ax.FontSize = 13; 
xlabel('DFA Hurst Exponent, Leader','FontSize',12);
ylabel('DFA Hurst Exponent, Follower','FontSize',12);
title('Unidirectional','FontSize',13);

% fit the combined data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues([H_Lall(leadingInd_LR(1:36)); H_Rall(followingInd_LR(1:36))], ...
    [H_Rall(leadingInd_LR(1:36)); H_Lall(followingInd_LR(1:36))],0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); 
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'k-','LineWidth',3);
text(max(xxx)+0.01,max(FitValues),sprintf('\\rho=%.2f',RHO),'Color',[0 0 0],'FontSize', 13)

% fit the synch data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues([H_Lall(leadingInd_LR(1:2:35)); H_Rall(followingInd_LR(1:2:35))], ...
    [H_Rall(leadingInd_LR(1:2:35)); H_Lall(followingInd_LR(1:2:35))],0.2,1.4);
A=polyfit(xxx, yyy,1);
[RHO,~]=corr(xxx, yyy);
Alpha1=A(1); 
FitValues=polyval(A,xxx);
hold on; 
plot(xxx,FitValues,'-','color',darkgreen,'LineWidth',3);
text(max(xxx)-0.1,max(FitValues)-0.03,sprintf('\\rho=%.2f',RHO),'Color',darkgreen,'FontSize', 13)

% fit the synco data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues([H_Lall(leadingInd_LR(2:2:36)); H_Rall(followingInd_LR(2:2:36))],...
    [H_Rall(leadingInd_LR(2:2:36)); H_Lall(followingInd_LR(2:2:36))],0.2,1.4);
A=polyfit(xxx, yyy,1);
[RHO,~]=corr(xxx, yyy);
Alpha1=A(1); 
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',pink,'LineWidth',3);
text(max(xxx)-0.1,max(FitValues)+0.1,sprintf('\\rho=%.2f',RHO),'Color',pink,'FontSize', 13)
xlim([0.2 1.4]);ylim([0.2 1.4]);
plot([0.5 1], [0.5 1],'m--');
plot([0.5 1], [0.5 1],'m--');
hold off;
grid on;

lg=legend('synch - subj pair 1','synch - subj pair 2','synch - subj pair 3','synch - subj pair 4', ...
    'synch - subj pair 5','synch - subj pair 6', ...
    'synco - subj pair 1','synco - subj pair 2','synco - subj pair 3','synco - subj pair 4',...
    'synco - subj pair 5','synco - subj pair 6','location','eastoutside');


subplot(3,1,3);
plot(H_Lall(mutualInd_LR(synchind(1:3))),H_Rall(mutualInd_LR(synchind(1:3))),'.','MarkerSize',30,'color',darkgreen);hold on;
plot(H_Lall(mutualInd_LR(synchind(4:6))),H_Rall(mutualInd_LR(synchind(4:6))),'square','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(7:9))),H_Rall(mutualInd_LR(synchind(7:9))),'^','MarkerSize',8,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(10:12))),H_Rall(mutualInd_LR(synchind(10:12))),'pentagram','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(13:15))),H_Rall(mutualInd_LR(synchind(13:15))),'*','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(synchind(16:18))),H_Rall(mutualInd_LR(synchind(16:18))),'diamond','MarkerSize',10,'color',darkgreen,'MarkerFaceColor',darkgreen);
plot(H_Lall(mutualInd_LR(syncoind(1:3))),H_Rall(mutualInd_LR(syncoind(1:3))),'.','MarkerSize',30,'color',pink);
plot(H_Lall(mutualInd_LR(syncoind(4:6))),H_Rall(mutualInd_LR(syncoind(4:6))),'square','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(7:9))),H_Rall(mutualInd_LR(syncoind(7:9))),'^','MarkerSize',8,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(10:12))),H_Rall(mutualInd_LR(syncoind(10:12))),'pentagram','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(13:15))),H_Rall(mutualInd_LR(syncoind(13:15))),'*','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
plot(H_Lall(mutualInd_LR(syncoind(16:18))),H_Rall(mutualInd_LR(syncoind(16:18))),'diamond','MarkerSize',10,'color',pink,'MarkerFaceColor',pink);
ax = gca;
ax.FontSize = 13; 
xlabel('DFA Hurst Exponent, Participant A','FontSize',12);
ylabel('DFA Husrt Exponent, Participant B','FontSize',12);
title('Bidirectional','FontSize',13);
xlim([0.2 1.4]);ylim([0.2 1.4]);

% fit all the data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(mutualInd_LR(1:36)),...
    H_Rall(mutualInd_LR(1:36)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); 
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'k-','LineWidth',3);
text(max(xxx)-0.2,max(FitValues)-0.2,sprintf('\\rho=%.2f',RHO),'Color',[0 0 0],'FontSize', 13)

% fit the synch data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(mutualInd_LR(1:2:35)),...
    H_Rall(mutualInd_LR(1:2:35)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); 
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',darkgreen,'LineWidth',3);
text(max(xxx)-0.1,max(FitValues)-0.1,sprintf('\\rho=%.2f',RHO),'Color',darkgreen,'FontSize', 13)

% fit the synco data
A=[];Alpha1=[];FitValues=[];RHO=[];xxx=[];yyy=[];
[xxx,yyy]=keepvalues(H_Lall(mutualInd_LR(2:2:36)),...
    H_Rall(mutualInd_LR(2:2:36)),0.2,1.4);
A=polyfit(xxx,yyy,1);
[RHO,~]=corr(xxx,yyy);
Alpha1=A(1); 
FitValues=polyval(A,xxx);
hold on; plot(xxx,FitValues,'-','color',pink,'LineWidth',3);
text(max(xxx)-0.5,max(FitValues)-0.25,sprintf('\\rho=%.2f',RHO),'Color',pink,'FontSize', 13)
xlim([0.2 1.4]);ylim([0.2 1.4]);
plot([0.5 1], [0.5 1],'m--');
hold off;
grid on;


lg=legend('synch - subj pair 1','synch - subj pair 2','synch - subj pair 3','synch - subj pair 4', ...
    'synch - subj pair 5','synch - subj pair 6', ...
    'synco - subj pair 1','synco - subj pair 2','synco - subj pair 3','synco - subj pair 4',...
    'synco - subj pair 5','synco - subj pair 6','location','eastoutside');
% lg.Position = [0.9475 0.15 0.01 0.15];

delete(findall(gcf,'type','annotation'))
set(gcf,'color','w'); 

