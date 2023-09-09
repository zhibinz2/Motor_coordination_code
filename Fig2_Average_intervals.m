clear
addpath(genpath('./util'))
run color_scheme.m

data_path= '../Cleaned_data/'; % Cleaned EEG data can be found at https://osf.io/rstpu/. 
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

% load conditions
condition_all=[];
for s=1:numSes
    clear conditions
    runid=num2str(seeds(s,:));
    load([data_path 'clean_' runid '.mat'],'conditions');
    condition_all(s,:)=conditions; % 12 session x 12 trials
end
% reshape into a vector in time sequence
condition_all=reshape(condition_all',[],1); % 144 x 1 

% load tapping intervals
intervals_L={}; % 12 session x 12 trials
intervals_R={}; % 12 session x 12 trials
for s=1:numSes
    runid=num2str(seeds(s,:));
    clear intervals
    load([data_path 'clean_' runid '.mat'],'intervals');
    for b=1:12
        intervals_L{s,b}=intervals{b}(:,1); % L subject 
        intervals_R{s,b}=intervals{b}(:,2); % R subject 
    end
end

% indicies in the synch and synco time sequence
synchind=[1:3 7:9 13:15 19:21 25:27 31:33]; % 3 trials x 6 sessions
syncoind=[4:6 10:12 16:18 22:24 28:30 34:36]; % 3 trials x 6 sessions
synind=[synchind;syncoind];


% indices for 4 states
states4names={'Independent','Leading','Following','Bidrectional'};
% find the indices for each condition in L R conbined sequence (4 states)
uncoupleInd_LR=[find(condition_all==1);12*numSes+find(condition_all==1)];
leadingInd_LR=[find(condition_all==2);12*numSes+find(condition_all==3)];
followingInd_LR=[find(condition_all==3);12*numSes+find(condition_all==2)];
mutualInd_LR=[find(condition_all==4);12*numSes+find(condition_all==4)];
Inds4_LR=[uncoupleInd_LR leadingInd_LR followingInd_LR mutualInd_LR];

% reorganize intervals_L and intervals_R;
% reorganize the cell sequence for intervals_all_L and intervals_all_R;
intervals_all_L=reshape(intervals_L',[],1);
intervals_all_R=reshape(intervals_R',[],1);
intervals_all_LR_all=[intervals_all_L; intervals_all_R]; % combine cell array

% organize intervals for 4 states in 2 syn types
syn2names={'Synch','Synco'};
intervals_2t_4ss={}; % 2 x 4 cells
for t=1:2
    for ss=1:4
        cat_stateIntervals=[];
        cat_stateIntervals=cat(1,intervals_all_LR_all{Inds4_LR([synind(t,:) 36+synind(t,:)],ss)});
        intervals_2t_4ss{t,ss}=cat_stateIntervals;
    end
end


% transpose to (4 states x 2 groups)
intervals_2t_4ss=intervals_2t_4ss';
% compute mean and std (4 states x 2 groups)
intervals_2t_4ss_mean=[];
intervals_2t_4ss_std=[];
% intervals_2t_4ss_ste=[]; % standard error
intervals_2t_4ss_25=[]; % percentile
intervals_2t_4ss_75=[]; % percentile
for t=1:2
    for ss=1:4
        intervals_2t_4ss_mean(ss,t)=mean(intervals_2t_4ss{ss,t}./2); % and convert to ms
        intervals_2t_4ss_std(ss,t)=std(intervals_2t_4ss{ss,t}./2);
    end
end

%% plot Fig2
% barplot with errorbar (published)
canvas(0.23, 0.4);
model_series = intervals_2t_4ss_mean;
model_error = intervals_2t_4ss_std;
b = bar(model_series, 'grouped');
b(1).FaceColor=darkgreen;b(2).FaceColor=pink;
hold on;
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars 
errorbar(x',model_series,model_error,'k','linestyle','none','LineWidth',2);

xticks(1:4);xticklabels({'Independent','Leader','Follower','Bidirectional'});
xl = get(gca,'XTickLabel');  
set(gca,'XTickLabel',xl,'fontsize',17,'FontWeight','bold')
xlim([0.5 4.5]);ylim([400 1000]);
ylabel('Mean tapping interval (ms)','FontSize',15); 
yl = get(gca,'YTickLabel');  
set(gca,'YTickLabel',yl,'fontsize',17,'FontWeight','bold')
set(gcf,'color','w'); 
lg=legend({'Synch','Synco'},'location','north');lg.FontSize=17;