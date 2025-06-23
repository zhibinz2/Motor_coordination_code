clear
addpath(genpath('./util'))
run color_scheme.m

data_path= '../Cleaned_data/'; % Cleaned EEG data can be found at https://osf.io/rstpu/. 
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

% d estimation from H (DFA method) based on the good matched intervals
Xcorr10Lag=nan(numSes,12,21);
for r=1:numSes;
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    load([data_path  'clean_' runid '.mat'],'conditions','intervals');
    % sort order and plot results
    [x,sortorder]=sort(conditions);
    for j = 1:12
        clear intL_good_dmean intR_good_dmean
        % remove the mean
        intL_good_dmean=intervals{sortorder(j)}(:,1)-mean(intervals{sortorder(j)}(:,1));
        intR_good_dmean=intervals{sortorder(j)}(:,2)-mean(intervals{sortorder(j)}(:,2));
        % Xcorr based on int_dmean_drm (before d removal)********
        r12=[];lags12=[];
        [r12,lags12]=xcorr(intL_good_dmean,intR_good_dmean,10,'normalized');
        Xcorr10Lag(r,j,1:21)=r12;
    end
end

syn2Ind={[1:2:11],[2:2:12]};syn2names={'Synch','Synco'};


%% Plot Fig6
Xcorr10Lag; % sorted order
direction3names={'Independent','Unidirectional','Bidirectional'};
sorted3inds={[1:3],[4:9],[10:12]};
sorted4inds=[1:3; 4:6; 7:9; 10:12];
syn2Ind;
canvas(0.13, 0.9); % figure;
for i=1:3
    for syn=1:2
    subplot(3,1,i)
    if i==2;
        % L-lead
        L_mat=squeeze(mean(Xcorr10Lag(syn2Ind{syn},sorted4inds(2,:),:),2));
        % R-lead
        R_mat=squeeze(mean(Xcorr10Lag(syn2Ind{syn},sorted4inds(3,:),:),2));
        % combine
        LR_mean=mean([L_mat;fliplr(R_mat)]);
        % plot Leader vs Follower
        hold on;
        plot(-10:1:10,LR_mean,'linewidth',4,'color',[syn2colors(syn,:) 1]);
        ax = gca;
        ax.FontSize = 13; 
        title(direction3names{i},'FontSize',13);
        ylabel('\rho','FontSize',15,'FontWeight','bold');
        xlabel({'Lag','    Leader -leading <- 0 -> Follower - following'},'FontSize',13);ylim([-0.1 0.8])
    elseif i==1 | i==3;
        hold on;
        plot(-10:1:10,mean(squeeze(mean(Xcorr10Lag(syn2Ind{syn},sorted3inds{i},:),2)),1),'linewidth',4,'color',[syn2colors(syn,:) 1]);
        ax = gca;
        ax.FontSize = 15; 
        title(direction3names{i},'FontSize',13);
        ylabel('\rho','FontSize',15,'FontWeight','bold');xlabel({'Lag', 'A -leading <- 0 -> B -leading'},'FontSize',13);ylim([-0.1 0.8])
    end    
    hold off;
    end
    yline(0,'color',[1 0.8 0.2]);xline(0,'color',[1 0.8 0.2]);
    xline(-1,'color',[1 0.8 0.2]);xline(1,'color',[1 0.8 0.2]);
    lg=legend({'synch','synco','','','',''});
    lg.FontSize=17;
end

set(gcf,'color','w'); % set background white for copying in ubuntu

%% save the Xcorr -1 0 +1 for corrleation with EEG
Xcorr10Lag; %12 ses x 12 sorted order x 21 
save('Xcorr10Lag.mat','Xcorr10Lag');

