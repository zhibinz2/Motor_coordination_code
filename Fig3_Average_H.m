clear
addpath(genpath('./util'))
run color_scheme.m

data_path= '../Cleaned_data/'; % Cleaned EEG data can be found at https://osf.io/rstpu/. 
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

% Extract the H
H=zeros(2,12,numSes);
for r=1:numSes
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    load([data_path  'clean_' runid '.mat'],'intervals','conditions');
    % sort order
    [x,sortorder]=sort(conditions);
    for j = 1:12
        clear intL intR intL_dmean intR_dmean
        intL = intervals{sortorder(j)}(:,1);
        intR = intervals{sortorder(j)}(:,2);
        % remove the mean
        intL_dmean=intL-mean(intL);
        intR_dmean=intR-mean(intR);
        [~,H(1,j,r)]=DFA_main(intL_dmean);
        [~,H(2,j,r)]=DFA_main(intR_dmean);
    end
end

%% average all trials of the same condition in synch/o (2 subplots)
H_L_uncouple_synch = reshape(H(1,[1:3],[1:2:11]),1,[]);  % L,uncouple,synch
H_L_Llead_synch    = reshape(H(1,[4:6],[1:2:11]),1,[]);  % L,L lead,synch
H_L_Rlead_synch    = reshape(H(1,[7:9],[1:2:11]),1,[]);  % L,R lead,synch
H_L_mutual_synch   = reshape(H(1,[10:12],[1:2:11]),1,[]);% L, mutual,synch
H_R_uncouple_synch = reshape(H(2,[1:3],[1:2:11]),1,[]);  % R,uncouple,synch
H_R_Llead_synch    = reshape(H(2,[4:6],[1:2:11]),1,[]);  % R,L lead,synch
H_R_Rlead_synch    = reshape(H(2,[7:9],[1:2:11]),1,[]);  % R,R lead,synch
H_R_mutual_synch   = reshape(H(2,[10:12],[1:2:11]),1,[]);% R, mutual,synch
H_LR_synch_mean = [mean(H_L_uncouple_synch) mean(H_R_uncouple_synch);...
    mean(H_L_Llead_synch) mean(H_R_Llead_synch);...
    mean(H_L_Rlead_synch) mean(H_R_Rlead_synch);...
    mean(H_L_mutual_synch) mean(H_R_mutual_synch)];

H_L_uncouple_synco = reshape(H(1,[1:3],[2:2:12]),1,[]);  % R,uncouple,synco
H_L_Llead_synco    = reshape(H(1,[4:6],[2:2:12]),1,[]);  % R,L lead,synco
H_L_Rlead_synco    = reshape(H(1,[7:9],[2:2:12]),1,[]);  % R,R lead,synco
H_L_mutual_synco   = reshape(H(1,[10:12],[2:2:12]),1,[]);% R, mutual,synco
H_R_uncouple_synco = reshape(H(2,[1:3],[2:2:12]),1,[]);  % R,uncouple,synco
H_R_Llead_synco    = reshape(H(2,[4:6],[2:2:12]),1,[]);  % R,L lead,synco
H_R_Rlead_synco    = reshape(H(2,[7:9],[2:2:12]),1,[]);  % R,R lead,synco
H_R_mutual_synco   = reshape(H(2,[10:12],[2:2:12]),1,[]);% R, mutual,synco
H_LR_synco_mean = [mean(H_L_uncouple_synco) mean(H_R_uncouple_synco);...
    mean(H_L_Llead_synco) mean(H_R_Llead_synco);...
    mean(H_L_Rlead_synco) mean(H_R_Rlead_synco);...
    mean(H_L_mutual_synco) mean(H_R_mutual_synco)];

%% separate synch and synco with error bar
H_LR_synch_mean=[mean([H_L_uncouple_synch H_R_uncouple_synch]);...
    mean([H_L_Llead_synch H_R_Rlead_synch]);... 
    mean([H_R_Llead_synch H_L_Rlead_synch]);...
    mean([H_L_mutual_synch H_R_mutual_synch])];
H_LR_synco_mean=[mean([H_L_uncouple_synco H_R_uncouple_synco]);...
    mean([H_L_Llead_synco H_R_Rlead_synco]);... 
    mean([H_R_Llead_synco H_L_Rlead_synco]);...
    mean([H_L_mutual_synco H_R_mutual_synco])];
H_LR_synch_std=[std([H_L_uncouple_synch H_R_uncouple_synch]);...
    std([H_L_Llead_synch H_R_Rlead_synch]);... 
    std([H_R_Llead_synch H_L_Rlead_synch]);...
    std([H_L_mutual_synch H_R_mutual_synch])];
H_LR_synco_std=[std([H_L_uncouple_synco H_R_uncouple_synco]);...
    std([H_L_Llead_synco H_R_Rlead_synco]);... 
    std([H_R_Llead_synco H_L_Rlead_synco]);...
    std([H_L_mutual_synco H_R_mutual_synco])];

%% Plot Fig3
canvas(0.23, 0.4);
model_series = [H_LR_synch_mean H_LR_synco_mean];
model_error = [H_LR_synch_std H_LR_synco_std];

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
hold off
xticks(1:4);xticklabels({'Independent','Leader','Follower','Bidirectional'});
xl = get(gca,'XTickLabel');  
set(gca,'XTickLabel',xl,'fontsize',17,'FontWeight','bold')
xlim([0.5 4.5]);
ylabel('DFA Hurst Exponent','FontSize',17); ylim([0.4 1]);
set(gcf,'color','w');
delete(findall(gcf,'type','annotation'))
sg=annotation('textbox',[0.05 0.01 0.5 0.07],'string',...
    {['mean H(matched int) ^{ *PLOT 6}' char(datetime('now'))]})
sg.Rotation=90
yline(0.5,'color', deepyellow,'LineWidth',5)
lg=legend({'Synch','Synco'},'location','north');lg.FontSize-17;
