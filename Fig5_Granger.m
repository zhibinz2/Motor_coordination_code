clear
addpath(genpath('./util'))
run color_scheme.m

% For d estimate and removal
% http://www.lucafaes.net/LMSE-MSE_VARFI.html
% addpath(genpath('./util/MSE-VARFI'))

% Granger Causality
addpath(genpath('./util/MVGC1'));
run ./util/MVGC1/startup.m

data_path= '../Cleaned_data/'; % Cleaned EEG data can be found at https://osf.io/rstpu/. 
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

% d estimation from H (DFA method) based on the good matched intervals
H=zeros(2,12,numSes);
for r=1:numSes
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    load([data_path  'clean_' runid '.mat'],'intervals','conditions');
    % sort order and plot results
    [x,sortorder]=sort(conditions);
    for j = 1:12
        clear intL_good_dmean intR_good_dmean
        % remove the mean
        intL_good_dmean=intervals{sortorder(j)}(:,1)-mean(intervals{sortorder(j)}(:,1));
        intR_good_dmean=intervals{sortorder(j)}(:,2)-mean(intervals{sortorder(j)}(:,2));
        % estimate the d
        [~,H(1,j,r)]=DFA_main(intL_good_dmean);
        [~,H(2,j,r)]=DFA_main(intR_good_dmean);
    end
end
d1=H-0.5;

% d removal
int_dmean_drm=cell(2,12,numSes);
for r=1:numSes
    clear intervals conditions sortorders
    runid = num2str(seeds(r,:));
    load([data_path  'clean_' runid '.mat'],'intervals','conditions');
    % sort order and plot results
    [x,sortorder]=sort(conditions);
    for j = 1:12
        clear intL_good_dmean intR_good_dmean
        % remove the mean
        intL_good_dmean=intervals{sortorder(j)}(:,1)-mean(intervals{sortorder(j)}(:,1));
        intR_good_dmean=intervals{sortorder(j)}(:,2)-mean(intervals{sortorder(j)}(:,2));
        % d removal  (method1: based on DFA, d=h-0.5)
        [int_dmean_drm{1,j,r}]=remove_d(intL_good_dmean,d1(1,j,r));
        [int_dmean_drm{2,j,r}]=remove_d(intR_good_dmean,d1(2,j,r));
    end
end

% orgainized the synchronization trials into 4 conditions
condi4Ind={[1:3],[4:6],[7:9],[10:12]};
syn2Ind={[1:2:11],[2:2:12]};syn2names={'Synch','Synco'};

%% MVGC
Fs=nan(2,4);
for condi=1
    for t=1:2
    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{condi},syn2Ind{t});
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),18,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),18,1));
    X=[];
    X=[L_mat R_mat; R_mat L_mat];
    X=permute(X,[2,1]);
    [F,~,~] = myGCautocov(X);
    Fs(t,1)=F(2,1);
    end
end

for condi=[2 3]
    for t=1:2
    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{2},syn2Ind{t});
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),18,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),18,1));
    X1=[];
    X1=[L_mat R_mat];

    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{3},syn2Ind{t});
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),18,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),18,1));
    X2=[];
    X2=[R_mat L_mat];
    
    X=[];
    X=[X1; X2];
    X=permute(X,[2,1]);
    [F,~,~] = myGCautocov(X);
    Fs(t,2)=F(2,1);Fs(t,3)=F(1,2);
    end
end

for condi=4
    for t=1:2
    test_data=[];
    test_data=int_dmean_drm(:,condi4Ind{condi},syn2Ind{t});
    L_mat=[];R_mat=[];
    L_mat=cell2mat(reshape(squeeze(test_data(1,:,:)),18,1));
    R_mat=cell2mat(reshape(squeeze(test_data(2,:,:)),18,1));
    X=[];
    X=[L_mat R_mat; R_mat L_mat];
    X=permute(X,[2,1]);
    [F,~,~] = myGCautocov(X);
    Fs(t,4)=F(2,1);
    end
end


%% Plot Fig5
canvas(0.4, 0.5);
model_series = Fs';
b = bar(model_series, 'FaceColor','flat');
b(1).FaceColor=darkgreen;b(2).FaceColor=pink;
xticks(1:4);xticklabels({'Uncoupled','Leader -> Follower','Follower -> Leader','Mutual'});
xl = get(gca,'XTickLabel');  
set(gca,'XTickLabel',xl,'fontsize',20,'FontWeight','bold')
legend({'Synch','Synco'},'location','northwest','FontSize', 30);
ylabel('GC','FontSize',17); 
yl = get(gca,'YTickLabel');  
set(gca,'YTickLabel',yl,'fontsize',20,'FontWeight','bold');
set(gcf,'color','w');
