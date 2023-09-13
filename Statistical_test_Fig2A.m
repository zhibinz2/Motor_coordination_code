clear
addpath(genpath('./util'))
run color_scheme.m

data_path= '../Cleaned_data/'; % Cleaned EEG data can be found at https://osf.io/rstpu/. 
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

% load tapping intervals
intervals_all=nan(2,12,12); % 2 subjects x 12 trials x 12 session
for s=1:numSes
    runid=num2str(seeds(s,:));
    clear intervals
    load([data_path 'clean_' runid '.mat'],'intervals','conditions');
    [x,sortorder]=sort(conditions);
    for t=1:12
        clear intL intR 
        intL = intervals{sortorder(t)}(:,1);
        intR = intervals{sortorder(t)}(:,2);
        intervals_all(1,t,s)=mean(intL)/2; % L subject 
        intervals_all(2,t,s)=mean(intR)/2; % R subject 
    end
end

intervals_L_uncouple_synch = reshape(intervals_all(1,[1:3],[1:2:11]),1,[])  ;% L,uncouple,synch
intervals_L_Llead_synch    = reshape(intervals_all(1,[4:6],[1:2:11]),1,[])  ;% L,L lead,synch
intervals_L_Rlead_synch    = reshape(intervals_all(1,[7:9],[1:2:11]),1,[])  ;% L,R lead,synch
intervals_L_mutual_synch   = reshape(intervals_all(1,[10:12],[1:2:11]),1,[]);% L, mutual,synch
intervals_R_uncouple_synch = reshape(intervals_all(2,[1:3],[1:2:11]),1,[])  ;% R,uncouple,synch
intervals_R_Llead_synch    = reshape(intervals_all(2,[4:6],[1:2:11]),1,[])  ;% R,L lead,synch
intervals_R_Rlead_synch    = reshape(intervals_all(2,[7:9],[1:2:11]),1,[])  ;% R,R lead,synch
intervals_R_mutual_synch   = reshape(intervals_all(2,[10:12],[1:2:11]),1,[]);% R, mutual,synch

intervals_L_uncouple_synco = reshape(intervals_all(1,[1:3],[2:2:12]),1,[])  ;% R,uncouple,synco
intervals_L_Llead_synco    = reshape(intervals_all(1,[4:6],[2:2:12]),1,[])  ;% R,L lead,synco
intervals_L_Rlead_synco    = reshape(intervals_all(1,[7:9],[2:2:12]),1,[])  ;% R,R lead,synco
intervals_L_mutual_synco   = reshape(intervals_all(1,[10:12],[2:2:12]),1,[]);% R, mutual,synco
intervals_R_uncouple_synco = reshape(intervals_all(2,[1:3],[2:2:12]),1,[])  ;% R,uncouple,synco
intervals_R_Llead_synco    = reshape(intervals_all(2,[4:6],[2:2:12]),1,[])  ;% R,L lead,synco
intervals_R_Rlead_synco    = reshape(intervals_all(2,[7:9],[2:2:12]),1,[])  ;% R,R lead,synco
intervals_R_mutual_synco   = reshape(intervals_all(2,[10:12],[2:2:12]),1,[]);% R, mutual,synco

%% Two-way ANOVA
[p,tbl,stats] = anova2([...
    [[intervals_L_uncouple_synch intervals_R_uncouple_synch]' [intervals_L_uncouple_synco intervals_R_uncouple_synco]'];...
    [[intervals_L_Llead_synch intervals_R_Rlead_synch]' [intervals_L_Llead_synco intervals_R_Rlead_synco]'];... 
    [[intervals_R_Llead_synch intervals_L_Rlead_synch]' [intervals_R_Llead_synco intervals_L_Rlead_synco]'];...
    [[intervals_L_mutual_synch intervals_R_mutual_synch]' [intervals_L_mutual_synco intervals_R_mutual_synco]']...
    ],36); 
