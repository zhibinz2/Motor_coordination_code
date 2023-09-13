clear
addpath(genpath('./util'))
run color_scheme.m

data_path= '../Cleaned_data/'; % Cleaned EEG data can be found at https://osf.io/rstpu/. 
seeds=[20220713;20220721;20220804;20220808;20220810;20220811;20220815;20220816;20221003;2022100401;
        2022100402;20221005];
numSes=size(seeds,1);

% Compute and extract H
H=zeros(2,12,numSes);
for s=1:numSes;
    clear intervals conditions
    runid=num2str(seeds(s,:));
    load([data_path 'clean_' runid '.mat'],'intervals','conditions')  
    [x,sortorder]=sort(conditions);
    for t = 1:12
        clear intL intR intL_dmean intR_dmean
        intL = intervals{sortorder(t)}(:,1);
        intR = intervals{sortorder(t)}(:,2);
        % remove the mean
        intL_dmean=intL-mean(intL);
        intR_dmean=intR-mean(intR);
        [~,H(1,t,s)]=DFA_main(intL_dmean);
        [~,H(2,t,s)]=DFA_main(intR_dmean);
    end
end

H_L_uncouple_synch = reshape(H(1,[1:3],[1:2:11]),1,[])  ;% L,uncouple,sync
H_L_Llead_synch    = reshape(H(1,[4:6],[1:2:11]),1,[])  ;% L,L lead,synch
H_L_Rlead_synch    = reshape(H(1,[7:9],[1:2:11]),1,[])  ;% L,R lead,synch
H_L_mutual_synch   = reshape(H(1,[10:12],[1:2:11]),1,[]);% L, mutual,synch
H_R_uncouple_synch = reshape(H(2,[1:3],[1:2:11]),1,[])  ;% R,uncouple,synch
H_R_Llead_synch    = reshape(H(2,[4:6],[1:2:11]),1,[])  ;% R,L lead,synch
H_R_Rlead_synch    = reshape(H(2,[7:9],[1:2:11]),1,[])  ;% R,R lead,synch
H_R_mutual_synch   = reshape(H(2,[10:12],[1:2:11]),1,[]);% R, mutual,synch

H_L_uncouple_synco = reshape(H(1,[1:3],[2:2:12]),1,[])  ;% R,uncouple,synco
H_L_Llead_synco    = reshape(H(1,[4:6],[2:2:12]),1,[])  ;% R,L lead,synco
H_L_Rlead_synco    = reshape(H(1,[7:9],[2:2:12]),1,[])  ;% R,R lead,synco
H_L_mutual_synco   = reshape(H(1,[10:12],[2:2:12]),1,[]);% R, mutual,synco
H_R_uncouple_synco = reshape(H(2,[1:3],[2:2:12]),1,[])  ;% R,uncouple,synco
H_R_Llead_synco    = reshape(H(2,[4:6],[2:2:12]),1,[])  ;% R,L lead,synco
H_R_Rlead_synco    = reshape(H(2,[7:9],[2:2:12]),1,[])  ;% R,R lead,synco
H_R_mutual_synco   = reshape(H(2,[10:12],[2:2:12]),1,[]);% R, mutual,synco
%% Two-way ANOVA
[p,~,stats] = anova2([...
    [[H_L_uncouple_synch H_R_uncouple_synch]' [H_L_uncouple_synco H_R_uncouple_synco]'];...
    [[H_L_Llead_synch H_R_Rlead_synch]' [H_L_Llead_synco H_R_Rlead_synco]'];... 
    [[H_R_Llead_synch H_L_Rlead_synch]' [H_R_Llead_synco H_L_Rlead_synco]'];...
    [[H_L_mutual_synch H_R_mutual_synch]' [H_L_mutual_synco H_R_mutual_synco]']...
    ],36); 
