syn2names={'synch','synco'};
nx3names={'degree ctr','efficiency','betweenness ctr'};
condi4names={'Independent','Leading','Following','Mutual'};
band7labels = {'\delta', '\theta', '\alpha', ...
               '\mu', '\beta_1', '\beta_2', ...
               '\gamma'};

AllchanNames={'Fp1','Fpz','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6','M1','T7','C3','Cz','C4','T8','M2','CP5','CP1','CP2','CP6','P7','P3','Pz','P4','P8','POz','O1','Oz','O2'};



%% correlate the betweenness centrality between A and B in hyperscanning at each electrode
% open /home/zhibinz2/Documents/GitHub/Motor_coordination_code/Fig7_bootstrap_plot.m
csyn2condi4; % 3nx x 2 syn x 2 condi (follow and mutual) x 7 freq x 32 chan
cd /home/zhibinz2/Documents/GitHub/Motor_coordination_code
load('nx144.mat'); % 144 trials  x 2 subjects x 7 freq x 32 chan
% dg_ctr144;
% efficiency144;
% bw_ctr144;

load('cc3_syn.mat'); %3nx x 2 syn x 4 condi x36 tr x 7 freq x 32 chan
btwAB=squeeze(cc3_syn(3,:,3:4,:,:,:)); % btw nx x 2syn x 2 condi(follow, mutual) x 36 tr x 7 freq x 32 chan
corr_btwAB=nan(2,2,7,32);
for syn=1:2
    for condi=1:2
        for freq=1:7
            for chan=1:32
                A=squeeze(btwAB(syn,condi,1:18,freq,chan));
                B=squeeze(btwAB(syn,condi,19:36,freq,chan));
                corr_btwAB(syn,condi,freq,chan)=corr(A,B);
            end
        end
    end
end

figure;clf
for syn=1:2
    for condi=1:2
        subplot(2,2,2*(syn-1)+condi);
        mat=squeeze(corr_btwAB(syn,condi,:,:));
        imagesc(mat);colorbar;colormap('jet');clim([-0.5 0.5])
        title([num2str(syn) '' num2str(condi)])
    end
end
set(gcf,'color','w'); % set backg
%% Correlate betweenness centrality to complexity score
load('H_syn.mat');
H_syn;% 2 syn x 4 condi x 36 tr
corr_btw_H=nan(2,2,7,32);
for syn=1:2
    for condi=1:2
        H_tmp=H_syn{syn}{condi+2}';
        for freq=1:7
            for chan=1:32
                btw_tmp=squeeze(btwAB(syn,condi,:,freq,chan));
                corr_btw_H(syn,condi,freq,chan)=corr(H_tmp,btw_tmp);
            end
        end
    end
end

figure;clf
for syn=1:2
    for condi=1:2
        subplot(2,2,2*(syn-1)+condi);
        mat=squeeze(corr_btw_H(syn,condi,:,:));
        imagesc(mat);colorbar;colormap('jet');clim([-0.5 0.5])
        title([num2str(syn) '' num2str(condi)])
    end
end
set(gcf,'color','w'); % set