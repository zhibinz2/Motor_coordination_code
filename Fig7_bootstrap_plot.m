% /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/hyperscanEEG_correlation/zscore_20250515/Copy_of_plots.m
clear
cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/hyperscanEEG_correlation/zscore_20250515

dcsynch=load('dcsynch.mat');

dcsynch4condi=nan(2,7,32);
dcsynch4condi(1,:,:)=dcsynch.follow;
dcsynch4condi(2,:,:)=dcsynch.mutual;

dcsynco=load('dcsynco.mat');
dcsynco4condi=nan(2,7,32);
dcsynco4condi(1,:,:)=dcsynco.follow;
dcsynco4condi(2,:,:)=dcsynco.mutual;

dcsyn2condi4=nan(2,2,7,32);
dcsyn2condi4(1,:,:,:)=dcsynch4condi;
dcsyn2condi4(2,:,:,:)=dcsynco4condi;


bcsynch=load('bcsynch.mat')
bcsynch4condi=nan(2,7,32);
bcsynch4condi(1,:,:)=bcsynch.follow;
bcsynch4condi(2,:,:)=bcsynch.mutual;

bcsynco=load('bcsynco.mat')
bcsynco4condi=nan(2,7,32);
bcsynco4condi(1,:,:)=bcsynco.follow;
bcsynco4condi(2,:,:)=bcsynco.mutual;

bcsyn2condi4=nan(2,2,7,32);
bcsyn2condi4(1,:,:,:)=bcsynch4condi;
bcsyn2condi4(2,:,:,:)=bcsynco4condi;

csyn2condi4=nan(3,2,2,7,32);
csyn2condi4(1,:,:,:,:)=dcsyn2condi4;
csyn2condi4(3,:,:,:,:)=bcsyn2condi4;

syn2names={'synch','synco'};
nx3names={'degree ctr','efficiency','betweenness ctr'};
condi4names={'Independent','Leading','Following','Mutual'};
band7labels = {'\delta', '\theta', '\alpha', ...
               '\mu', '\beta_1', '\beta_2', ...
               '\gamma'};

AllchanNames={'Fp1','Fpz','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6','M1','T7','C3','Cz','C4','T8','M2','CP5','CP1','CP2','CP6','P7','P3','Pz','P4','P8','POz','O1','Oz','O2'};


addpath(genpath('/home/zhibinz2/Documents/GitHub/eeglab'))

cd /home/zhibinz2/Documents/GitHub/Motor_cordination/1_over_f/data_analysis/channels_info
load('chaninfo.mat')


% Extract coordinates
theta = [chaninfo.theta]; % Polar angle in degrees
radius = [chaninfo.radius]; % Radius (0-1)

% Convert to Cartesian coordinates
[x, y] = pol2cart(deg2rad(theta), radius);
electrode_coords = [x', y'];

% Extract original coordinates (example)
% Perform transformation
theta = -90; % degrees
rotation_matrix = [cosd(theta) -sind(theta); 
                 sind(theta)  cosd(theta)];
rotated_coords = (rotation_matrix * electrode_coords')';
final_coords = [rotated_coords(:,1), -rotated_coords(:,2)];

run /home/zhibinz2/Documents/GitHub/Motor_coordination_code/util/color_scheme.m

%% color xylabels 
addpath('/home/zhibinz2/Documents/GitHub/matlab-archive/hnlcode/common/gen_code/color')
% ===== Custom Blue-White-Red Colormap =====
% ===== Custom Blue-White-Red Colormap =====
N = 256; % Number of color steps (even number for symmetry)
% Deep blue [0 0 1] to white [1 1 1]
blue_white = [
    linspace(0, 1, N/2)' ...    % R: 0→1
    linspace(0, 1, N/2)' ...    % G: 0→1
    ones(N/2, 1)                % B: stays 1
];
% White [1 1 1] to deep red [1 0 0]
white_red = [
    ones(N/2, 1) ...           % R: stays 1
    linspace(1, 0, N/2)' ...   % G: 1→0
    linspace(1, 0, N/2)' ...    % B: 1→0
];
% Combine the two halves
custom_map = [blue_white; white_red];

% get collapsed topo values (-2 blue filled, -1 blue circle, 0 mixed black
% circle, +1 red circle, +2 red filled)
topocollapseV=nan(2,2,2,32);% 3 nx x 2 syn type x 2 contrast condition x 32 chan

for ccc=[1 3]
    figure('Position',[53         477        2498         600]);% get(gcf, 'outerposition');figure('Position',[ans]);
    for syn=1:2
        for condi=1:2
            subplot(2,2,2*(condi-1)+syn)
            mat_display = squeeze(csyn2condi4(ccc,syn,condi,:,:));
            imagesc(mat_display);
           
            xticks([1:32]); xticklabels(AllchanNames);
            colormap(custom_map);colorbar;clim([-3 3])
            yticks([1:7]);yticklabels(band7labels);
            ax = gca;
            
            
            % collect the collapsed values
            for col = 1:32
                column_data = mat_display(:, col);
                if any(column_data > 1.96)
                    topocollapseV(ccc,syn,condi,col) = 2;
                elseif all(column_data >= 0 & column_data <= 1.96)
                    topocollapseV(ccc,syn,condi,col) = 1;
                elseif any(column_data < -1.96)
                    topocollapseV(ccc,syn,condi,col) = -2;
                elseif all(column_data < 0 & column_data >= -1.96)
                    topocollapseV(ccc,syn,condi,col) = -1;
                else
                    topocollapseV(ccc,syn,condi,col) = 0;
                end
            end

            % Set up color thresholds
            upper_threshold = 1.96;
            lower_threshold = -1.96;
            
            % Process X-tick labels (columns)
            xTicks = ax.XTick;
            xTickLabels = ax.XTickLabel;
            for idx = 1:length(xTicks)
                % Get current column data
                col_data = mat_display(:, idx);
                
                % Check conditions for column
                has_high = any(col_data > upper_threshold);
                has_low = any(col_data < lower_threshold);
                
                % Determine color
                if has_high && has_low
                    color = darkgreen;
                elseif has_high
                    color = red;
                elseif has_low
                    color = blue;
                else
                    color = 'black'; % default color
                end
                
                % Clear default label
                ax.XTickLabel{idx} = '';
                
                % Create colored text label
                text(xTicks(idx), ax.YLim(2) + 0.5, xTickLabels{idx}, ...
                    'Color', color, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'top', ...
                    'FontSize',9,'FontWeight','bold',...
                    'Parent', ax);
            end
            
            % Process Y-tick labels (rows)
            yTicks = ax.YTick;
            yTickLabels = ax.YTickLabel;
            for idy = 1:length(yTicks)
                if idy <= size(mat_display, 1) % Ensure we don't exceed matrix rows
                    % Get current row data
                    row_data = mat_display(idy, :);
                    
                    % Check conditions for row
                    has_high = any(row_data > upper_threshold);
                    has_low = any(row_data < lower_threshold);
                    
                    % Determine color
                    if has_high && has_low
                        color = darkgreen;
                    elseif has_high
                        color = red;
                    elseif has_low
                        color = blue;
                    else
                        color = 'black'; % default color
                    end
                    
                    % Clear default label
                    ax.YTickLabel{idy} = '';
                    
                    % Create colored text label
                    text(ax.XLim(1) - 0.5, yTicks(idy), yTickLabels{idy}, ...
                        'Color', color, ...
                        'HorizontalAlignment', 'right', ...
                        'VerticalAlignment', 'middle', ...
                        'FontSize',12,'FontWeight','bold',...
                        'Parent', ax);
                    title([syn2names{syn}   '   '  condi4names{condi+2}] )
                end
            end
            
            % Adjust figure margins if needed
            set(ax, 'TickLength', [0 0]); % Remove tick marks
        end
    end
    sgtitle([nx3names{ccc}])
    set(gcf,'color','w'); % set backg
end

%% collapsed topo

for ccc=[1 3]
    figure('Position',[53          41        1000         435]);% get(gcf, 'outerposition');figure('Position',[ans]);
    for syn=1:2
        for condi=1:2
            subplot(2,2,2*(condi-1)+syn)
            array_display = squeeze(topocollapseV(ccc,syn,condi,:));
            % clf
            topoplot([], chaninfo', ...
            'style', 'blank', ...          % No interpolation
            'electrodes', 'off', ...       % Show electrodes
            'emarker', {'.', 'k', 12, 1}, ...  % Black dots, size 12
            'whitebk', 'on', ...          % White background
            'plotrad', 0.5, ...% Head radius
            'hcolor',[0.8 0.8 0.8]);              
            hold on;
            title([syn2names{syn}   '   '  condi4names{condi+2}] )

            for n=1:32
                % set alpha for sig vs non sig roi
                if topocollapseV(ccc,syn,condi,n)==2
                    alpha_tmp=1; dotcolor=red; edgecolor=red; alpha_edge=1;
                elseif topocollapseV(ccc,syn,condi,n)==1
                    alpha_tmp=0; dotcolor=red; edgecolor=red;alpha_edge=0.5;
                elseif topocollapseV(ccc,syn,condi,n)==0
                    alpha_tmp=0; dotcolor=grey; edgecolor=grey; alpha_edge=0.5;
                elseif topocollapseV(ccc,syn,condi,n)==-1
                    alpha_tmp=0; dotcolor=blue; edgecolor=blue;alpha_edge=0.5;
                else topocollapseV(ccc,syn,condi,n)==-2
                    alpha_tmp=1; dotcolor=blue; edgecolor=blue;alpha_edge=1;

                end
            
                scatter(final_coords(n,1), final_coords(n,2), 50, dotcolor, ...
                    'filled', 'MarkerFaceAlpha',alpha_tmp, ...
                    'MarkerEdgeColor',edgecolor,'MarkerEdgeAlpha',alpha_edge,'LineWidth', 2.5);
                axis equal; sizescale=0.58; xlim([-1*sizescale sizescale]); ylim([-1*sizescale sizescale]);
                % drawnow limitrate;
            end
        end
    end
    sgtitle([nx3names{ccc}])
    set(gcf,'color','w'); % set backg
end


