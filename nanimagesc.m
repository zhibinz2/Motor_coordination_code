function nanimagesc(mat,nanval,cmp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
data = mat;
data(isnan(data)) = nanval; % replace NaN values with -0.01 for grey

% Display the matrix using imagesc
% figure;clf
imagesc(data);

% Create a custom colormap
jetMap = cmp;
grayColor = [1, 1, 1]; % Light Gray color for NaN values
customMap = [grayColor; jetMap]; % Add gray as the first color
% Apply the custom colormap
colormap(customMap);

% Adjust color axis to map -1 to gray
if min(mat(:)) == max(data(:))
    % do nothing
else
    clim([min(mat(:)), max(mat(:))]);
end

% % Add a colorbar
colorbar
% colorbarHandle = colorbar; % Create a colorbar
% % Get current tick labels
% currentTicks = colorbarHandle.Ticks;
% % Replace the tick label corresponding to 0 with 'NaN'
% tickLabels = string(currentTicks); % Convert numeric ticks to strings

% if min(mat(:))<0
%     % do nothing
% elseif min(mat(:))==0
%     % do nothing
% else
%     tickLabels(currentTicks == 0) = "NaN"; % Replace 0 with 'NaN'
% end

% Apply the modified tick labels to the colorbar
% colorbarHandle.TickLabels = tickLabels;
end