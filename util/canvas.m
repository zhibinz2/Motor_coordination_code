function [FigureXpixels,FigureYpixels] = canvas(xr,yr);
% This function draw a canvas with specified size for the figure
% based on your combined screen size
% xr yr are ratio to the size of the screen
ScreenSize=get(0,'MonitorPositions');
FigureXpixels=ScreenSize(3);FigureYpixels=ScreenSize(4);
figure('units','pixels','position',[0 0 FigureXpixels*xr FigureYpixels*yr]);
end