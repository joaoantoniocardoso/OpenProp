% ------------------------------------------------------  plotaxes(H,style)
% 
% This function plots the axes on the figure with axis handle H.
%
% Example:
%           plotaxes(gca,'r')
%
% -------------------------------------------------------------------------


function [ ] = plotaxes(H,style)

if nargin == 0 
    H = gca;
    style = 'k-';
elseif nargin == 1 
    style = 'k-';
end



D = get(H);

plot([0 0] ,D.YLim,style)
plot(D.XLim,[0 0] ,style)

