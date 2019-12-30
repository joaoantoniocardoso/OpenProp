% ----------------------------------------------------------------------- %
%                                                                         %
%                              0111000                                    %
%                           100 1 100 001                                 %
%                         10    1  1  1 00                                %
%                        01  1  1  1      0                               %
%                       0 1  1  1   1  1 1 0                              %
%                       0   1   1   1  1  1 0                             %
%                       0 1     1   1  1  1 0                             %
%                       0 1  1  1   1  0  1 0                             %
%                       0 1  1  1   0  1    0                             %
%                       01 1        1  1 1 0                              %
%                        0    0  1  0 1   0                               %
%                         0         1    0                                %
%                    10010 0 1101111110111                                %
%                  10 1 1  1111111111 11 11                               %
%                 0 1 1 1 11111111101011010111                            %
%                01 11    11111111 1  1    1 110                          %
%               011    1 1 111111110011  1 1 1 110                        %
%               0   11 1 1 1 111      0  1 1 1   10                       %
%               0 1   11  1  0         1 1 1 1 1 1 0                      %
%               1  11 1 1   11          0  1 1 1 1 11                     %
%                0     1 1  0           011  1 1 1 10                     %
%                10 1   1  0             0  1 1 1  11                     %
%                 10     01               01      10                      %
%                   10001                   001 100                       %
%                                             111                         %
%                                                                         %
%             ____                   _____                                %
%            / __ \                 |  __ \                               %
%           | |  | |_ __   ___ _ __ | |__) | __ ___  _ __                 %
%           | |  | | '_ \ / _ \ '_ \|  ___/ '__/ _ \| '_ \                %
%           | |__| | |_) |  __/ | | | |   | | | (_) | |_) |               %
%            \____/| .__/ \___|_| |_|_|   |_|  \___/| .__/                %
%                  | |                              | |                   %
%                  |_|                              |_|                   %
%                                                                         %
%             An integrated rotor design and analysis tool.               %
%                                                                         %
%                                                                         %
% Copyright (C) 2011, Brenden Epps.                                       %
%                                                                         %
% This program is free software; you can redistribute it and/or modify it %
% under the terms of the GNU General Public License as published by the   %
% Free Software Foundation.                                               %
%                                                                         %
% This program is distributed in the hope that it will be useful, but     %
% WITHOUT ANY WARRANTY; without even the implied warranty of              %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    %
% See the GNU General Public License for more details.                    %
%                                                                         %
% ----------------------------------------------------------------------- %

function [] = Make_2D_Blade_Image(RG,x2Dr,y2Dr, GUI_flag,Plots,PlotPanels)

if nargin == 3
    GUI_flag = 0;
end

Mp = size(x2Dr,1) - 1;
Np = size(x2Dr,2)/2;


if GUI_flag 
            
    set(0,'CurrentFigure',Plots);
    h = axes('parent',PlotPanels(13));
        axes(h);

else    
    Fig2_S = figure('units','normalized','position',[0.31 .06 .4 .3],'name','Blade Image','numbertitle','off');
end
            hold on;
            axis equal;     
            grid on;
            box on,
            title('2D Blade Image');  
            xlabel('X (2D) [m]');  
            ylabel('Y (2D) [m]');

            style      = ['r' 'g' 'b' 'm' 'k'];
            str_prefix = {'r/R = '};            

flag = 1;

plot(  [min(x2Dr),max(x2Dr)],0*[min(y2Dr),max(y2Dr)],'k')
plot(0*[min(x2Dr),max(x2Dr)],  [min(y2Dr),max(y2Dr)],'k')

for i = 1:ceil(Mp/5):Mp     % for five radial sections from root to tip
    handle_legend(flag) = plot(x2Dr(i,:),y2Dr(i,:),style(flag),'linewidth',2);
    
    
%     plot(x2Dr(i,1:Np),y2Dr(i,1:Np),'g','linewidth',2)
%     plot(x2Dr(i,Np+1:Np+Np),y2Dr(i,Np+1:Np+Np),'r','linewidth',2)


    plot([0.5*(x2Dr(i,1)+x2Dr(i,2*Np)),0.5*(x2Dr(i,Np)+x2Dr(i,Np+1))],[0.5*(y2Dr(i,1)+y2Dr(i,2*Np)),0.5*(y2Dr(i,Np)+y2Dr(i,Np+1))],style(flag),'linewidth',1);
    
    str_legend(flag) = strcat(str_prefix,num2str(RG(i)));
    
    flag = flag+1;
end

legend(handle_legend,str_legend,'location','northwest');
