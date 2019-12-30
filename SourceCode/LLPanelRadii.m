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

% =========================================================================
% =================================================== LLPanelRadii Function
% Last Modified: 10/21/2011, Brenden Epps
%
% Lifting Line panel radii layout
%
% RC == radius of control points / propeller radius
% RV == radius of vortex  points / propeller radius
% DR == difference in vortex radii / propeller radius

function [RC,RV,DR] = LLPanelRadii(Mp,Rhub_oR,Hub_flag,Duct_flag)

if     Duct_flag == 0 && Hub_flag == 0

   % Constant spacing -- 1/4 panel inset at hub and tip  (use with no hub or duct image)
    RV = Rhub_oR + (1-Rhub_oR) * ((0:Mp)+0.25)/(Mp+0.50);
    RC = Rhub_oR + (1-Rhub_oR) * ((1:Mp)-0.25)/(Mp+0.50);
    
elseif Duct_flag == 1 && Hub_flag == 0 
    
   % Constant spacing -- 1/4 panel inset at hub only  (use with no hub image but yes duct image)
    RV = Rhub_oR + (1-Rhub_oR) * ((0:Mp)+0.25)/(Mp+0.25);
    RC = Rhub_oR + (1-Rhub_oR) * ((1:Mp)-0.25)/(Mp+0.25);
    
elseif Duct_flag == 1 && Hub_flag == 1    

    % Constant spacing -- no inset  (use with both duct image and hub image)
    RV = Rhub_oR + (1-Rhub_oR) * ((0:Mp)     )/Mp;
    RC = Rhub_oR + (1-Rhub_oR) * ((1:Mp)-0.50)/Mp;

elseif Duct_flag == 0 && Hub_flag == 1
    
    % Constant spacing -- 1/4 panel inset at tip only  (use with hub image but no duct image)
    RV = Rhub_oR + (1-Rhub_oR) * ((0:Mp)     )/(Mp+0.25);
    RC = Rhub_oR + (1-Rhub_oR) * ((1:Mp)-0.50)/(Mp+0.25);

end

DR = diff(RV);   

