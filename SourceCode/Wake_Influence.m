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
% ======================================================== Wake_Influence.m
% Last update: 11/2/2011 Brenden Epps
%
% -------------------------------------------------------------------------
% This function evaluates the influence of a unit strength vortex line.
%
% [XH,YH ZH] = equal length vectors of vortex segment endpoint positions
% Xpf        = [x,y,z] position of field point at which to evaluate velocities
% epsilon    = vortex core size
% Z          = number of propeller blades
% -------------------------------------------------------------------------

function [UA UT UR] = Wake_Influence(XH,YH,ZH,Xpf,epsilon,Z)

% function [VELhelix] = Helix_Influence(XH,YH,ZH,Xpf,epsilon,Z)

%%
% XH = WX(end,:);
% YH = WX(end,:);
% ZH = WX(end,:);
% Xpf = [0 0 0.21];
% epsilon = 0.001;

% Check that XH,YH,ZH are all same size and are vectors
if length(XH) ~= length(YH) | length(XH) ~= length(ZH) | length(Xpf) ~= 3
    display('WARNING: Helix_Influence input data is not formatted correctly.')
    return
end

if size(XH,2) == 1
    XH = XH';
    YH = YH';
    ZH = ZH';
end

if size(Xpf,2) == 1
    Xpf = Xpf';
end


S = length(XH) - 1;  % number of segments in each helix


HPData = [XH;YH;ZH]; % Helix position data, 
% HPData is formatted such that:  XH(s) == HPData(1,s), and a helix can be
% plotted using:
%              plot3(HPData(1,:),HPData(2,:),HPData(3,:),'b','linewidth',2) 


theta_Z = 0:360/Z:360*(Z-1)/Z;            % angle between blades [deg]
        
UA = 0;
UT = 0;
UR = 0;        
        
for k = 1:Z % for each blade
    
    Rotation = [1,                0,                 0; ...
                0, cosd(theta_Z(k)), -sind(theta_Z(k)); ...
                0, sind(theta_Z(k)),  cosd(theta_Z(k))];
    
            
    HPDr = Rotation*HPData;  % rotated helix position data
    
    % plot3(HPDr(1,:),HPDr(2,:),HPDr(3,:),'r','linewidth',2)
    
    for s = 1:S  % for each vortex segment
        % Helix point (1) is on the lifting line and point (S+1) is far
        % downstream.  Therefore, computing the effect of a helix with the
        % direction of the circulation pointing AWAY FROM the lifting line
        % amounts to evaluating segments for wich each segment starts at 
        % pont (s) and ends at point (s+1).   
                                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % [VELseg] = VortexSegment(Xp1,Xp2,Xpf,epsilon)
        % Xp1 = (x1,y1,z1) = 3-element position vector of start of segment 
        % Xp2 = (x2,y2,z2) = 3-element position vector of end   of segment
        % Xpf = (xf,yf,zf) = 3-element position vector of field point     
        % epsilon = vortex core radius         
        %
        % % The 2*pi coefficient is part of the non-dimensionalization.
        % VELseg = 2*pi*VortexSegment([HPDr(1,s+1),HPDr(2,s+1),HPDr(3,s+1)], ...
        %                             [HPDr(1,s)  ,HPDr(2,s),  HPDr(3,s)],Xpf,epsilon);
        

        Xp1 = [HPDr(1,s)  ,HPDr(2,s),  HPDr(3,s)];
        Xp2 = [HPDr(1,s+1),HPDr(2,s+1),HPDr(3,s+1)];

        X0 = Xp2 - Xp1; % vortex from point p1 to point p2
        X1 = Xpf - Xp1; % vortex from point p1 to point pf
        X2 = Xpf - Xp2; % vortex from point p2 to point pf

        % ------ Biot-Savart Law for straight line segment with vortex core
               crossX1X2 = [X1(2)*X2(3)-X1(3)*X2(2) , ...
                            X1(3)*X2(1)-X1(1)*X2(3) , ...
                            X1(1)*X2(2)-X1(2)*X2(1)];

        SQnorm_crossX1X2 = crossX1X2(1)^2 + crossX1X2(2)^2 + crossX1X2(3)^2;
        SQnorm_X0        =        X0(1)^2 +        X0(2)^2 +        X0(3)^2;
        SQnorm_X1        =        X1(1)^2 +        X1(2)^2 +        X1(3)^2;
        SQnorm_X2        =        X2(1)^2 +        X2(2)^2 +        X2(3)^2;

        dot_X0X1         =  X0(1)*X1(1)   +  X0(2)*X1(2)   +  X0(3)*X1(3);
        dot_X0X2         =  X0(1)*X2(1)   +  X0(2)*X2(2)   +  X0(3)*X2(3);

        epsilonSQ        = epsilon^2;

        % Biot-Savart Law.  Note: the leading factor of 0.5 is 2*pi times
        % the front factor in VortexSegment(...), which accounts for the
        % non-dimensaionalization of UAHIF, UTHIF, URHIF
        VELseg = 0.5                                                  ...
                *(1 / (SQnorm_crossX1X2 + epsilonSQ * SQnorm_X0) )    ...
                *( dot_X0X1 / sqrt(SQnorm_X1+epsilonSQ)               ...
                  -dot_X0X2 / sqrt(SQnorm_X2+epsilonSQ) )             ...
                *crossX1X2;      
        % -----------------------------------------------------------------                               
         
        UA = UA - VELseg(1);  % UASTAR positive in  -X directon (-axial      directon)
        UT = UT + VELseg(2);  % UTSTAR positive in  +Y directon (-tangential directon)
        UR = UR + VELseg(3);  % URSTAR positive in  +Z directon (+radial     directon)        
    end
end