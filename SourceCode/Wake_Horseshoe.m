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
% ================================================= Wake Horseshoe Function
% Last update: 11/2/2011 Brenden Epps
%
% -------------------------------------------------------------------------
% This function computes the vortex Horseshoe Influence Functions given
% in Kerwin, p.179.
%
% UAHIF = 2*pi*R*(HIF in Kerwin)
% UTHIF = 2*pi*R*(HIF in Kerwin)
%
% UAHIF(n,m) = influence of mth horseshoe vortex on nth control point
%
% -------------------------------------------------------------------------

function [UAHIF,UTHIF,URHIF] = Wake_Horseshoe(Mp,Z,TANBIV,RC,RV,SCF,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR,WX,WY,WZ,epsilon)
% 
% RC      = 0.25:0.1:0.95;
% Mp      = 8;
% SCF     = 1;
% epsilon = 0.01;

% function [UAHIF,UTHIF,URHIF] = Wake_Horseshoe(Mp,Z,RC,SCF,WX,WY,WZ,epsilon)

                                
UAHIF = zeros(Mp,Mp);
UTHIF = zeros(Mp,Mp);
URHIF = zeros(Mp,Mp);
UAW   = zeros( 1,Mp+1);
UTW   = zeros( 1,Mp+1);
URW   = zeros( 1,Mp+1);

for n = 1:Mp                 % for each control point, n     (FOR LOOP MF2)
    for m = 1:Mp+1           % for each vortex  point, m     (FOR LOOP MF3)

        % --- Find velocity induced at RC(n) by a unit vortex shed at RV(m)
        % %     (Wrench returns 2*pi*R*u_bar)
        % [UAW(m),UTW(m)] = Wrench(Z,TANBIV(m),RC(n),RV(m));

        XH  = WX(m,:); % X-position of segment endpoints of trailing vortex
        YH  = WY(m,:); % Y-position of segment endpoints of trailing vortex
        ZH  = WZ(m,:); % Z-position of segment endpoints of trailing vortex
        Xpf = [0, 0, RC(n)];  % [X,Y,Z] position of control point
        
        [UAW(m), UTW(m), URW(m)] = Wake_Influence(XH,YH,ZH,Xpf,epsilon,Z);
                   

        % ---------------------------- Find hub-image effects, Kerwin p.181
        if Hub_flag == 1
            RVW    = Rhub_oR^2/RV(m);               % Kerwin eqn 259, p.183
            TANBIW = TANBIV(1)*RV(1)/RVW;

            [UAWh,UTWh] = Wrench(Z,TANBIW,RC(n),RVW); 
            
            UAW(m) = UAW(m) - UAWh; 
            UTW(m) = UTW(m) - UTWh;

        end

        % ----------------------------------------- Find duct-image effects
        if Duct_flag == 1
            % RVD = radius of duct-image trailer (Coney p.69, eqn 3.21)
            % (mirrored the method used above for the hub image)
            RVD    = Rduct_oR^2/RV(m);
            TANBID = TANBIV(end)*RV(end)/RVD;
                
            [UAWd,UTWd] = Wrench(Z,TANBID,RC(n),RVD);

            UAW(m) = UAW(m) - UAWd;
            UTW(m) = UTW(m) - UTWd;
        end
    end                                                % (END FOR LOOP MF3)
    
    % ---------------------------- Apply the swirl cancellation factor, SCF
    UTW = UTW*SCF;              % SCF = Swirl Cancellation Factor
    

    % -------------------------- Determine the Horseshoe Influence Function
    % The Horseshoe Influence Function for vortex panel m is the
    % effect of the induction by a helical trailing vortex at
    % vortex point m   with circulation -Gamma(m) and another at
    % vortex point m+1 with circulation +Gamma(m).
    % UAHIF(n,m) = u_barA horseshoe influence function in eqn 254.
    % UAW(m)     = u_barA Wrench velocity given in eqn 202-203.
    for m = 1:Mp                                % for each vortex  panel, m
        UAHIF(n,m) = UAW(m+1)-UAW(m);           % 2*pi*R*(HIF)
        UTHIF(n,m) = UTW(m+1)-UTW(m);           % 2*pi*R*(HIF)
        URHIF(n,m) = URW(m+1)-URW(m);           % 2*pi*R*(HIF)
    end

end                                                    % (END FOR LOOP MF2)
% ================================================== END Horseshoe Function
% =========================================================================