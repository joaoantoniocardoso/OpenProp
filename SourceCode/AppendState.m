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
% MIT OpenProp_v3.2.0
% Last modified: 11/15/2011 Brenden Epps
% =========================================================================
%
% This function appends a state onto states data structure "s"
%
%--------------------------------------------------------------------------


function s = AppendState(s,L,G,UASTAR,UTSTAR,VSTAR,TANBC,TANBIC,ALPHA,CL,CD, KT,KQ,CT,CQ,CP,EFFY,VMWV,  Duct_flag,Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT,TAU,CTD)

if nargin == 17
    Duct_flag = 0;
end

% s.part1 '------ Off-design states ------';
s.L       = [s.L;     L];                       % [NL x 1] tip speed ratio 
s.Js      = [s.Js; pi/L];                       % [NL x 1] advance coefficient 

s.G       = [s.G;           G'];                % [NL x Mp] circulation distribution
s.UASTAR  = [s.UASTAR; UASTAR];                 % [NL x Mp] axial      induced velocity distribution
s.UTSTAR  = [s.UTSTAR; UTSTAR];                 % [NL x Mp] tangential induced velocity distribution

s.VSTAR   = [s.VSTAR;  VSTAR];                  % [NL x Mp] total inflow       velocity distribution
s.TANBC   = [s.TANBC;  TANBC];                  % [NL x Mp] tangent of free-stream inflow angle
s.TANBIC  = [s.TANBIC; TANBIC];                 % [NL x Mp] tangent of total inflow angle
s.ALPHA   = [s.ALPHA;  ALPHA];                  % [NL x Mp] ALPHA = alpha - alphaI, RADIANS
% % s.deltaALPHA1 = [];                         % [NL x Mp] 
% % s.deltaALPHA2 = [];                         % [NL x Mp]
s.CL      = [s.CL; CL];                         % [NL x Mp] lift coefficient distribution
s.CD      = [s.CD; CD];                         % [NL x Mp] drag coefficient distribution

% s.part3 '------ Duct parameters ------';                         
if Duct_flag == 1  % Duct is present 
    s.Gd      = [s.Gd; Gd];                     % [NL x 1]  duct circulation
    s.UARINGq = [s.UARINGq; UARINGq];           % [NL x 1]  UA    at duct quarter chord
    s.URRINGq = [s.URRINGq; URRINGq];           % [NL x 1]  UR    at duct quarter chord 
    s.VSRINGq = [s.VSRINGq; VSRINGq];           % [NL x 1]  VSTAR at duct quarter chord
    s.BetaIDq = [s.BetaIDq; BetaIDq];           % [NL x 1]  angle at duct quarter chord 
    s.CLd     = [s.CLd;         CLd];           % [NL x 1]  duct lift coefficient
    
    s.UARING  = [s.UARING;   UARING];           % [NL x Nd], induced velocity at duct
    s.URRING  = [s.URRING;   URRING];           % [NL x Nd], induced velocity at duct
    
    s.UADUCT  = [s.UADUCT;   UADUCT];           % [NL x Mp], induced velocity at propeller
    
    s.TAU     = [s.TAU; TAU];                   % [NL x 1]  thrust ratio 
    s.CTD     = [s.CTD; CTD];                   % [NL x 1]  CT for duct
end

% s.part4 '------ Performance metrics ------';
s.KT      = [s.KT;   KT];     % [NL x 1]  thrust coefficient
s.KQ      = [s.KQ;   KQ];     % [NL x 1]  torque coefficient
s.CT      = [s.CT;   CT];     % [NL x 1]  thrust coefficient
s.CQ      = [s.CQ;   CQ];     % [NL x 1]  torque coefficient
s.CP      = [s.CP;   CP];     % [NL x 1]  power  coefficient
s.EFFY    = [s.EFFY; EFFY];   % [NL x 1]  efficiency
s.VMWV    = [s.VMWV; VMWV];   % [NL x 1]  volumetric mean wake velocity
end

