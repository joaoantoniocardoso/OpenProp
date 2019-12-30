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
% ========================================================= Forces Function
% Last Modified: 11/2/2011, Brenden Epps
%
% This function computes the thrust, torque, and power coefficients, and it
% computes the efficiency of the propeller.
%
% -------------------------------------------------------------------------
% Input Variables:
%       RC                radius of control point / R
%       DR                vortex panel width / R
%       VAC,VTC           axial/tangential inflow  velocity / Vs
%       UASTAR,UTSTAR    	axial/tangential induced velocity / Vs
%       UADUCT            axial velocity induced by duct    / Vs
%       CD                section drag coefficient
%       CoD               section chord length / D
%       G                 circulation / (2*pi*R*Vs)
%
%       Z                 number of blades
%       Js == Vs/(n*D)    advance coefficient 
%       VMIV              volumetric mean inflow velocity / Vs
%       Hub_flag          hub image flag (1 = yes, 0 = no)
%       Rhv               hub vortex radius / hub radius
%       CTD               CT for the duct (CTD == 0 for no duct)
%
% Define:
%       Va == VMIV*Vs     volumetric mean inflow velocity (dimensional)
%       Vw == VMWV*Vs     volumetric mean wake   velocity (dimensional)
%
%
% Output variables:
%
%       Ja == Va/(n*D)    inflow-adapted advance coefficient
%       Jw == Vw/(n*D)      wake-adapted advance coefficient
%
%       VMWV              volumetric mean wake velocity
%
%       CT                thrust coefficient  == thrust / (0.5*rho*Vs^2*pi*R^2)
%       CQ                torque coefficient  == torque / (0.5*rho*Vs^2*pi*R^3)
%       CP                power  coefficient  == power  / (0.5*rho*Vs^3*pi*R^2)
%
%       KT                thrust coefficient  == thrust / (rho*n^2*D^4)
%       KQ                torque coefficient  == torque / (rho*n^2*D^5)
%
%       EFFYo                 open-water efficiency  == thrust*Vs / torque*(2*pi*n)
%       EFFY == EFFYa     inflow-adapted efficiency  == thrust*Va / torque*(2*pi*n)
%       EFFYw               wake-adapted efficiency  == thrust*Vw / torque*(2*pi*n)
%
%       ADEFo               open-water   actuator disk efficiency
%       ADEFa             inflow-adapted actuator disk efficiency
%       ADEFw               wake-adapted actuator disk efficiency
% 
%       QFo == EFFYo/ADEFo,   open-water   quality factor
%       QFa == EFFYa/ADEFa, inflow-adapted quality factor
%       QFw == EFFYw/ADEFw,   wake-adapted quality factor
%
%       CTH               hub image thrust coefficient (negative, since actually drag) 
%                               CTH == (-hub drag) / (0.5*rho*Vs^2*pi*R^2)
%       CTP,KTP           thrust coefficent of rotor plus hub
%       TAU               thrust ratio == rotor thrust / total thrust 
%
% -------------------------------------------------------------------------
 
function [CT,CQ,CP,KT,KQ, CTH,TAU, Ja,Jw,VMWV, EFFYo, EFFY,ADEFFY,QF,QFo,QFw] = ...
          Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,UADUCT,CD,CoD,G,Z,Js,VMIV,Hub_flag,Rhub_oR,Rhv,CTD) 
            
% ------------------------------------ Volumetric Mean Wake Velocity (VMWV)       
XRtemp   = linspace(Rhub_oR,1,100);
XVAtemp  = pchip(RC, VAC+UADUCT+UASTAR,  XRtemp );      
VMWV     = 2*trapz(XRtemp,XRtemp.*XVAtemp)/(1-Rhub_oR^2);  
% -------------------------------------------------------------------------

% --------------------------------------- Wake-adapted advance coefficients
Ja  = Js*VMIV;
Jw  = Js*VMWV;
% -------------------------------------------------------------------------
      
% -------------------------------------------------------------------------
VSTAR     = sqrt((VAC + UADUCT + UASTAR).^2 + (pi*RC/Js + VTC + UTSTAR).^2);  

sin_BetaI = (UADUCT   + VAC + UASTAR)./VSTAR;
cos_BetaI = (pi*RC/Js + VTC + UTSTAR)./VSTAR;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
CTP = 4*Z*sum((VSTAR.*G'.*cos_BetaI - (1/(2*pi)).*VSTAR.^2.*CoD.*CD.*sin_BetaI)    .*DR);
CQ  = 4*Z*sum((VSTAR.*G'.*sin_BetaI + (1/(2*pi)).*VSTAR.^2.*CoD.*CD.*cos_BetaI).*RC.*DR);
% -------------------------------------------------------------------------

% ----------------- Compute hub effect on thrust coefficient (Kerwin p.181)
if Hub_flag == 1
    CTH = -0.5*(log(1/Rhv)+3)*(Z*G(1))^2;           % Kerwin eqn 260, p.184
else
    CTH = 0;
end
% -------------------------------------------------------------------------


CT    =  CTP+CTH + CTD;     % total thrust coefficient

TAU   = (CTP+CTH)/CT;        % thrust ratio (where CTP+CTH is the net rotor thrust)

CP    = CQ*pi/Js;           % power coefficient based on torque

% KTP = CTP*Js^2*pi/8;      % propeller thrust coefficient
KT    = CT *Js^2*pi/8;      % total     thrust coefficient
KQ    = CQ *Js^2*pi/16;     % total     torque coefficient


% -------------------------------------- Meaningful for propeller case only

  EFFYo = CT/CP;              %   open-water   efficiency
  EFFY  = EFFYo*VMIV;         % inflow-adapted efficiency, EFFY === EFFYa by convention (Kerwin)
% EFFYw = EFFYo*VMWV;         %   wake-adapted efficiency


  ADEFo  =              2*Js/(Js+sqrt(Js^2+TAU*(8/pi)*KT));     % == 2/(1+sqrt(1+TAU*CT       ));     open-water   actuator disk efficiency
  ADEFFY =              2*Ja/(Ja+sqrt(Ja^2+TAU*(8/pi)*KT));     % == 2/(1+sqrt(1+TAU*CT/VMIV^2));   inflow-adapted actuator disk efficiency
% ADEFw  =              2*Jw/(Jw+sqrt(Jw^2+TAU*(8/pi)*KT));     % == 2/(1+sqrt(1+TAU*CT/VMWV^2));     wake-adapted actuator disk efficiency

  QFo   = (1/(4*pi))*(KT/KQ)*(Js+sqrt(Js^2+TAU*(8/pi)*KT)); % == EFFYo/ADEFo,   open-water   quality factor
  QF    = (1/(4*pi))*(KT/KQ)*(Ja+sqrt(Ja^2+TAU*(8/pi)*KT)); % == EFFYa/ADEFa, inflow-adapted quality factor
  QFw   = (1/(4*pi))*(KT/KQ)*(Jw+sqrt(Jw^2+TAU*(8/pi)*KT)); % == EFFYw/ADEFw,   wake-adapted quality factor

% ===================================================== END Forces Function
% =========================================================================
