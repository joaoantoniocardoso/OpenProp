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
% ======================================================== Duct_Influence.m
% Created: J.M. Stubblefield, 2008    (originally named ductVort)
% Last Modified: 11/2/2011, Brenden Epps
%
% -------------------------------------------------------------------------
% This function:
%   1) Forms a discrete vortex ring representation of a duct with 
%           -- chord length Cduct_oR and axial offset Xduct_oR,
%           -- radius Rduct_oR, and
%           -- a NACA a=0.8 meanline 
%
%   2) Calculates axial velocity duct influence function (UADIF)
%      on propeller lifting line control points
%
% -------------------------------------------------------------------------
% Model parameters:
%
%     Rduct_oR                  duct radius / propeller radius
%     Cduct_oR                  duct chord  / propeller radius
%     Xduct_oR                  location of duct mid-chord downstream of propeller (Xduct = 0 by default)
%
%     XdRING            [1,Nd], x/R location of each vortex ring downstream of propeller
%     VARING            [1,Nd], (axial free-stream velocity at duct)/Vs
%
%     GdRING            [1,Nd],  fraction of total non-dimensional duct circulation, sucht that sum(GdRING) = 1 
%                                (i.e. non-dimensional circulation per unit Gd)
%     Gd                         total non-dimensional circulation about the duct, Gd == Gamma_d/(2*pi*R*Vs),
%                                such that the non-dimensional circulation of ring n is Gd*GdRING(n).
%     Gamma_d [m^2/s]            total     dimensional circulation about the duct [m^2/s]
%
%     CDd                       section drag coefficient for the duct
%     CTDDES                    desired duct thrust coefficient
%     CTD                       duct thrust coeff (viscous drag included) with total duct circulation of Gd 
%
% Influence of propeller on duct:
%
%   (DAHIF ,DRHIF)     [Nd,Mp], (axial/radial horseshoe influence functions of prop on duct)*(2*pi*R)
%   (UARING,URRING)    [1,Nd],  (axial/radial velocity induced at duct (XdRING,Rduct_oR) by prop) / Vs
%
%
% Influence of duct on propeller:
%
%   UADUCT      [1,Mp]  (axial velocity induced on PROP by duct)/Vs
%   UADIF       [1,Mp]   axial velocity induced on PROP by duct per unit Gd, 
%                        i.e. non-dimensional Duct Influence Function
%                         == 2*pi*R * (duct influence function with unit dimensional duct circulation)
%
% -------------------------------------------------------------------------
%
% Inputs:
%     Rduct_oR,Cduct_oR,Xduct_oR
%     RC                        radius of control points on lifting line / R
%
% Outputs:
%     XdRING,GdRING,UADIF
%
% =========================================================================

function [XdRING,GdRING,UADIF] = Duct_Influence(Rduct_oR,Cduct_oR,Xduct_oR,RC)

Mp = length(RC);  % number of panels

% % Test code:
% R = 1;
% Rduct_oR = 1;
% Mp = 10;
% RV = linspace(0.2*R,R,Mp+1);
% RC = 0.5*(RV(1:Mp)+RV(2:Mp+1));
% Cduct_oR = R;
% Nd = 10;

% ----------------------------------------- Setup vortex ring axial spacing
% Note that the code is numerically unstable for Nd > 20 and generally 
% slow for Nd > 10. Nd must be even.
Nd = 12;             % number of vortex rings 
% if rem(Nd,2)~=0      % ensures Nd is even
%     Nd = Nd+1;
% end


dS  = 1/Nd;          % (even spacing between vortex rings) / Cduct_oR
hdS = 0.5*dS;        % half of dS

% -------------------------------------------------------------------------
% Compute the circulation on the vortex rings which each represent 
% a section of length dS located at position XdRING of a NACA a=0.8 
% mean line (L.E.=0.0, T.E.=1.0) that has unit circulation 1 [m^2/s].
%
XdRING  = zeros(1,Nd);
GdRING  = zeros(1,Nd);

% Note that since Gamma_d*GdRING(n) is the circulation of ring (n), then
%    the circulation distribution of the rings, GdRING, is actually unitless!
%    If Gamma_d is scaled or normalized, then GdRING remains unitless and
%    always gives the correct circulation distribution.
for n=1:Nd
    XdRING(n) = (n-1)*dS+hdS;

    x2 = XdRING(n) + hdS;
    x1 = XdRING(n) - hdS;
    if x2 <= 0.8
        GdRING(n) = dS/0.9;
    elseif x1 >= 0.8
        y1 = 1.0 - (x1 - 0.8)/0.2;
        y2 = 1.0 - (x2 - 0.8)/0.2;
        GdRING(n) = dS*0.5*(y1 + y2)/0.9;
    else
        y2 = 1.0 - (x2 - 0.8)/0.2;
        front    = 0.8 - x1;
        back     = 0.5*(1.0 + y2)*(x2 - 0.8);
        GdRING(n) = (front + back)/0.9;
    end
end

LED      = -(Nd/2)*dS;                  % X position of "leading edge" vortex ring
XdRING = XdRING + LED;                  % X position of vortex rings / Cduct_oR
XdRING = XdRING * Cduct_oR  + Xduct_oR; % X position of vortex rings / R


% % plot(XdRING,GdRING,'*')

% ------------------------------------------------------------------------- 
% Calculate duct influence function (UADIF) on the UASTAR
% axial velocity at propeller lifting line control points
    % Note: No tangential influence
    % Note: Radial influence does not create a force on radial lifting line
    
UADIF = zeros(1,Mp);            % axial influence of unit strength

for m=1:Mp                      % for each control point on lifting line    
    for n=1:Nd                  % cycle thru all vortex rings on duct
        UAD      = vRing(XdRING(n),Rduct_oR,0,RC(m),GdRING(n));
        
        UADIF(m) = UADIF(m) + UAD;
    end
end

% The discrete vortex ring formulation breaks down near the duct, so 
% extrapolate UADIF for RC > 0.9*Rduct/R
indices = find( RC <= 0.9*Rduct_oR );
UADIF   = interp1(RC(indices),UADIF(indices),RC,'linear','extrap');


% Since GdRING is non-dimensionalized with respect to Gd, then
% UADIF represents the axial flow velocity induced per unit Gd, which
% is the correct influence function.

end  % ========================================= END ductInfluence Function
% =========================================================================



% =========================================================================
% ===================================================== vRing Function
% Created: J.M. Stubblefield, 2008
% Edited:  B.P. Epps, 2010
%
% Returns axial and radial velocity at field point (XF,RF) induced by 
% vortex ring of NON-DIMENSIONAL circulation Gring = Gamma/(2*pi*R*Vs) 
% at x-axis location (Xring) and radius (Rring).
%
% NOTE: The inputs MUST be non-dimensionalized as follows.  For dimensional
%       inputs, the output velocity must be multiplied by 2*pi*Vs
%
% Example: for R = 1 m, Vs = 1 m/s, Gamma = 1 m^2/s
%              Xring = XF = 0, Rring = 1, RF = 0 (i.e. center of ring)
%          then Gring = 1/(2*pi)
%          so  [UA,UR] = vRing(0,1,0,0,1/(2*pi))
%          check dimensional output: UA*Vs = Gamma / (2*Rring*R)
%
% Variables:
%     R  [m]     propeller radius
%     Vs [m/s]   ship speed
%     Gring      vortex ring circulation / (2*pi*R*Vs)
%     Xring      x/R, axial location of vortex ring / R
%     Rring      r/R, radius         of vortex ring / R
%     XF         x/R, field point at which velocity is induced
%     RF         r/R, field point at which velocity is induced
%
% Returns: 
%     UA         axial  velocity / Vs     at field point
%     UR         radial velocity / Vs     at field point
%
% Ref: Kuchemann and Weber, Aerodynamics of Propulsion p 305.
% =========================================================================

function [UA,UR] = vRing(Xring,Rring,XF,RF,Gring)

if Rring == 0               %stops function if Rring = 0
    UA = 0;
    UR = 0;
    return
end

if XF == Xring && RF == Rring          % stop if field point on vortex ring
    UA = 0;
    UR = 0;
    return
end

% --------------------------------------- Non-dimensional coordinates (x,r)
x = (XF-Xring)/Rring;                       %x/r' from Kuchemann
r =  RF       /Rring;                       %r/r' from Kuchemann

%Elliptic integral method (Kuchemann p. 305)
%uses parameter k where k^2 = m for elliptic integrals

k     = sqrt(4*r/(x^2+(r+1)^2));
% [K,E] = ellipke(k^2);
[K,E] = elliptic12(pi/2,k^2);

UA = Gring/(Rring)/sqrt(x^2+(r+1)^2)*(K-(1+2*(r-1)/(x^2+(r-1)^2))*E);

if r==0
    UR = 0;
else
    UR = Gring/(Rring)*(-x)/r/sqrt(x^2+(r+1)^2)*(K-(1+2*r/(x^2+(r-1)^2))*E);
end

end  % ================================================= END vRing Function
% =========================================================================