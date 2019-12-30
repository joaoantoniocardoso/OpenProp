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
% =========================================================== Duct_Thrust.m
% Created by: J.M. Stubblefield, 2008
% Last Modified: 11/2/2011, Brenden Epps
%
% -------------------------------------------------------------------------
% Model parameters:
%
%     Rduct_oR                  duct radius / propeller radius
%     Cduct_oR                  duct chord  / propeller radius
%     Xduct_oR                  location of duct mid-chord downstream of propeller (Xduct = 0 by default)
%
%     XdRING            [1,Nd],  x/R location of each vortex ring downstream of propeller
%     VARING            [1,Nd],  (axial free-stream velocity at duct)/Vs
%    (UARING,URRING)    [1,Nd],  (axial/radial velocity induced at duct (XdRING,Rduct_oR) by prop) / Vs
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
%     Vs [m/s]          ship speed (reference speed)
%     R [m]             radius of propeller
%
% Returns:
%   (1) The duct thrust coefficient, CTD, given duct circulaiton, Gd, and
%   (2) The required circulation, GdDES, to provide the desired duct thrust
%       coefficient, CTDDES.
% -------------------------------------------------------------------------

function [CTD,GdDES] = Duct_Thrust(XdRING,Rduct_oR,VARING,UARING,URRING,GdRING,Gd,CDd,CTDDES)

% ---------------------------------------------------------- Compute forces
delS = abs(XdRING(2) - XdRING(1));  % (vortex ring spacing)/R,   (linear spacing assumed)


% CTD_inviscid     = sum( 4 * (-URRING) .* (Gd*GdRING) *  (2*pi*Rduct_oR) ); % inviscid Kutta-Joukowski thrust
  CTD_inviscid_oGd = sum( 4 * (-URRING) .* (   GdRING) *  (2*pi*Rduct_oR) );

  CTD_inviscid     =     CTD_inviscid_oGd * Gd;
  
  CTD_viscous      = sum(- (VARING+UARING).^2 * CDd * delS * (2*Rduct_oR) ); % viscous drag on duct so negative thrust

  CTD              = CTD_inviscid + CTD_viscous;                             % total duct thrust coefficient


% ------------ Scale duct circulation so that duct provides required thrust
% CTD    = CTD_inviscid    + CTD_viscous
% CTDDES = CTD_inviscidDES + CTD_viscous
% GdDES  = CTD_inviscidDES / (CTD_inviscid/Gd)

GdDES = (CTDDES - CTD_viscous)/CTD_inviscid_oGd;

% ================================================= END ductThrust Function
% =========================================================================

