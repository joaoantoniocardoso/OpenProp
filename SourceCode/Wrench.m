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
% ================================================================ Wrench.m
% Last update: 11/2/2011 Brenden Epps
%
% -------------------------------------------------------------------------
% This function evaluates the velocity induced on a point on the key
% lifting line due to a unit-strength helical trailing vortex shed from
% each of the Z lifting lines. This function returns non-dimensional 
% induction factors: 2*pi*R*(dimensional u_bar)
%
% Reference:
%
% Lerbs, H.W. ``Moderately Loaded Propellers with a Finite Number of Blades 
%               and an Arbitrary Distribution of Circulation."  
%               Trans. SNAME, v. 60, 1952.
% Wrench, J. W. ``The Calculation of Propeller Induction Factors," Technical 
%                 Report 1116, David Taylor Model Basin, February, 1957.
%
% -------------------------------------------------------------------------
% Sign convention:
%
% The circulation of a trailing vortex is defined positive when it is 
% directed AWAY from the lifting line (i.e. downstream).
%
% UASTAR is defined positive in the direction direction
% UTSTAR is defined positive in the direction of the apparent inflow, omega*r
%
% -------------------------------------------------------------------------
% Variables:
%       Z         [ ],    number of blades
%       tan_betaW [ ],    tangent of the pitch angle of helical wake trail
%       rc        [ ],    radius of control point  / propeller radius
%       rv        [ ],    radius of helical vortex / propeller radius
%
%       u_barA    [ ],    2*pi*R*u_bar velocity in the axial      direction
%       u_barT    [ ],    2*pi*R*u_bar velocity in the tangential direction
%       y,y0,U,F1,F2,     auxilary variables
%
% -------------------------------------------------------------------------
    

function [u_barA, u_barT] = Wrench(Z,tan_betaW,rc,rv)

% % ------------------- Enable this to check for infinite bladed propellers
% if Z > 20   % Return infinite blade result if Z > 20.
%     if rc == rv
%         u_barA = 0;
%         u_barT = 0;
% 
%     elseif rc < rv
%         u_barA = Z/(2*rv*tan_betaW);    
%         u_barT = 0;                    
% 
%     else % rc > rv
%         u_barA = 0;               
%         u_barT = Z/(2*rc);             
%     end
%     return;
% end

y  = rc/(rv*tan_betaW);
y0 = 1 /    tan_betaW;
U  = ((y0*(sqrt(1+y ^2)-1))*exp(sqrt(1+y^2)-sqrt(1+y0^2))/...
      (y *(sqrt(1+y0^2)-1)))^Z;

if rc == rv
    u_barA = 0;
    u_barT = 0;

elseif rc < rv

    F1     = -(1/(2*Z*y0)) * ((1+y0^2)/(1+y^2))^0.25 * ((1/(U^-1-1)) + (1/(24*Z)) * (((9*y0^2+2)/(1+y0^2)^1.5)+((3*y^2-2)/(1+y^2)^1.5)) * log(1+(1/(U^-1-1))) );

    u_barA = Z/(2*rv*tan_betaW) -  y*Z^2*y0*F1/rc;   
    u_barT =                         Z^2*y0*F1/rc;        

else % rc > rv

    F2     =  (1/(2*Z*y0)) * ((1+y0^2)/(1+y^2))^0.25 * ((1/(U-1))    - (1/(24*Z)) * (((9*y0^2+2)/(1+y0^2)^1.5)+((3*y^2-2)/(1+y^2)^1.5)) * log(1+(1/(U-1))) );

    u_barA =           - y*Z^2*y0*F2/rc;       
    u_barT =  Z/(2*rc) +   Z^2*y0*F2/rc; 

end
% ===================================================== END Wrench Function
% =========================================================================