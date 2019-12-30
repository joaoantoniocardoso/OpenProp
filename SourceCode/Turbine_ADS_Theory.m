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
%
% =========================================================================
% ==================================================== Turbine_ADS_Theory.m
% Last Modified: 11/2/2011, Brenden Epps  
%
% -------------------------------------------------------------------------
%
% Compute turbine actuator disc theory performance, with swirl and viscous 
% losses (Stewart, 1976).
%
% Inputs:
%   L  (lambda)             tip speed ratio
%   Z                       number of blades
%   CDoCL = CD / CL         viscous effects         
%   RC    = r/R             radii for load distributions
%
% Outputs:
%   CPBetz                  power coefficient
%   BetzGRC                 circulation at given RC
%   BetzRC                  r/R for turbine theory values below
%   BetzG                   circulation at BetzRC
%   BetzUA,BetzUT           (axial,tangential) velocities at BetzRC
%   BetzTAN                 tangent of inflow angle at BetzRC
% -------------------------------------------------------------------------
%
%        [CPBetz,BetzRC,BetzG,BetzUA,BetzUT,BetzTAN, BetzGRC] = Turbine_ADS_Theory(L,Z,CDoCL,RC)
%
% -------------------------------------------------------------------------

function [CPBetz,BetzRC,BetzG,BetzUA,BetzUT,BetzTAN, BetzGRC] = Turbine_ADS_Theory(L,Z,CDoCL,RC)

%%
% clr,
% L=10; Z=100;
% CDoCL = 0.10;
% RC=0.2:0.05:1;
 
    phi = [pi/3       :-0.01   :pi/3-0.6,...
           pi/3-0.601 :-0.001  :pi/3-0.7,...
           pi/3-0.7001:-0.0001 :0.0001];  % inflow angle  (min phi for CDoCL > 0.1 is 0.0997)

    a = 1./(2 + CDoCL*sec(phi).^2./(2*tan(phi)-CDoCL) + sec(phi).*( tan(phi)./(tan(phi)-CDoCL)  + CDoCL^2*sec(phi).^2./(4*(tan(phi)-CDoCL).^2)  ).^0.5 );
    
    x = (1 - a.*sec(phi).^2)./tan(phi);


    % Error checking to remove bad values
    ind = min([ length(x), find( diff(x) < 0 , 1) ]);

    phi = phi(1:ind);
    a   =   a(1:ind);
    x   =   x(1:ind);
    % ---------------

    
    % ap = a.*tan(phi)./x;
    % BetzRC_times_lambda = x;
    
    BetzUA     = -a;
    BetzUT     =  a.*tan(phi);    % = ap.*x;
    BetzTAN    =     tan(phi);    % = ap.*x./a; 
    BetzVSTAR  = sqrt((1+BetzUA).^2+(x+BetzUT).^2);
    BetzRC     = x/L;
    BetzG      = (-1/Z) * 2*BetzRC.*BetzVSTAR.*(  a./(1-a) .* sin(phi).*tan(phi) );  % (-1/Z) factor to be consistent with OpenProp sign convention
  %
  % BetzGwrong = 2*BetzRC.*BetzVSTAR.*(1-cos(phi)) * (-1/Z);  % NOTE BUG FIX:  This is Stewart eqn 14, which is only correct for CDoCL == 0 case!!!
  %  
%     fig; 
%         plot(phi, 1-cos(phi), 'k')
%         plot(phi, a./(1-a) .* sin(phi).*tan(phi), 'r')
%    
%     fig;
%         plot(BetzRC,BetzG,'k',BetzRC,BetzG2,'r')
        
        
    
    BetzGRC   = interp1(BetzRC,BetzG,RC);
    
    CPp     = 4*a.*(1-a).*x.*(tan(phi)-CDoCL);
    CP      = (2./x.^2).*cumtrapz(x,CPp.*x);

    xtemp   =  x(find(CP > 0 & CP < 0.6));
    CPtemp  = CP(find(CP > 0 & CP < 0.6));
    CPBetz  = -interp1(xtemp,CPtemp,L); % - sign to be consistent with OpenProp sign convention
    
    
    
    % Only report values of BetzRC less than 1
    indices = find(BetzRC <= 1 & BetzRC > 0);
    
    BetzRC  =  BetzRC(indices);
    BetzG   =   BetzG(indices);
    BetzUA  =  BetzUA(indices);
    BetzUT  =  BetzUT(indices);
    BetzTAN = BetzTAN(indices);
    
%     plot(BetzRC,BetzG,'b')

    % fig; plot(BetzRC,BetzUA,BetzRC,BetzUT)
    % fig; plot(BetzRC,BetzG,'k')    
    % fig; plot(BetzRC,BetzTAN,'r--');
 
end
% -------------------------------------------------------------------------