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

% -------------------------------------------------------------------------
% Created: 3/20/11, Brenden Epps
%
%
% Augment thickness and chord distributions to satisfy ABS thickness rule, 
% while keeping the t0oc distribution the same.
%
% Find required thickness at r/R = 0.25
% Fit a spline to given thickness distribution.
% Find point (R0,T0) on spline where line connecting (R0,T0) and (0.25,t0oD25) is tangent to spline at (R0,T0).
% Use linear distribution between (0.25,t0oD25) and (R0,T0), and use spline for radii larger than R0.
%
% INPUTS:
%   RC          r/R   radii
%   t0oD        t0/D  thickness distribution   at RC
%   t0oc        t0/c  thickness to chord ratio at RC
%   TANtheta    [ ]   pitch angle distribution at RC,  TANtheta = tand( atand(pt.design.TANBIC) + pt.design.alphaI );  % tan(pitch angle)
%   R           [m]   radius
%   rho         [kg]  water density
%   Z           [ ]   number of blades
%   CP          [ ]   power coefficient at rated speed
%   Vs          [m/s] ship speed        at rated speed
%   N           [RPM] rotation rate     at rated speed
%
% OUTPUTS:
%   CoD2         c/D  modified chord     distribution at RC
%   t0oD2       t0/D  modified thickness distribution at RC
%
% -------------------------------------------------------------------------

function [CoD2, t0oD2] = ABS_Blade(RC,t0oD,t0oc,TANtheta,  CP,R,rho,Vs,N,Z)

D = 2*R;

% -------------------------------------------------------------------------
% ABS thickness: compute required CoD at RC = 0.25 for given t0oc distribution

        % Compute approximate t0/D at blade root using approximation to the ABS rule 
        ABS_w = 7.5; ABS_f = 2.62;
        
        % PkW = (CP * 0.5 * rho * Vs^3 * pi * R^2) / 1000; % power, in kW
        P = (CP * 0.5 * rho * Vs^3 * pi * R^2); % power, in W
        
        PoD70 = pi * 0.70 * interp1(RC,TANtheta,0.70,'pchip','extrap');  % blade pitch angle approximately the wake pitch angle
        PoD25 = pi * 0.25 * interp1(RC,TANtheta,0.25,'pchip','extrap');  % blade pitch angle approximately the wake pitch angle
        
        t0oc25 = interp1(RC,t0oc,0.25,'pchip','extrap');  % thickness to chord ratio at r/R = 0.25
        
        % Constants in ABS rule calculation: T0 == 1.025 * 337 / sqrt(0.1) == 1.0923e+03
        % CoD25 = ( (1.0923e+03)^2 * ( (1 + 6/PoD70 + 4.3*PoD25)/(1 + 1.5*PoD25) * (PkW/1000)/(ABS_f*D*N*Z) ) / 1000^2 / D^2 / t0oc25^2 )^(1/3);
        % 
        % CoD25 = ( (1.0923e+03)^2 * ( (1 + 6/PoD70 + 4.3*PoD25)/(1 + 1.5*PoD25) * (P)/(ABS_f*D^3*N*Z) ) / 1000^4 /  t0oc25^2 )^(1/3);
        % 
        % 
        % T1 = (1.0923e+03/1000^2) * sqrt( (1 + 6/PoD70 + 4.3*PoD25)/(1 + 1.5*PoD25) * (P)/(ABS_f*D^3*N*Z) ) 
        % CoD25 = ( T1^2 /  t0oc25^2 )^(1/3);
        % 
        % 
        % T1 = (1.0923e+03/1000^2) * sqrt( (1 + 6/PoD70 + 4.3*PoD25)/(1 + 1.5*PoD25) * (P)/(ABS_f*D^3*N*Z) ) 
        % CoD25 = ( T1 /  t0oc25 )^(2/3);
        %
        % t0oD25 = t0oc25*CoD25;   % thickness at RC = 0.25 for given t0oc
        %
        % ---
        % T1 = (1.0923e+03/1000^2) * sqrt( (1 + 6/PoD70 + 4.3*PoD25)/(1 + 1.5*PoD25) * (P)/(ABS_f*D^3*N*Z) ) 
        % t0oD25 = T1^(2/3) * t0oc25^(1/3);   % thickness at RC = 0.25 for given t0oc
        %
        % ---
        t0oD25 = ( (1.0923e+03/1e6)^2 * ( (1 + 6/PoD70 + 4.3*PoD25)/(1 + 1.5*PoD25) * (P)/(ABS_f*D^3*N*Z) ) * t0oc25)^(1/3);   % thickness at RC = 0.25 for given t0oc
         
% -------------------------------------------------------------------------


%   fig; 
%     plot(RC,t0oD,'k.-')
% 
%     plot(0.25,t0oD25,'k.')      

%
% Given t0oD25, expand the t0oD distribution to fit that:
% Expand t0oD as opposed to CoD, because the maximum t0oD occurrs inboard of the maximum CoD   


% -------------------------------------------------------------------------
% Fit spline to all t0oD data,
% Find point (R0,T0) where tangent vector points through (0.25,t0oD25)
% Linear interpolate between (R0,T0) and (0.25,t0oD25)
% ------------------------------------------------------------------------- 
% ----------------------------------------------------------- Spline inputs   
Md = length(RC);  % number of  spanwise data sites, Mc == m+1
m  = Md - 1;      % Mc spline basis functions / control points
k  = 4;           % polynomial order (k == 4 for cubic spline)
Mk = k+m+1;       % number of  spanwise knots

% ------------------------- Find spline parameters using centripital method
dseg = ( diff(RC').^2 + diff(t0oD').^2 ).^(1/4);

dtot = sum(dseg);

ubar = [0;cumsum(dseg)]/dtot;    % size [Md,1]


% -------------------------------- Find spline knots using averaging method
uknot         = zeros(Mk,1);  % spanwise knot sequence

uknot(m+2:Mk) = ones(k,1);

for j = 1:m+1-k
    uknot(k+j) = sum(ubar(1+(j:j+k-2)))/(k-1);
end


% ----------------------------------- Evaluate the B-spline basis functions
[BC, DC] = Bspline_basis(ubar,uknot,k);    % size [Md,Md], BC(ubar,uspline)


% ------------------------------- Solve linear system for spline amplitudes
% Let: 
%       xr2(ubar) = BC(ubar,uspline) * AR(uspline)
% and
%       xc2(ubar) = BC(ubar,uspline) * AT(uspline)
%
% where:
%
%   size(RC)   = [Md,1], data matrix 
%   size(t0oD) = [Md,1], data matrix 
%   size(BC)   = [Md,Md], Md == length(ubar) == number of u-splines 
%   size(AR)   = [Md,1], xr2 spline amplitudes 
%   size(AT)   = [Md,1], xc2 spline amplitudes 
%
% --------------------------------- Solve for spline amplitudes, P(uspline)
AR = linsolve(BC,RC');
AT = linsolve(BC,t0oD');

%
% % -------------------------------------------------------- Check solution
% RC'   - BC * AR
% t0oD' - BC * AT
% -------------------------------------------------------------------------



% % ---------------------------------------------------- Find tangent vectors
% % Tangent vectors (TR,TT) at (ubar)
% % Note: negative sign reverses the direction of the vector!
% TR  = -DC * AR;
% TT  = -DC * AT;
% 
% [TR,TT,junk] = NormalizeVector(TR,TT,0*TT);
% 
% % quiver(RC',t0oD',TR,TT)
% % -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Find spline parameter u0 such that minus tangent vector points through (0.25,t0oD25), 
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ------------------------------------- Evaluate spline on finer resolution
Mdd = 1001;

ubar2 = linspace(0,1,Mdd)';

[BC, DC] = Bspline_basis(ubar2,uknot,k);    % size [Md,Md], BC(ubar,uspline)

xr2 = BC * AR;
xt2 = BC * AT;

% plot(xr2,xt2,'g')

% save temp
% disp('paused...')
% pause,

if pchip(xr2,xt2,0.25) < t0oD25

    shift = find(xr2 > 0.25, 1) - 1;

    % plot(xr2(shift+1:end), xt2(shift+1:end), 'b')

    % Minus tangent vectors (TR,TT) at (ubar2)
    TR2  = -DC * AR;
    TT2  = -DC * AT;

    [TR2,TT2,junk] = NormalizeVector(TR2,TT2,0*TT2);


    i2 = min( [ find((t0oD25 - xt2(shift+1:end))./(0.25 - xr2(shift+1:end)) - TT2(shift+1:end)./TR2(shift+1:end) > 0, 1)+ shift , Mdd] );

     R0 = xr2(i2);
     T0 = xt2(i2);
    TR0 = TR2(i2);
    TT0 = TT2(i2);

    % plot(R0,T0,'.r')
    % quiver(R0,T0,TR0,TT0,'r')

    % -------------------------------------------------------------------------
    % Form new thickness distribution (RC,t0oD2)
    % Keep t0oD data for RC > R0, linearly interpolate for RC < R0
    t0oD2 = t0oD;

    t0oD2(RC<R0) = interp1([0.25,R0],[t0oD25,T0],RC(RC<R0),'linear','extrap');

else
    % t0oD at r/R = 0.25 is greater than required t0oD25, so keep as is
    t0oD2 = t0oD;

end



% plot(RC,t0oD2,'m')

CoD2 = t0oD2 ./ t0oc;



% RG = linspace(0.25,1,100);
% 
% CoDG2 = InterpolateChord(RC,CoD2,RG);



% % fig: 
% PlotBlade(XR,XCoD,'r')
% hold on,
% plot(RG,CoDG/2,'g'), plot(RG,-CoDG/2,'g')
% 
% plot(0.25,CoD25/2,'k.'), plot(0.25,-CoD25/2,'k.')
% 
% plot(RC,CoD2/2,'b.'), plot(RC,-CoD2/2,'b.')
% plot(RG,CoDG2/2,'b'), plot(RG,-CoDG2/2,'b')

end



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% This function scales the given vector to have a norm of 1.

function [vx vy vz] = NormalizeVector(vx,vy,vz)


veclength = sqrt(vx.^2 + vy.^2 + vz.^2);

vx = vx./veclength;
vy = vy./veclength;
vz = vz./veclength;

end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------