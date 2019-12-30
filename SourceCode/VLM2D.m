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


% ================================================================= VLM2D.m
% Created  by: J.E. Kerwin
% Modfied  by: H.-L. Chung and R. W. Kimball, 01/28/2007
% Modified by: B. Epps, 3/10/2010
%
% 2D Vortex/Source Lattice with Lighthill Correction (VLM)
%
% -------------------------------------------------------------------------
%
% Inputs:
%   CLI         ideal lift coefficient of foil section
%   ALPHA       == Alpha - AlphaIdeal, where: 
%                                      Alpha      == actual angle of attack
%                                      AlphaIdeal == ideal  angle of attack
%           
%   t0oc        thickness / chord  ratio
%   N           number of panels along the chord, default == 16
%
%   Note: ALPHA must be in RADIANS!
%
%  All length scales are normalized by chord length, c
%  Velocity and vortex sheet strength is normalized by free-stream speed, U
%  Circulation is normalized by U*c
%
%
% Outputs:
%   xp          x/c location of pressure outputs
%   MCPU        -Cp on upper foil surface at x/c 
%   MCPL        -Cp on lower foil surface at x/c 
% 
% -------------------------------------------------------------------------

function [xp, MCPU, MCPL] = VLM2D(N,CLI,ALPHA,t0oc)
%%
% clc,
% N     = 16;
% CLI    = 1;
% ALPHA = 5*pi/180;
% t0oc  = 0.2;
%  
%================================================== Cosine spacing geometry
xv = 1/2 * (1-cos(((1:N)'-1/2)*pi/N));     % vortex  point position, x/c
xc = 1/2 * (1-cos( (1:N)'     *pi/N));     % control point position, x/c
dx = pi  * sqrt(xv.*(1-xv)) / N;           % interval between vrotices

xt = [0;xc];                               % xt(1) == leading edge
xp = [0;xv];                               % xp(1) == leading edge

% ============================================================= Camber Term
% Influence Matrix A (i:control point; j:vortex), 
% for circulation positive clockwise.  A(i,j) = 1/(2*pi*(xv(j)-xc(i))) 
A = zeros(N,N); 

for i = 1:N         
    A(i,:) = 1./(2*pi*(xv-xc(i)));        
end

[B] = MeanLine(xc);               % For 'NACA a=0.8' meanline with unit CLI
%[B,F,Gexact] = MeanLine(xv,xc);  % For 'NACA a=0.8' meanline with unit CLI
% F(j)   == camber            at xv(j) including AlphaIdeal, for unit CLI
% B(i)   == slope of meanline at xc(i) including AlphaIdeal, for unit CLI
% Gexact == exact vortex sheet strength at xv(i)

% for i = 1:N
%    B(i) = CLI*B(i) - ALPHA;        % (Slope  F' scaled for CLI)  -  Alpha
%    F(i) = CLI*F(i);                % (Camber F  scaled for CLI)
% end

B = CLI*B - ALPHA;                   % (Slope  F' scaled for CLI)  -  Alpha

% B = 0*xc - ALPHA; % flat plate case

% ----- Solve the linear system of equations for the point vortex strengths
Gamma = linsolve(A,B);   % Point Vortex Strength (positive clockwise) at xv
G     = Gamma./dx;       % Vortex Sheet Strength (positive clockwise) at xv
% CLnum = 2*sum(Gamma);    % Numerical Lift Coeff via Kutta Joukowski


% % Exact circulation distribution for a flat plate
% xx  = 0:0.001:1;
% GFP = 2*sin(ALPHA)*sqrt((1-xx)./xx);
% figure, plot(xx,GFP,'k',xv,G,'b.'), axis([0 1 0 1])

% ========================================================== Thickness Term
[SourceStrength,dydxSRC,RLE] = Thickness(t0oc, xt);

UT = zeros(N,1);                    % u/U at xc, due to thickness

for i = 1:N                         % (i:control point; j:source point)
    UT(i) = sum(SourceStrength./(2*pi*(xc(i)-xv)));
end

UT = interp1(xc,UT,xv,'pchip','extrap');  % u/U at xv, due to thickness

% =========================================================== Find Pressure
MCPU = zeros(N+1,1);
MCPL = zeros(N+1,1);

% ------------------------------------------------------------ Leading Edge
if RLE < 1e-6
    QU  = ALPHA*sqrt(2/1e-6);              % surface velocity         at LE
else
    QU  = ALPHA*sqrt(2/RLE);               % surface velocity         at LE
end
MCPU(1) = QU^2-1;                          % -Cp on the upper surface at LE
MCPL(1) = QU^2-1;                          % -Cp on the lower surface at LE

% ------------------------------------------------------ Remainder of chord
for i = 1:N
    if dydxSRC(i)>0
        % Scherer/Riegels modified Lighthill correction:
        FLH = 1/sqrt(1+dydxSRC(i)^2);      
    else
        % No leading edge correction beyond point of max thickness:
        FLH = 1;
    end

        % % ORIGINAL Lighthill correction:
        % FLH = sqrt(xp(i+1)/(xp(i+1)+RLE/2));
    
    QU    = (1+UT(i)+1/2*G(i))*FLH;       % Velocity on Upper Surface at xv
    QL    = (1+UT(i)-1/2*G(i))*FLH;       % Velocity on Lower Surface at xv
    
    MCPU(i+1) = QU^2-1;                   % -Cp at xv    
    MCPL(i+1) = QL^2-1;                   % -Cp at xv
end

% figure, 
%     plot(xp,MCPU,'b',xp,MCPL,'b'), hold on, axis([0 1 -1 2])


% % ============================================== Lift and drag coefficients
% Cxv = N * Gamma / pi;           % C/U, leading edge suction parameter at xv
% 
% C0  = interp1(xv,Cxv,0,'linear','extrap');  % C/U LE suc parameter at x = 0
% CS  = C0^2 * pi/2;                          % suction force coefficient
% 
% % % Check versus exact values for flat plate:
% % C0 / (2*sin(ALPHA))
% % CS / (2*pi*(sin(ALPHA))^2)
% % figure, plot(xv,C,'.')
% 
% 
% % Check flat plate case:
% CN = trapz(xv,MCPU(2:end) - MCPL(2:end));
% CL = CN*cos(ALPHA) + CS*sin(ALPHA);
% CD = CN*sin(ALPHA) - CS*cos(ALPHA);
% % CL / (2*pi*sin(ALPHA));
% % CL / CLnum

%%
end % function VLM2D ======================================================
% =========================================================================
% =========================================================================

% ============================================== MeanLine function 01/20/07
% Returns the 'NACA a=0.8' meanline shape and loading for unit CLI for 
% a foil inclined at its ideal angle of attack, AlphaIdeal, where the 
% foild is pinned at the trailing edge.
%
% xv and xc must be non-dimensionalized by chord length, as in x/c
% F(j) == camber            at xv(j) including AlphaIdeal, for unit CLI
% B(i) == slope of meanline at xc(i) including AlphaIdeal, for unit CLI
% Gexact == exact vortex sheet strength at xv(i)
% function [B,F,Gexact] = MeanLine(xv,xc)
function [B] = MeanLine(xc)

a  = 0.8;    % For NACA a=0.8
g  = -1/(1-a) * (a^2*(log(a)/2-1/4)+1/4);
h  =  1/(1-a) * ((1-a)^2*log(1-a)/2 - (1-a)^2/4) + g;

N  = length(xc);
% F      = 0*xc;
% Gexact = 0*xv;

AlphaIdeal = abs(h) / (2*pi*(a+1));

% % --------------------- Camber distribution: foil pinned at trailing edge
% for j = 1:N
%    C1 = max(1 - xv(j),1e-6);
%    CA =     a - xv(j);
%    
%    if abs(CA) < 1e-6
%        CA = CA+1e-5;
%    end
%    
%    P    = 1/2*CA^2*log(abs(CA))-1/2*C1^2*log(C1)+1/4*(C1^2-CA^2);
%    
%    % Camber distribution: foil pinned at trailing edge
%    F(j) = (P/(1-a)-xv(j)*log(xv(j))+g-h*xv(j))/(2*pi*(a+1)) ...
%          + C1*AlphaIdeal;
%     
%    if (xv(j)<=a)
%        Gexact(j) = 1/(a+1);
%    else
%        Gexact(j) = 1/(a+1) * (1-xv(j))/(1-a);
%    end
% end

% ---------------------- Slope of meanline, including ideal angle of attack
C1 = max(1 - xc,1e-6);
CA =     a - xc;

for i = 1:N
    if abs(CA(i))<1e-6
       CA(i) = CA(i)+1e-5;
    end
end

% Slope of meanline: foil pinned at trailing edge
B = ((-CA.*log(abs(CA))+C1.*log(C1))/(1-a) - log(xc) -1-h)/(2*pi*(a+1))...
      - AlphaIdeal; 
     
end % function meanline ===================================================

% ============================================= Thickness function 01/20/07
function [SourceStrength,dydxSRC,RLE] = Thickness(t0oc, xt)

% % -------------------------------------------------------------------------
% % --------------------------------------------- Use "NACA66" thickness form
% PC  = [0.000, 0.010, 0.025, 0.050, 0.100, 0.200, 0.300, 0.400, 0.450,...
%        0.500, 0.600, 0.700, 0.800, 0.900, 0.950, 0.975, 0.990, 1.000];
% TC = [0.0000, 0.1870, 0.2932, 0.4132, 0.5814, 0.8000, 0.9274,...
%        0.9904, 1.0000, 0.9917, 0.9256, 0.7934, 0.5950, 0.3306,...
%        0.1736, 0.0888, 0.0360, 0.0000];
% 
% RLE_CONST = 0.448;               % leading edge radius / chord for t0oc==1
% RLE       = RLE_CONST*t0oc^2;    % leading edge radius / chord
% TRLE      = 2*sqrt(2*RLE_CONST);
% 
% PSQ = sqrt(PC); % square root stretched coordinate for spline interpolation
% XSQ = sqrt(xt); % square root stretched coordinate for spline interpolation
% 
% thickness      = t0oc*spline(PSQ,TC,XSQ);       % thickness at xt == [0;xc]
% % -------------------------------------------------------------------------
  
% -------------------------------------------------------------------------
% ------------------------------------------ Use NACA 65A010 thickness form
% x-position in percent chord
% PC = [0 .5 .75 1.25 2.5 5 7.5 10 15 20 25 30 35 40 45 50 ...
%       55 60 65 70 75 80 85 90 95 100]./100;
PC = [  ...     
     0  0.0050  0.0075  0.0125    0.0250    0.0500    0.0750    0.1000 ...
0.1500  0.2000  0.2500  0.3000    0.3500    0.4000    0.4500    0.5000 ...
0.5500  0.6000  0.6500  0.7000    0.7500    0.8000    0.8500    0.9000 ...
0.9500  1.0000];

% thickness/2 (in percent chord) normalized by the maximum thickness/2
% TC = [0 .765 .928 1.183 1.623 2.182 2.65 3.04 3.658 4.127 ...
%           4.483 4.742 4.912 4.995 4.983 4.863 4.632 4.304     ...
%           3.899 3.432 2.912 2.352 1.771 1.188 .604 .021]/4.9950;
TC = [          0 0.153153153153153 0.185785785785786 0.236836836836837 ...
0.324924924924925 0.436836836836837 0.530530530530530 0.608608608608609 ...
0.732332332332332 0.826226226226226 0.897497497497497 0.949349349349349 ...
0.983383383383383 1.000000000000000 0.997597597597597 0.973573573573574 ...
0.927327327327327 0.861661661661662 0.780580580580581 0.687087087087087 ...
0.582982982982983 0.470870870870871 0.354554554554555 0.237837837837838 ...
0.120920920920921 0.004204204204204];

RLE = 0.00639 * t0oc^2/(2*0.04995)^2;
% -------------------------------------------------------------------------
      
PSQ = sqrt(PC); % square root stretched coordinate for spline interpolation
XSQ = sqrt(xt); % square root stretched coordinate for spline interpolation


% 6/1/2012 BEPPS: pchip(PSQ,TC,XSQ) and pchip(PC,TC,xt) agree to witin 0.3 percent.
%                 Jake uses spline(PSQ,TC,XSQ), but spline produces "wiggles" in the curve.

thickness      = t0oc*pchip(PSQ,TC,XSQ);       % thickness at xt == [0;xc]
% -------------------------------------------------------------------------


% % -------------------------------------------------------------------------
% % ------------------------------------------------ Parabolic thickness form
% tot0 = (1-(2*(xt-0.5)).^2);
% 
% fig; plot(xt,tot0)
% 
% RLE = 0.00639/(2*0.04995)^2 * t0oc^2;
%     
% thickness = tot0*t0oc;
% % -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
SourceStrength = diff(thickness);    % Point Source strength (source at xv)

% ---------- Half of dtdx is used in forming the Scherer/Riegels correction
dydxSRC = 0.5*diff(thickness)./diff(xt);       
          

% figure, hold on
% plot(PC,  t0oc*TC/2,'b.-',PC,-  t0oc*TC/2,'b.-') 
% plot(xt,thickness/2,'.-r',xt,-thickness/2,'.-r') 
% theta = 0:0.01:2*pi;
% plot(RLE+RLE*cos(theta),RLE*sin(theta),'k')
% axis equal
% %%
end % function thickness ==================================================

