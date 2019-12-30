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
% Last modified: 11/17/2011 Brenden Epps
% =========================================================================

%--------------------------------------------------------------------------
% 10/21/2011 Brenden Epps
%
% Given: Vs, n, and propeller characteristics, estimate maximum blade stress
%
%   Vs      [m/s], ship speed
%   n       [rev/s], propeller rotation speed
%
%   D       [m], propeller diameter
%   Dhub    [m], propeller hub diameter
%
%   Propeller geometry: r/R, c/D, theta (pitch), skew, rake, t0/D
%   
%--------------------------------------------------------------------------

function [MAX_STRESS, blade_stress] = Stress_Analysis(pt)

warning off
if isfield(pt,'i'), i = pt.i; else i = pt.input; end
if isfield(pt,'d'), d = pt.d; else d = pt.design; end

input = i;

% -------------------------------------------------------- Unpack variables

% --------------------------------------------------------- Geometry inputs
Z = input.Z;           % number of blades

if isfield(input,   'Mp'), Mp    = input.Mp;    else  Mp    = 20;   end
if isfield(input,   'Np'), Np    = input.Np;    else  Np    = 20;   end
if isfield(input,   'Vs'), Vs    = input.Vs;    else  Vs    = 1;    end % m/s

if     isfield(input,'D'),   R = input.D/2;  
elseif isfield(input,'R'),   R = input.R;
else                         R = 1;
end

if     isfield(input,'Dhub'),   Rhub = input.Dhub/2;
elseif isfield(input,'Rhub'),   Rhub = input.Rhub;
else                            Rhub = 0.2*R;
end

if isfield(input,'Rcirc'), Rcirc = input.Rcirc; else  Rcirc = Rhub; end % m
if isfield(input,'Rroot'), Rroot = input.Rroot; else  Rroot = Rhub; end % m

if Rcirc < Rhub,  disp('ERROR: Rcirc must be >= Rhub.  Setting Rcirc = Rhub.' ), Rcirc = Rhub;  end
if Rroot < Rcirc, disp('ERROR: Rroot must be >= Rcirc. Setting Rroot = Rcirc.'), Rroot = Rcirc; end
% -------------------------------------------------------------------------

% ------------------------------------------------------ Tip speed ratio, L
if i.Propeller_flag == 1
    Js    = input.Js;
    L     = pi/Js;
else
    L     = input.L;
    Js    = pi/L;
end

D = 2*R;          % [m]
n = Vs/(Js*D);    % [rev/s]
N = 60*n;         % [RPM]
% -------------------------------------------------------------------------

% ------------------------------------------------------- Cavitation inputs
if isfield(input,'rho'),  rho = input.rho;  else rho  = 1000;  end % kg/m^3
% -------------------------------------------------------------------------


rho_p           = 8200;     % [kg/m^3] density of Ni-Al-Bronze propeller
allowable_stress= 0.5 * 430e6; % [N/m^2] 0.5 * yield stress of Ni-Al-Bronze propeller


XR              = i.XR;
skew0           = i.skew0;  % [deg]
rake0           = i.rake0;


% '------ Design geometry ------'
RC         = d.RC;          % [1 x Mp] control point radii
RV         = d.RV;          % [1 x Mp+1] vortex point radii
DR         = d.DR;          % [1 x Mp] difference in vortex point radii
CoD        = d.CoD;         % [1 x Mp] 
t0oD       = d.t0oD;
Rhub_oR    = d.Rhub_oR;     % [1 x 1]

Mp = length(RC);

% '------ Design state ------'
G          = d.G';          % [1 x Mp] circulation distribution
VAC        = d.VAC;         % [1 x Mp] 
VTC        = d.VTC;         % [1 x Mp] 
UASTAR     = d.UASTAR;      % [1 x Mp] 
UTSTAR     = d.UTSTAR;      % [1 x Mp] 
VSTAR      = d.VSTAR;       % [1 x Mp] 
TANBIC     = d.TANBIC;      % [1 x Mp] 
VMIV       = d.VMIV;        % [1 x 1]

% '------ 2D section performance ------'
CL         = d.CL;          % [1 x Mp] 
CD         = d.CD;          % [1 x Mp] 
CP         = d.CP;


% 0 == do not display plots, 1 == display plots
if isfield(i,'Plot_flag'),  Plot_flag = i.Plot_flag;  
                     else   Plot_flag = 0;  end 

   
D       = 2*R;    % [m]
Dhub    = 2*Rhub; % [m]
Rhub_oR = Rhub/R;
N       = 60*Vs/(Js*D); % [RPM]
n       = N/60;         % [rev/s]
omega   = 2*pi*n;       % [rad/s]   angular rotation rate

Vs_dot = 0;
 n_dot = 0;
 

skew   = pchip(XR,skew0,RC);     % [deg], angular translation along mid-chord helix
rake   = pchip(XR,rake0,RC);     % [],   translation along propeller axis (3D X-axis) / diameter
 

% --------------------------------------------- Compute Expanded Area Ratio
EAR = (2*Z/pi)*trapz(  linspace(Rhub_oR,1,100)  ,  interp1(RC,CoD,  linspace(Rhub_oR,1,100)  ,'spline','extrap')  );  

% -------------------------------------------------------------------------
if isfield(i, 'Meanline'), Meanline  = i.Meanline;   else  Meanline  = 'NACA a=0.8 (modified)';       end
if isfield(i,'Thickness'), Thickness = i.Thickness;  else  Thickness = 'NACA 65A010 (Epps modified)'; end

% CLItilde    [ ]   ideal lift coefficient (should match with Meanline type)
% alphaItilde [deg] ideal angle of attack  (should match with Meanline type)
% If no values are given, assume 'NACA a=0.8' meanline
% -------------------------------------------------------------------------
% ---------------------- Find normalized 2D foil geometry (at x0 positions)
%   f0octilde = zeros(Mp+1,1);  % can either be scalars or given versus RG(1:Mp+1)
%    CLItilde = zeros(Mp+1,1);  % can either be scalars or given versus RG(1:Mp+1)
% alphaItilde = zeros(Mp+1,1);  % can either be scalars or given versus RG(1:Mp+1)
if iscell(Meanline)  % Assume Meanline is given versus XR
    
      f0octilde = zeros(length(XR),1);  % can either be scalars or given versus RG(1:Mp+1)
       CLItilde = zeros(length(XR),1);  % can either be scalars or given versus RG(1:Mp+1)
    alphaItilde = zeros(length(XR),1);  % can either be scalars or given versus RG(1:Mp+1)
    
    for i = 1:length(XR)
        [f0octilde(i), CLItilde(i), alphaItilde(i)] = GeometryFoil2D(Meanline{i},Thickness{i},0);
    end

      f0octilde = pchip(XR,   f0octilde, RC);
       CLItilde = pchip(XR,    CLItilde, RC);
    alphaItilde = pchip(XR, alphaItilde, RC);

else
    [f0octilde, CLItilde, alphaItilde] = GeometryFoil2D(Meanline,Thickness,0);
end
% -------------------------------------------------------------------------


alphaI = alphaItilde.*CL./CLItilde; % [deg]
f0oc   =   f0octilde.*CL./CLItilde;   % [ ],   camber ratio scaled with lift coefficient

BetaIC = atand(TANBIC); % [deg]

theta  = BetaIC+alphaI; % [deg]

PoD    = tand(theta).*pi.*RC;           % Pitch / propeller diameter, [ ]


theta  =  theta*pi/180; % [rad]
BetaIC = BetaIC*pi/180; % [rad]
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Blade geometry
[xl,xu,yl,yu, X3D,Y3D,Z3D] = Stress_BladeGeometry(D,RV,RC,CoD,t0oD,f0oc,skew,rake,PoD,Mp,Np,Meanline,Thickness);

% -------------------------------------------------------------------------
% Section Centroid and Moments of Inertia Calculation (functions of RC)
[Ixc, Iyc, Ixyc, Area, Xbar, Ybar, xl, yl, xu, yu] = Stress_Moment_of_Inertia(xl,xu,yl,yu);

% -------------------------------------------------------------------------
% Compute blade stress
[MAX_STRESS, MAX_n_dot, blade_stress] = Stress_Calculation(rho, rho_p, allowable_stress, Vs, Vs_dot, n, n_dot, R, Mp, Np, RC, DR, CD, CL, VSTAR, CoD, BetaIC, theta, ...
                        Ixc, Iyc, Ixyc, Area, Xbar, Ybar, xl, yl, xu, yu);

                    
                    
if Plot_flag == 0                    
                    
    [I J] = find(blade_stress == max(blade_stress(:)));

    % Make Stress Patch figure
    Plot_Blade_Contours2(X3D, Y3D, Z3D, blade_stress, 'Blade Stress (Pa)')
    
    plot3(X3D(I,J), Y3D(I,J), Z3D(I,J),'rs','MarkerFaceColor','r')

end


% MAX_STRESS
% 
% MAX_n_dot
%%


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

