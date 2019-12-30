% --------------------------------------------------------- Example_input.m
% Created: 5/28/09, Brenden Epps, bepps@mit.edu
% 
% This script creates an "input." data structure for use in OpenProp.
%
% To design a propeller using these inputs, run:  design = EppsOptimizer(input)
%
% -------------------------------------------------------------------------
clear, close all, clc,

filename   = 'PropStudy';   % filename prefix
notes      = '';            % design notes

% ------------------------------------------------------- Design parameters
Z         = 5;   		   % number of blades   
D         = 1;          % (approx 10 in) propeller diameter [m] (Note: 39.37 in/m ) 
Dhub      = 0.2;       % (3.3 in) hub diameter [m] (must be greater than 0.15*D) 

rho       = 1000;          % water density [kg/m^3]
R         = D/2;
Vs        = 1;             % ship velocity [m/s]
THRUST    = 0.512 * 0.5*rho*Vs^2*pi*R^2;            % (11.240 lb) required thrust [N] (0.2248 lb/N)


Js        = 0.1;            % advance ratio == Vs/nD

n         = Vs/(Js*D);      % [rev/s]
N         = 60*n;           % propeller speed [RPM]


Mp        = 40;            % number of vortex panels over the radius
Np        = 20;            % Number of points over the chord for geometry plots [ ]
ITER      = 50;            % number of iterations



% --------------------------------------------- Blade 2D section properties
Meanline   = 'NACA a=0.8 (modified)';           % Meanline type  (1 == NACA a=0.8, 2 == parabolic)
Thickness  = 'NACA 65A010 (modified)';           % Thickness form (1 == NACA 65A010, 2 == elliptical, 3 == parabolic)

XR         = [0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9    0.95   0.9750 0.9875 1.0];     % radius / propeller radius
XCoD       = [0.2002 0.2274 0.2531 0.2743 0.2878 0.2913 0.2817 0.2410 0.1962 0.1480 0.0988 0.0259];  % chord / diameter (factor of 1.25 was chosen to make t0oc look good)
XCD        = [1      1      1      1      1      1      1      1      1      1      1      1]*0.008; % section drag coefficient
XVA        = [1      1      1      1      1      1      1      1      1      1      1      1]; % axial      inflow velocity / ship velocity
XVT        = [0      0      0      0      0      0      0      0      0      0      0      0]; % tangential inflow velocity / ship velocity
skew0      = [0      0      0      0      0      0      0      0      0      0      0      0]; % skew [deg]
rake0      = [0      0      0      0      0      0      0      0      0      0      0      0]; % rake / diameter


% Thickness profile (see notes below)
TTRF  = 0.5;  % Tip Thickness Reduction Factor == modified thickness at tip / baseline thickness at tip
XRmax = 0.8;  % maximum XR for which thickness reduction is less than 1%
HTTR  = 3.5;                       % Hub-Tip Thickness Ratio = t0(hub) / t0(tip)
t0tip = 0.00254;               % [m] == 0.254 mm = 0.1 inch, max thickness at tip section
t0oc0 = t0tip*(HTTR - (HTTR-1).*(XR-Dhub/D)/(1-Dhub/D))./(XCoD*D) .* (1-(1-TTRF)*exp(-4.6*(1-XR)/(1-XRmax)));

Xt0oD = t0oc0 .* XCoD;

% ------------------------------------------------------------------- Flags
Propeller_flag  = 1;      % 0 == turbine, 1 == propeller
  Viscous_flag  = 1;      % 0 == viscous forces off (CD = 0), 1 == viscous forces on
      Hub_flag  = 0;      % 0 == no hub, 1 == hub
     Duct_flag  = 0;      % 0 == no duct, 1 == duct
    Chord_flag  = 0;    
     Plot_flag  = 0;      % 0 == do not display plots, 1 == display plots
 

% ---------------------------------------------- Compute derived quantities
R       = D/2;                        % propeller radius [m]
Rhub    = Dhub/2;                     % hub radius [m]
Rhub_oR = Rhub/R;
L       = pi/Js;                        % tip-speed ratio
CTDES   = THRUST/(0.5*rho*Vs^2*pi*R^2); % CT thrust coefficient required          
    




% =========================================================================       
% ================================================= Pack up input variables
input.part1      = '------ Performance inputs ------';
input.Z          = Z;           % [1 x 1], [ ] number of blades
input.N          = N;           % propeller speed [RPM]
input.D          = D;           % propeller diameter [m]  
input.Vs         = Vs;          % [1 x 1], [m/s] ship speed
input.Js         = Js;          % [1 x 1], [ ] advance coefficient, Js = Vs/nD = pi/L
input.L          = L;           % [1 x 1], [ ] tip speed ratio, L = omega*R/V
input.CT         = CTDES;       % [1 x 1], [ ] desired thrust coefficient

input.part2      = '------ Geometry inputs ------';
input.Mp         = Mp;          % [1 x 1], [ ] number of blade sections
input.Np         = Np;          % [1 x 1], [ ] number of points along the chord
input.R          = R;           % [1 x 1], [m] propeller radius
input.Rhub       = Rhub;        % [1 x 1], [m] hub radius
input.XR         = XR;          % [length(XR) x 1], [ ] input radius/propeller radius
input.XVA        = XVA;         % [length(XR) x 1], [ ] input axial inflow velocity  at XR
input.XVT        = XVT;         % [length(XR) x 1], [ ] input swirl inflow velocity  at XR
input.XCD        = XCD;         % [length(XR) x 1], [ ] input drag coefficient       at XR
input.XCoD       = XCoD;        % [length(XR) x 1], [ ] input chord / diameter       at XR
input.Xt0oD      = Xt0oD;       % [length(XR) x 1], [ ] input thickness / diameter   at XR 
input.skew0      = skew0;       % [length(XR) x 1], [ ] input skew  [deg]      at XR 
input.rake0      = rake0;       % [length(XR) x 1], [ ] input rake X/D       at XR 
input.Meanline   = Meanline;    % 2D section meanline  flag
input.Thickness  = Thickness;   % 2D section thickness flag 


input.part3      = '------ Computational inputs ------';
input.Propeller_flag  = Propeller_flag; % 0 == turbine, 1 == propeller
input.Viscous_flag    = Viscous_flag;   % 0 == viscous forces off (CD = 0), 1 == viscous forces on
input.Hub_flag        = Hub_flag;       % 0 == no hub, 1 == hub
input.Duct_flag       = Duct_flag;      % 0 == no duct, 1 == duct
input.Plot_flag       = Plot_flag;      % 0 == do not display plots, 1 == display plots
input.Chord_flag      = Chord_flag;     % 0 == do not optimize chord lengths, 1 == optimize chord lengths
input.ITER            = ITER;           % [ ] number of iterations


% ---------------------------- Pack up propeller/turbine data structure, pt
pt.filename = filename; % (string) propeller/turbine name
pt.date     = date;     % (string) date created
pt.notes    = notes;    % (string or cell matrix)   notes
pt.i        = input;    % (struct) input parameters
pt.d        = [];       % (struct) design conditions
pt.g        = [];       % (struct) design geometry
pt.s        = [];       % (struct) off-design state analysis

% --------------------------------------------------------- Save input data
save OPinput pt 

clear, clc, load OPinput pt

pt.i