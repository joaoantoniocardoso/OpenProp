% -------------------------------------------------------------------------
clear, close all, clc,

filename   = 'OpenProp';   % filename prefix
notes      = 'Ducted propeller from Sutbblefield (2008) M.S. thesis';           

% ------------------------------------------------------- Design parameters
Z         = 5;   		   % number of blades   
N         = 150;           % propeller speed [RPM]
D         = 3.048;         % propeller diameter [m]
   
THRUST    = 94328;         % required thrust [N]
Vs        = 4.572;         % ship velocity [m/s]
Dhub      = 0.6096;        % hub diameter [m] (must be greater than 0.15*D) 

Mp        = 20;            % number of vortex panels over the radius
Np        = 20;            % number of points along the chord

rho       = 1031;          % sea water density [kg/m^3]

% --------------------------------------------------------- Duct parameters
% Inputs for no duct: Duct_flag = 0; TAU = 1; Rduct_oR = 1; CDd = 0;
TAU        = 0.8;          % thrust ratio == prop thrust / total thrust
Rduct      = D/2;          % duct radius [m]
Cduct      = D/2;          % duct chord length [m]
CDd        = 0.008;        % duct viscous drag coefficient


% --------------------------------------------- Blade 2D section properties
Meanline   = 'NACA a=0.8';           % Meanline type  (1 == NACA a=0.8, 2 == parabolic)
Thickness  = 'NACA 65A010';           % Thickness form (1 == NACA 65A010, 2 == elliptical, 3 == parabolic)

XR         = [0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9    0.95   1.0];    % radius / propeller radius
XCoD       = [0.1600 0.1818 0.2024 0.2196 0.2305 0.2311 0.2173 0.1806 0.1387 0.0010]; % chord / diameter
XCD        = [0.0080 0.0080 0.0080 0.0080 0.0080 0.0080 0.0080 0.0080 0.0080 0.0080]; % section drag coefficient
XVA        = [1      1      1      1      1      1      1      1      1      1     ]; % axial      inflow velocity / ship velocity
XVT        = [0      0      0      0      0      0      0      0      0      0     ]; % tangential inflow velocity / ship velocity
t0oc0      = [0.2056 0.1551 0.1181 0.0902 0.0694 0.0541 0.0419 0.0332 0.0324 0.0000]; % max section thickness / chord
skew0      = [0      0      0      0      0      0      0      0      0      0     ]; % skew [deg]
rake0      = [0      0      0      0      0      0      0      0      0      0     ]; % rake / diameter


% ------------------------------------------------------------------- Flags
Propeller_flag  = 1;      % 0 == turbine, 1 == propeller
  Viscous_flag  = 1;      % 0 == viscous forces off (CD = 0), 1 == viscous forces on
      Hub_flag  = 1;      % 0 == no hub, 1 == hub
     Duct_flag  = 1;      % 0 == no duct, 1 == duct
     Plot_flag  = 1;      % 0 == do not display plots, 1 == display plots
    Chord_flag  = 0;      % 0 == do not optimize chord lengths, 1 == optimize chord lengths

% ---------------------------------------------- Compute derived quantities
n       = N/60;                       % revolutions per second [rps]
R       = D/2;                        % propeller radius [m]
Rhub    = Dhub/2;                     % hub radius [m]
Rhub_oR = Rhub/R;
Js      = Vs/(n*D);                   % advance coefficient
L       = pi/Js;                      % tip-speed ratio
CTDES   = THRUST/(0.5*rho*Vs^2*pi*R^2); % CT thrust coefficient required          
 
dCLdALPHA  = 2*pi;      % d(CL)/d(alpha)

% =========================================================================       
% ================================================= Pack up input variables
input.part1      = '------ Performance inputs ------';
input.Z          = Z;           % [1 x 1], [ ] number of blades
input.N          = N;           % propeller speed [RPM]
input.D          = D;           % propeller diameter [m]  
input.Vs         = Vs;          % [1 x 1], [m/s] ship speed
input.Js         = Js;          % [1 x 1], [ ] advance coefficient, Js = Vs/nD = pi/L
input.L          = L;           % [1 x 1], [ ] tip speed ratio, L = omega*R/V
input.CTDES      = CTDES;       % [1 x 1], [ ] desired thrust coefficient

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
input.t0oc0      = t0oc0;       % [length(XR) x 1], [ ] input thickness / chord      at XR 
input.skew0      = skew0;       % [length(XR) x 1], [ ] input skew  [deg]      at XR 
input.rake0      = rake0;       % [length(XR) x 1], [ ] input rake X/D       at XR 
input.Meanline   = Meanline;    % 2D section meanline  flag
input.Thickness  = Thickness;   % 2D section thickness flag 
input.dCLdALPHA  = dCLdALPHA;   % d(CL)/d(alpha)

input.part3      = '------ Computational inputs ------';
input.Propeller_flag  = Propeller_flag; % 0 == turbine, 1 == propeller
input.Viscous_flag    = Viscous_flag;   % 0 == viscous forces off (CD = 0), 1 == viscous forces on
input.Hub_flag        = Hub_flag;       % 0 == no hub, 1 == hub
input.Duct_flag       = Duct_flag;      % 0 == no duct, 1 == duct
input.Plot_flag       = Plot_flag;      % 0 == do not display plots, 1 == display plots
input.Chord_flag      = Chord_flag;     % 0 == do not optimize chord lengths, 1 == optimize chord lengths


input.part4      = '------ Duct inputs ------';
input.TAU        = TAU;         % [1 x 1], [ ] propeller thrust / total thrust
input.Rduct      = Rduct;       % [1 x 1], [m] duct radius
input.Cduct      = Cduct;       % [1 x 1], [m] duct chord length
input.CDd        = CDd;         % [1 x 1], [ ] duct drag coefficient


% ---------------------------- Pack up propeller/turbine data structure, pt
pt.filename = filename; % (string) propeller/turbine name
pt.date     = date;     % (string) date created
pt.notes    = notes;    % (string or cell matrix)   notes
pt.input    = input;    % (struct) input parameters
pt.design   = [];       % (struct) design conditions
pt.geometry = [];       % (struct) design geometry
pt.states   = [];       % (struct) off-design state analysis

% --------------------------------------------------------- Save input data
save OPinput pt input

clear, clc, 
pause(0.01),
pause(0.01),

load OPinput, 

pause(0.01),
pause(0.01),

pt.input