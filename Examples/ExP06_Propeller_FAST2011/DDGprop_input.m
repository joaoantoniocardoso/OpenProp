% -------------------------------------------------------------------------
% Example case in (Epps et al., FAST 2011)
% -------------------------------------------------------------------------
clear, close all, clc,

filename   = 'DDGprop';   % filename prefix
notes      = '';           

% ------------------------------------------------------- Design parameters
Z         = 5;   		     % number of blades 
D         = 17.00*0.3048;    % propeller diameter [m] = [ft] * [0.3048 m/ft]
Dhub      = 0.2819*D;        % hub diameter [m]   
 

% thrust requirement with 9% thrust margin -- see ShipResistance100611.xls
Tmargin   = 1.09;          % thrust margin
Vs        = 20*0.51444;    % m/s == 20 knots * (0.5144 (m/s)/knot);         % ship velocity [m/s]   
THRUST    = 381529*Tmargin;% required thrust [N] = [lbf] * [4.448 N/lbf]
N         = 93.8;          % propeller speed [RPM] to give design advance coefficient of (Js = 1.27)
rho       = 1025;          % water density [kg/m^3] = [slug/ft^3] * (515.38 [kg/m^3]/[slug/ft^3])


% DRAFT = 20.7; % [ft]
H = (20.7)*0.3048;   % Shaft centerline depth [m] = [ft] * [0.3048 m/ft]


% "high speed" requirements with margin 
Vh      = 30*0.51444;            % [m/s] == 30 knots * (0.5144 (m/s)/knot);  
THRUSTh = 1266007*Tmargin;       % [N] 30 kt thrust requirement

  

Mp        = 20;            % number of vortex panels over the radius
Np        = 20;            % number of points along the chord
ITER      = 40;            % number of iterations in wake alignment


% --------------------------------------------- Blade 2D section properties
Meanline   = 'NACA a=0.8';                                  % Meanline type  
Thickness  = 'NACA 65A010';                                 % Thickness form


% USN propeller 5168
%   r/R         c/D     P/D     i/D  pitch [deg]   skew [deg] t0/c  f0/c
P5168 = [
    0.2819    0.18910  1.02029  0.01709  49.042  9.870  0.30984  -0.05404
    0.3000    0.20000  1.08750  0.01110  49.086  6.280  0.26936  -0.03800 
    0.3500    0.23700  1.24479 -0.00694  48.545 -0.754  0.17804  -0.00860
    0.4000    0.27600  1.36525 -0.02370  47.372 -4.824  0.11423   0.00903
    0.4500    0.30800  1.45791 -0.03579  45.882 -7.336  0.08897   0.01779
    0.5000    0.33410  1.54131 -0.04367  44.457 -8.865  0.07546   0.02789
    0.6000    0.38540  1.67347 -0.04825  41.599 -9.838  0.05858   0.03655
    0.7000    0.43500  1.63334 -0.04195  36.602 -8.108  0.04874   0.03323
    0.8000    0.47450  1.50246 -0.03025  30.871 -3.784  0.04108   0.02473
    0.9000    0.46500  1.33483 -0.01645  25.272  3.784  0.04376   0.01191
    0.9500    0.39000  1.18919 -0.00960  21.725  9.297  0.05883   0.00539
    0.9800    0.27696  1.03965 -0.00538  18.659 13.344  0.06402   0.00196
    1.0000    0.00000  0.90000 -0.00245  15.986 16.400  0.06222   0.00000];

XR     = P5168(:,1)';
XCoD   = P5168(:,2)';
XPoD   = P5168(:,3)'; % == pi*XR.*tand(theta0)
rake0  = P5168(:,4)';
pitch0 = P5168(:,5)';
skew0  = P5168(:,6)';
t0oc0  = P5168(:,7)';
f0oc0  = P5168(:,8)';
Xt0oD  = t0oc0 .* XCoD;

% Xt0oD = 0.02 + 0.08*(1-XR);


    XCD        = 0.008*ones(size(XR));   % CD3                       % section drag coefficient
    XVA        =  ones(size(XR));  
    XVT        = zeros(size(XR));                               
    
    
% ------------------------------------------------------------------- Flags
Propeller_flag  = 1;      % 0 == turbine, 1 == propeller
  Viscous_flag  = 1;      % 0 == viscous forces off (CD = 0), 1 == viscous forces on
      Hub_flag  = 1;      % 0 == no hub, 1 == hub
     Duct_flag  = 0;      % 0 == no duct, 1 == duct
     Plot_flag  = 1;      % 0 == do not display plots, 1 == display plots
    Chord_flag  = 1;      % 0 == do not optimize chord lengths, 1 == optimize chord lengths
 
ChordMethod = 'FAST2011dCTP';
    
  
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
input.THRUST     = THRUST;      % required thrust [N]

input.Vh         = Vh;          % [1 x 1], [m/s] ship speed -- high speed
input.THRUSTh    = THRUSTh;     % required thrust [N]       -- high speed

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
% input.t0oc0      = t0oc0;       % [length(XR) x 1], [ ] input thickness / chord      at XR 
input.Xt0oD      = Xt0oD;       % [length(XR) x 1], [ ] input thickness / D      at XR 
input.skew0      = skew0;       % [length(XR) x 1], [ ] input skew  [deg]            at XR 
input.rake0      = rake0;       % [length(XR) x 1], [ ] input rake X/D               at XR 
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
input.ChordMethod     = ChordMethod;
input.ITER            = ITER;           % [ ] number of iterations


input.part4      = '------ Cavitation inputs ------';
input.rho        = rho;         % [1 x 1], [kg/m^3] fluid density
input.H          = H;           % [1 x 1]


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

clear, clc, load OPinput pt input

pt.input