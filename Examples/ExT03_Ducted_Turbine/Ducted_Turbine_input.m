% -------------------------------------------------------------------------
% Ducted turbine design example:
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
clear, close all, clc,

filename   = 'turbine';   % filename prefix
notes      = 'Tom Lokocz ducted turbine';           

% -------------------------------------------------------------------------
i.part1      = '------ Performance inputs ------';

i.Z         = 3;   		   % number of blades   
i.N         = 650;           % propeller speed [RPM]
i.Vs        = 1.25;             % free-stream speed [m/s]

i.D         = 0.254;          % rotor diameter [m] (Note: 39.37 in/m ) 
i.Dhub      = 0.25*i.D;       % hub diameter [m] (must be greater than 0.15*D) 

i.L         = pi*(i.N/60)*i.D/i.Vs;                      % tip-speed ratio

% -------------------------------------------------------------------------
input.part2      = '------ Geometry inputs ------';
i.Mp         = 20;            % number of vortex panels over the radius

i.XR         = [0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9    0.95   1.0];    % radius / propeller radius

% XCoD       = [0.1600 0.1818 0.2024 0.2196 0.2305 0.2311 0.2173 0.1806 0.1387 0.000001]; % chord / diameter unducted
i.XCoD       = [0.1600 0.1818 0.2024 0.2196 0.2305 0.2311 0.2173 0.19 0.17 0.15]; %(use this one) chord / diameter ducted
% XCoD       = [0.2600 0.2321 0.2109 0.1957 0.1900 0.1845 0.1800 0.1800 0.1800 0.1800]; %(old) chord / diameter ducted

i.XCD        = [0.0080 0.0080 0.0080 0.0080 0.0080 0.0080 0.0080 0.0080 0.0080 0.0080]; % section drag coefficient
i.XVA        = [1      1      1      1      1      1      1      1      1      1     ]; % axial      inflow velocity / ship velocity
i.XVT        = [0      0      0      0      0      0      0      0      0      0     ]; % tangential inflow velocity / ship velocity
i.t0oc0      = [0.2056 0.1551 0.1181 0.0902 0.0694 0.0541 0.0419 0.0332 0.0324 0.0000]; % max section thickness / chord
i.skew0      = [0      0      0      0      0      0      0      0      0      0     ]; % skew [deg]
i.rake0      = [0      0      0      0      0      0      0      0      0      0     ]; % rake / diameter

i.Meanline   = 'NACA a=0.8';           % Meanline type  (1 == NACA a=0.8, 2 == parabolic)
i.Thickness  = 'NACA 65A010';           % Thickness form (1 == NACA 65A010, 2 == elliptical, 3 == parabolic)

i.ALPHAstall = 8*pi/180;  % [rad], stall angle of attack - ideal angle of attack
i.dCLdALPHA  = 2*pi;      % d(CL)/d(alpha)

i.XCLmax = 1;

% -------------------------------------------------------------------------
i.part3      = '------ Computational inputs ------';

i.Propeller_flag  = 0;      % 0 == turbine, 1 == propeller
  i.Viscous_flag  = 0;      % 0 == viscous forces off (CD = 0), 1 == viscous forces on
      i.Hub_flag  = 0;      % 0 == no hub, 1 == hub
     i.Duct_flag  = 1;      % 0 == no duct, 1 == duct
     i.Plot_flag  = 1;      % 0 == do not display plots, 1 == display plots
    i.Chord_flag  = 1;      % 0 == do not optimize chord lengths, 1 == optimize chord lengths
    
i.ITER      = 50;            % number of iterations in wake alignment

    

% -------------------------------------------------------------------------
i.part4      = '------ Duct inputs ------';
i.Rduct      = i.D/2;          % duct radius [m]
i.Cduct      = i.D/2;          % duct chord length [m]
i.Xduct      = 0;              % duct axial displacement downstream [m]
i.CDd        = 0.008;        % duct viscous drag coefficient
i.CTD        = -0.4;          % duct thrust coefficient



% -------------------------------------------------------------------------
i.part5      = '------ Cavitation inputs ------';
i.rho       = 1000;          % water density [kg/m^3]
i.H         = 1;             % Shaft centerline depth [m]
i.dV        = 0.2;           % Inflow variation [m/s]





% =========================================================================       
% ---------------------------- Pack up propeller/turbine data structure, pt
pt.name     = filename; % (string) propeller/turbine name
pt.date     = date;     % (string) date created
pt.notes    = notes;    % (string or cell matrix)   notes
pt.i        = i;        % (struct) input parameters
pt.d        = [];       % (struct) design conditions
pt.g        = [];       % (struct) design geometry
pt.s        = [];       % (struct) off-design state analysis

% --------------------------------------------------------- Save input data
save OPinput pt 

clear, clc, 
pause(0.01),
pause(0.01),

load OPinput pt 

pause(0.01),
pause(0.01),

pt.i