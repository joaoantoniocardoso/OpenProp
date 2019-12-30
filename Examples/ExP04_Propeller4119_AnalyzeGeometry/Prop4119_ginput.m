% -------------------------------------------------------------------------
% DTMB propeller 4119 geometry.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
clear, close all, clc,

filename   = 'DTMB4119';   
notes      = 'inputs for geometry analysis'; 

% -------------------------------------------------------------------------
i.part1     = '------ Performance inputs ------';

i.Z         = 3;   	     % number of blades 
i.Js        = 0.833;     % advance coefficient

i.D         = 1.00;      % propeller diameter [m] = [ft] * [0.3048 m/ft]
i.Dhub      = 0.2 ;      % hub diameter [m] = [ft] * [0.3048 m/ft]


% -------------------------------------------------------------------------
i.part2      = '------ Geometry inputs ------';

i.Meanline    = 'NACA a=0.8';              % Meanline type
i.Thickness   = 'NACA66 (DTRC Modified)';  % Thickness form 


i.XR         = [0.2000  0.3000  0.4000  0.5000  0.6000  0.7000  0.8000  0.9000  0.9500  1.0000];  % radius / propeller radius
i.XCoD       = [0.3200  0.3625  0.4048  0.4392  0.4610  0.4622  0.4347  0.3613  0.2775  0.0020];  % chord / diameter  --- NOTE FINITE CHORD AT TIP FOR NUMERICAL STABILITY ONLY AND NOT IN ORIGINAL PROP 4119 DESIGN
i.XCD        = [0.0080  0.0080  0.0080  0.0080  0.0080  0.0080  0.0080  0.0080  0.0080  0.0080];  % section drag coefficient
i.XVA        = [1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000  1.0000];  % axial      inflow velocity / ship velocity
i.XVT        = [0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000];  % tangential inflow velocity / ship velocity
i.f0oc0      = [0.01429 0.02318 0.02303 0.02182 0.02072 0.02003 0.01967 0.01817 0.01631 0.01175]; % max section camber    / chord
i.XPoD       = [1.105   1.102   1.098   1.093   1.088   1.084   1.081   1.079   1.077   1.075  ]; % pitch / diameter
i.t0oc0      = [0.2055  0.1553  0.1180  0.09016 0.06960 0.05418 0.04206 0.03321 0.03228 0.03160]; % max section thickness / chord
i.skew0      = zeros(size(i.XR));   % skew
i.rake0      = zeros(size(i.XR));   % rake

i.Xt0oD      = i.t0oc0 .* i.XCoD;   % thickness/diameter

        
% ------------------------------------------------------------------- Flags
i.part3      = '------ Computational inputs ------';

i.Propeller_flag  = 1;      % 0 == turbine, 1 == propeller
  i.Viscous_flag  = 1;      % 0 == viscous forces off (CD = 0), 1 == viscous forces on
      i.Hub_flag  = 1;      % 0 == no hub, 1 == hub
     i.Duct_flag  = 0;      % 0 == no duct, 1 == duct
     i.Plot_flag  = 0;      % 0 == do not display plots, 1 == display plots
    i.Chord_flag  = 0;      % 0 == do not optimize chord lengths, 1 == optimize chord lengths


% =========================================================================       
% ---------------------------- Pack up propeller/turbine data structure, pt
pt.filename = filename; % (string) propeller/turbine name
pt.date     = date;     % (string) date created
pt.notes    = notes;    % (string or cell matrix)   notes
pt.ginput   = i;        % (struct) input parameters
pt.states   = [];       % (struct) off-design state analysis

% --------------------------------------------------------- Save input data
save OPginput pt 

clear, clc, 
pause(0.01),
pause(0.01),

load OPginput pt 

pause(0.01),
pause(0.01),

pt.ginput