% --------------------------------------------------------- Example_input.m
% Created: 3/12/2010, Brenden Epps, bepps@mit.edu
% 
% This script creates an "input." data structure for use in OpenProp.
%
% To design a propeller using these inputs, run:  design = EppsOptimizer(input)
%
% -------------------------------------------------------------------------
clear, close all, clc,

% ------------------------------------------------------------------------- 
i.filename   = 'MyTurbine';    % filename prefix for output files
i.date       = date;              % (string) date created
notes        = '';                  % design notes

% ------------------------------------------------------------------------- 
i.part1      = '------ Performance inputs ------';

i.L    = 5;              % tip speed ratio

% ------------------------------------------------------------------------- 
i.part2      = '------ Geometry inputs ------';

i.Z    = 3;       % number of blades     (Z == Inf for actuator disk)  
i.R    = 1;       % turbine radius
i.Rhub = 0.001;   % hub radius
i.Mp   = 20;      % vortex panels

% ----------------------------------------------- Chord length optimization    
% NOTE:  CD/CL == CDoCL = XCD/XCLmax
%
i.XCLmax = 1; % max allowable lift coefficient
i.XCD    = 0; % section drag coefficient:

% ------------------------------------------------------------------------- 
i.part3      = '------ Computational inputs ------';

i.Propeller_flag  = 0;    % 0 == turbine, 1 == propeller
      i.Hub_flag  = 0;    % 0 == no hub, 1 == hub
  i.Viscous_flag  = 1;    % 0 == viscous forces off (CD = 0), 1 == viscous forces on
    i.Chord_flag  = 1;    % 0 == do not optimize chord lengths, 1 == optimize chord lengths
    

  

% =========================================================================       
% ================================================= Pack up input variables
pt.filename = i.filename; % (string) propeller/turbine name
pt.date     = i.date;     % (string) date created
pt.notes    = notes;      % (string or cell matrix)   notes
pt.input    = i;          % (struct) input parameters
pt.design   = [];         % (struct) design conditions
pt.geometry = [];         % (struct) design geometry
pt.states   = [];         % (struct) off-design state analysis

% --------------------------------------------------------- Save input data
save OPinput pt 

clear, clc, pause(0.01),

load OPinput pt 

pause(0.01),

pt.input