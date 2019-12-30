% --------------------------------------------------------- Example_input.m
% Created: 3/12/2010, Brenden Epps, bepps@mit.edu
% 
% This script creates an "input." data structure for use in OpenProp.
%
% To design a propeller using these inputs, run:  design = EppsOptimizer(input)
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
clear, close all, clc,

% ------------------------------------------------------------------------- 
filename   = 'ActuatorDisk';    % filename prefix for output files
notes        = '';                % design notes

% ------------------------------------------------------------------------- 
i.part1      = '------ Performance inputs ------';

i.Js    = 0.001;           % advance ratio        (Js == 0 for actuator disk)
i.CT    = 0.01;            % thrust coefficient   (may be anything)


% ------------------------------------------------------------------------- 
i.part2      = '------ Geometry inputs ------';

i.Z     = 100;             % number of blades     (Z == Inf for actuator disk)  
i.D     = 1;      
i.Dhub  = 0.001;

% ------------------------------------------------------------------------- 
i.part3      = '------ Computational inputs ------';

i.Propeller_flag  = 1;    % 0 == turbine, 1 == propeller
      i.Hub_flag  = 0;    % 0 == no hub, 1 == hub
  i.Viscous_flag  = 0;    % 0 == viscous forces off (CD = 0), 1 == viscous forces on

  

% =========================================================================       
% ================================================= Pack up input variables
pt.filename = filename;   % (string) propeller/turbine name
pt.date     = date;       % (string) date created
pt.notes    = notes;      % (string or cell matrix)   notes
pt.input    = i;          % (struct) input parameters
pt.design   = [];         % (struct) design conditions
pt.geometry = [];         % (struct) design geometry
pt.states   = [];         % (struct) off-design state analysis

% --------------------------------------------------------- Save input data
save OPinput pt 

clr, pause(0.01),

load OPinput pt 

pause(0.01),

pt.input
