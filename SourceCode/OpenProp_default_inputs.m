% -------------------------------------------------------------------------
% OpenProp default inputs for EppsOptimizer, Geometry, and cavitation analysis
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
clear, close all, clc,

filename   = 'OpenProp';   
notes      = '999 means required input'; 


% -------------------------------------------------------------------------
i.part1    = '------ Performance inputs ------';
i.Js       = 999;  % advance coefficient (for turbine,   Js = pi/L )
i.L        = 999;  % tip speed ratio     (for propeller, L  = pi/Js)

i.CT       = 999;  % desired total thrust coefficient (propeller plus duct)

%     % CTdes == desired total thrust coefficient == CTPdes + CTDdes
%     if     isfield(input,'THRUST'),
%             CTdes = input.THRUST / (0.5*rho*Vs^2*pi*R^2);
%              
%     elseif isfield(input,'CTDES'), 
%             CTdes = input.CTDES; 
%             
%     elseif isfield(input,'CT'), 
%             CTdes = input.CT; 
%             
%     elseif isfield(input,'KTDES'),
%             CTdes = input.KTDES * (8/pi)/Js^2;
%             
%     elseif isfield(input,'KT'),
%             CTdes = input.KT    * (8/pi)/Js^2;
%     else 
%             CTdes = 0;




% -------------------------------------------------------------------------
i.part2      = '------ Geometry inputs ------';

i.Z     = 999;           % number of blades



i.Meanline  = 'NACA a=0.8 (modified)';         
i.Thickness = 'NACA 65A010'; 

% -------------------------------------------------------------------------
i.part3      = '------ Blade/Inflow inputs ------';

i.XR    = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];    % r/R
X1      =  ones(size(i.XR));
X0      = zeros(size(i.XR));
i.XCoD  = [0.1600 0.1818 0.2024 0.2196 0.2305 0.2311 0.2173 0.1806 0.1387 0.0100];          % c/D
i.t0oc0 = [0.2056 0.1551 0.1181 0.0902 0.0694 0.0541 0.0419 0.0332 0.0324 0.0100];          % t0/c  
i.Xt0oD = i.t0oc0 .* i.XCoD;              
i.XCD   = 0.008*X1;  
i.XVA   = X1;                                      % Va/Vs at XR
i.XVT   = X0;                                      % Vt/Vs at XR
i.dXVA  = X0;
i.XCLmax      = 0.5 + (0.2-0.5)/(1-i.XR(1)) * (i.XR-i.XR(1));  % XCLmax == 0.5 at root and 0.2 at tip



% Inflow velocity profiles (optional)
i.ri   = i.XR;                % m, radii of inflow profile (i.e. i.XR * i.R)
i.VAI  =  ones(size(i.ri));   % Va/Vs  at ri
i.VTI  = zeros(size(i.ri));   % Vt/Vs  at ri
i.dVAI = zeros(size(i.ri));   % delta(Va)/Vs at ri      
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
i.part4 = '------ Computational inputs ------';
i.Propeller_flag = 999;  % 0 == turbine, 1 == propeller
i.Viscous_flag   = 999;  % 0 == viscous forces off (CD = 0), 1 == viscous forces on
i.Hub_flag       = 0;    % 0 == no hub, 1 == hub
i.Duct_flag      = 0;    % 0 == no duct, 1 == duct
i.Chord_flag     = 0;    % 0 == do not optimize chord lengths, 1 == optimize chord lengths
i.Plot_flag      = 0;    % 0 == do not display plots, 1 == display plots

i.Make2Dplot_flag = 1;     
i.Make3Dplot_flag = 1;    
i.Make_Rhino_flag = 0;     
i.Make_SWrks_flag = 0;  

i.QuarterChord_flag  = 0;  % 0 == lifting line at mid-chord, 1 == lifting line at quarter-chord (i.e. skew offset == 0.25*c/r)
i.LSGeoCorr = 'none';      % propeller Lifting Surface Geometry Corrections



i.ChordMethod = 'CLmax';     % Method of chord optimization {CLmax, ConeyPLL, FAST2011dCTP, FAST2011dVAC}

i.EARspec = 0;    % Specified Expanded Area Ratio (EAR)

% % For FAST2011dCTP case: 
% i.Vh  = i.Vs;   % m/s
% i.CTh = 999; 

i.Mp    = 20;   % number of panels radially
i.Np    = 20;   % number of panels chordwise
i.ITER  = 50;   % number of iterations in solvers 
i.HUF   = 0;    % hub unloading factor
i.TUF   = 0;    % tip unloading factor
i.Rhv   = 0.5;  % hub vortex radius / hub radius
i.ALPHAstall  = 8*pi/180; % [rad] stall angle of attack
% i.dCLdALPHA = 0; % Lift curve slope: (2*pi / (1 + 2/PAR)) * ones(size(RC));

% -------------------------------------------------------------------------

% ------------------------------------------------------- Cavitation inputs
i.part5 = '------ Dimensional inputs (SI units) ------';
i.Vs   = 1;        % m/s
i.R    = 1;        % m, rotor radius
i.Rhub = 0.2*i.R;  % m, hub radius
i.Dm   =   2*i.R;  % m, model diameter (for geometry outputs only)
i.rho  = 1000;     % kg/m^3
i.H    = 3.048;    % m
i.g    = 9.81;     % m/s^2
i.Patm = 101325;   % Pa
i.Pv   = 2500;     % Pa
% -------------------------------------------------------------------------


% --------------------------------------------------------------- Duct_flag
i.part6 = '------ Duct inputs ------';
i.TAU      = 1;         % thrust ratio == propeller thrust / total thrust 
i.Rduct_oR = 1;         % duct radius
i.Cduct_oR = 1;         % duct chord length
i.Xduct_oR = 0;         % duct axial position downstream
i.CDd      = 0.008;     % duct drag coefficient


%     if Propeller_flag == 1
% 
%     i.TAU      = 1;         % thrust ratio == propeller thrust / total thrust 
% 
%     %  Allocate CTdes between propeller and duct
%         CTPdes = CTdes*   TAU;                       % CT desired for the propeller
%         CTDdes = CTdes*(1-TAU);                      % CT desired for the duct   
% 
%     else
% 
%         % CTDdes == desired duct thrust coefficient
%         if     isfield(input,'THRUSTduct'),
%                CTDdes = input.THRUSTduct / (0.5*rho*Vs^2*pi*R^2); 
% 
%         elseif isfield(input,'CTD'), 
%                CTDdes = input.CTD;       
% 
%         elseif isfield(input,'KTD'),
%                CTDdes = input.KTD   * (8/pi)/Js^2;
%         else 
%                CTDdes = 0;
% 
% 
%         TAU    = 0;  % (not used) assume given CTdes is for the duct alone, not turbine   + duct
%         CTPdes = 0;  % (not used)
%         end
%     end
  
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Derived quantities:
%     D        = 2*R;    % [m]
%     Dhub     = 2*Rhub; % [m]
%     N        = 60*Vs/(Js*D); % [RPM]

% =========================================================================       
% ---------------------------- Pack up propeller/turbine data structure, pt
pt.filename = filename; % (string) propeller/turbine name
pt.date     = date;     % (string) date created
pt.notes    = notes;    % (string or cell matrix)   notes
pt.input    = i;        % (struct) input parameters
pt.design   = [];       % (struct) design conditions
pt.geometry = [];       % (struct) design geometry
pt.states   = [];       % (struct) off-design state analysis
  
