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
% Last modified: 11/15/2011 Brenden Epps
% =========================================================================

% =========================================================================
% =============================================================== Analyze.m
%
% This function computes the "states" of a given propeller/turbine "design"
% at given off-design tip speed ratios, "LAMBDAall".
%
% -------------------------------------------------------------------------
%
% Inputs:
%       pt          propeller/turbine data structure with fields:
%                       pt.i    or     pt.input
%                       pt.d    or     pt.design
% 
%       LAMBDAall   tip-speed ratios, omega*R/Vs
% 
%       feather     [deg], feathering angle, i.e. pitch angle displacement
%                          relative to design geometry, such that pitch 
%                          angle with feathering == pitch angle + feather
%                          (for controllable-pitch rotors)
% Outputs:   
%       states      data structure with operating states at LAMBDAall
%
% -------------------------------------------------------------------------
% Method:
%       This is a continuation of OpenProp v3.1.0 Analyze110718.m
%
%       Newton solver state variable: (VSTAR,ALPHA,CL,G), size [4*Mp,1] 
%       with (UASTAR,UTSTAR,TANBIC,...) updated between iterations
%
%       dCLdALPHA is estimated from elliptic wing unless specified in
%                   pt.input.dCLdALPHA  or  pt.design.dCLdALPHA
%
% -------------------------------------------------------------------------
            
function [states] = Analyze(pt,LAMBDAall,feather)

if nargin == 2
    feather = 0;
end

% =========================================================================
feather = feather * pi/180;  % [rad] feathering angle
% =========================================================================


warning off
if isfield(pt,'i'), pt.input  = pt.i; end
if isfield(pt,'d'), pt.design = pt.d; end

design = pt.design;
input  = pt.input;

% ========================================================== Error checking
% 0 == Horseshoe(...,Wrench(...)), 1 == Wake_Horseshoe(...)      
if isfield(input,'Wake_flag'),  Wake_flag = input.Wake_flag;  
                         else   Wake_flag = 0;  end
                         
if Wake_flag == 1, disp('ERROR: Wake_Horseshoe(...) is no longer supported in AnalyzeGeometry110628.'), return, end

% =========================================================================



% ------------------------------------------------------------------------- 
if isfield(input,'Duct_flag'),  Duct_flag = input.Duct_flag;  
                         else   Duct_flag = 0;  end

if Duct_flag == 1 
    Rduct_oR = design.Rduct_oR;
    Cduct_oR = design.Cduct_oR;
    Xduct_oR = design.Xduct_oR;    
    RC       = design.RC;
    RV       = design.RV;
    Z        =  input.Z;
    Hub_flag =  input.Hub_flag;
    Rhub_oR  = design.Rhub_oR;
    
    % Find propeller influence at the duct quarter chord, and
    % put these in pt.design, i.e. where they can be accessed from Analyze_some(...)
    disp(' '), disp('Computing rotor-duct interaction...be patient...'), disp(' '),
    [pt.design.DAHIFq_times_TANBIC,junk,pt.design.DRHIFq_times_TANBIC] = Horseshoe_intr_110830(Xduct_oR-Cduct_oR/4,Rduct_oR ,RC,ones(size(RC)),RV,Z,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR); 
end
% ------------------------------------------------------------------------- 


%%
% =========================================================================
% -------------------------- Initialize output data structure "states.____"
% Sort LAMBDAall into asending order
LAMBDAall = sort(LAMBDAall);

NL = length(LAMBDAall);
Mp = length(pt.design.RC);
    
states.part1   = '------ Off-design states ------';
states.L       = LAMBDAall(:);    % [NL x 1] tip speed ratio 
states.Js      = pi./states.L;    % [NL x 1] advance coefficient 

states.RC      = pt.design.RC;    % [ 1 x Mp] radius / propeller radius
states.CoD     = pt.design.CoD;   % [ 1 x Mp] chord  / propeller diameter 

states.G       = zeros(NL,Mp);    % [NL x Mp] circulation distribution
states.UASTAR  = zeros(NL,Mp);    % [NL x Mp] axial      induced velocity distribution
states.UTSTAR  = zeros(NL,Mp);    % [NL x Mp] tangential induced velocity distribution
if Wake_flag == 1
states.URSTAR  = zeros(NL,Mp);    % [NL x Mp] radial     induced velocity distribution
end
states.VSTAR   = zeros(NL,Mp);    % [NL x Mp] total inflow       velocity distribution
states.TANBC   = zeros(NL,Mp);    % [NL x Mp] tangent of free-stream inflow angle
states.TANBIC  = zeros(NL,Mp);    % [NL x Mp] tangent of total inflow angle

states.ALPHA       = zeros(NL,Mp);    % [NL x Mp] ALPHA = alpha - alphaI
% states.deltaALPHA1 = zeros(NL,Mp);    % [NL x Mp] 
% states.deltaALPHA2 = zeros(NL,Mp);    % [NL x Mp]
states.CL          = zeros(NL,Mp);    % [NL x Mp] lift coefficient distribution
states.CD          = zeros(NL,Mp);    % [NL x Mp] drag coefficient distribution



states.part2  = '------ Numerical metrics ------';
states.converged  = zeros(NL,1);
states.iteration  = zeros(NL,1);


if Duct_flag == 1  % Duct is present 
    states.part3   = '------ Duct parameters ------';

    Nd = length(design.Xduct_oR);
    states.Gd      = zeros(NL,1);     % [NL x 1]  duct circulation
    states.UARINGq = zeros(NL,1);     % [NL x 1]  UA    at duct quarter chord
    states.URRINGq = zeros(NL,1);     % [NL x 1]  UR    at duct quarter chord 
    states.VSRINGq = zeros(NL,1);     % [NL x 1]  VSTAR at duct quarter chord
    states.BetaIDq = zeros(NL,1);     % [NL x 1]  angle at duct quarter chord 
    states.CLd     = zeros(NL,1);     % [NL x 1]  duct lift coefficient
    states.UARING  = zeros(NL,Nd);    % [NL x Nd], induced velocity at duct
    states.URRING  = zeros(NL,Nd);    % [NL x Nd], induced velocity at duct
    states.UADUCT  = zeros(NL,Mp);    % [NL x Mp], induced velocity at propeller
    states.TAU     = zeros(NL,1);     % [NL x 1]  thrust ratio 
    states.CTD     = zeros(NL,1);     % [NL x 1]  CT for duct

    states.part4  = '------ Performance metrics ------';
else
    states.part3  = '------ Performance metrics ------';

end

states.KT      = zeros(NL,1);     % [NL x 1]  thrust coefficient
states.KQ      = zeros(NL,1);     % [NL x 1]  torque coefficient
states.CT      = zeros(NL,1);     % [NL x 1]  thrust coefficient
states.CQ      = zeros(NL,1);     % [NL x 1]  torque coefficient
states.CP      = zeros(NL,1);     % [NL x 1]  power  coefficient

if abs(design.VMIV - 1) > 1e-8,  % i.e. if VMIV is not equal to 1 
    states.EFFYo  = zeros(NL,1); % open water efficiency
    states.Ja     = zeros(NL,1); % == VMIV*Js
end
states.EFFY    = zeros(NL,1);     % [NL x 1]  efficiency == VMIV*EFFYo

states.VMWV    = zeros(NL,1);               % [1 x 1]

% -------------------------------------------------------------------------

%%              
% ---------------------- Analyze lower LAMBDA states first, then higher values
if isfield(input,'L')
    LAMBDA0 =    input.L;
else
    LAMBDA0 = pi/input.Js;
end
LAMBDA_low  = LAMBDAall(find(LAMBDAall <= LAMBDA0,1,'last'):-1:1  );
LAMBDA_high = LAMBDAall(find(LAMBDAall <= LAMBDA0,1,'last'): 1:end);


% If all LAMBDA are greater than the design LAMBDA, then 
% analyze all starting from the lowest LAMBDA ...
if isempty(LAMBDA_low)
    LAMBDA_low  = [];
    LAMBDA_high = [LAMBDAall(1);LAMBDAall(:)];
end

% % ----- 12/8/08 TEST CODE: (to anlayze all starting from the highest L...)
% LAMBDA_low  = LAMBDAall(end:-1:1);
% LAMBDA_high = [];
% % ----- END TEST CODE

% % ----- 12/8/08 TEST CODE: (to anlayze all starting from the lowest L...)
% LAMBDA_low  = [];
% LAMBDA_high = [LAMBDAall(1);LAMBDAall(:)];
% % ----- END TEST CODE

% ------------------- Analyze lower LAMBDA states first, then higher values
    [s_low]  = Analyze_some(pt,LAMBDA_low ,feather);
    [s_high] = Analyze_some(pt,LAMBDA_high,feather);
     
   
% ------------------------------- Output results in increasing LAMBDA order
states.G         = [     s_low.G(end:-1:1,:);     s_high.G(2:end,:)];   % [NL x Mp]
states.UASTAR    = [s_low.UASTAR(end:-1:1,:);s_high.UASTAR(2:end,:)];   % [NL x Mp] 
states.UTSTAR    = [s_low.UTSTAR(end:-1:1,:);s_high.UTSTAR(2:end,:)];   % [NL x Mp] 
if Wake_flag == 1
states.URSTAR    = [s_low.URSTAR(end:-1:1,:);s_high.URSTAR(2:end,:)];   % [NL x Mp] 
end
states.VSTAR     = [ s_low.VSTAR(end:-1:1,:); s_high.VSTAR(2:end,:)];   % [NL x Mp]     
states.TANBC     = [ s_low.TANBC(end:-1:1,:); s_high.TANBC(2:end,:)];   % [NL x Mp] 
states.TANBIC    = [s_low.TANBIC(end:-1:1,:);s_high.TANBIC(2:end,:)];   % [NL x Mp] 

states.ALPHA     = [ s_low.ALPHA(end:-1:1,:); s_high.ALPHA(2:end,:)];   % [NL x Mp] 
% states.deltaALPHA1     = [ s_low.deltaALPHA1(end:-1:1,:); s_high.deltaALPHA1(2:end,:)];   % [NL x Mp] 
% states.deltaALPHA2     = [ s_low.deltaALPHA2(end:-1:1,:); s_high.deltaALPHA2(2:end,:)];   % [NL x Mp]
states.CL        = [    s_low.CL(end:-1:1,:);    s_high.CL(2:end,:)];   % [NL x Mp] 
states.CD        = [    s_low.CD(end:-1:1,:);    s_high.CD(2:end,:)];   % [NL x Mp] 

if Duct_flag == 1  % Duct is present 
    states.Gd        = [     s_low.Gd(end:-1:1,:);      s_high.Gd(2:end,:)];   % [NL x 1]
    states.UARINGq   = [s_low.UARINGq(end:-1:1,:); s_high.UARINGq(2:end,:)];   % [NL x 1]
    states.URRINGq   = [s_low.URRINGq(end:-1:1,:); s_high.URRINGq(2:end,:)];   % [NL x 1]
    states.VSRINGq   = [s_low.VSRINGq(end:-1:1,:); s_high.VSRINGq(2:end,:)];   % [NL x 1]
    states.BetaIDq    = [ s_low.BetaIDq(end:-1:1,:);  s_high.BetaIDq(2:end,:)];   % [NL x 1]
    states.CLd       = [    s_low.CLd(end:-1:1,:);     s_high.CLd(2:end,:)];   % [NL x 1]
    states.UARING    = [ s_low.UARING(end:-1:1,:);  s_high.UARING(2:end,:)];   % [NL x Nd]
    states.URRING    = [ s_low.URRING(end:-1:1,:);  s_high.URRING(2:end,:)];   % [NL x Nd]
    states.UADUCT   = [s_low.UADUCT(end:-1:1,:); s_high.UADUCT(2:end,:)];   % [NL x Mp]
    states.TAU       = [    s_low.TAU(end:-1:1,:);     s_high.TAU(2:end,:)];   % [NL x 1]
    states.CTD       = [    s_low.CTD(end:-1:1,:);     s_high.CTD(2:end,:)];   % [NL x 1]
end

states.converged = [    s_low.converged(end:-1:1,:);    s_high.converged(2:end,:)];   % [NL x 1]
states.iteration = [    s_low.iteration(end:-1:1,:);    s_high.iteration(2:end,:)];   % [NL x 1]

states.KT        = [    s_low.KT(end:-1:1,:);    s_high.KT(2:end,:)];   % [NL x 1]
states.KQ        = [    s_low.KQ(end:-1:1,:);    s_high.KQ(2:end,:)];   % [NL x 1]

states.CT        = [    s_low.CT(end:-1:1,:);    s_high.CT(2:end,:)];   % [NL x 1]
states.CQ        = [    s_low.CQ(end:-1:1,:);    s_high.CQ(2:end,:)];   % [NL x 1]
states.CP        = [    s_low.CP(end:-1:1,:);    s_high.CP(2:end,:)];   % [NL x 1]

if abs(design.VMIV - 1) > 1e-8,  % i.e. if VMIV is not equal to 1 
    
states.EFFYo     = [ s_low.EFFYo(end:-1:1,:); s_high.EFFYo(2:end,:)];   % [NL x 1]
states.Ja        = [    s_low.Ja(end:-1:1,:);    s_high.Ja(2:end,:)];   % [NL x 1]
end

states.EFFY      = [  s_low.EFFY(end:-1:1,:);  s_high.EFFY(2:end,:)];   % [NL x 1]
states.VMWV      = [  s_low.VMWV(end:-1:1,:);  s_high.VMWV(2:end,:)];   % [NL x 1]




end % of the function [states] = Analyze(design,LAMBDAall)
% =========================================================================
% =========================================================================

% =========================================================================
% ============================================== FUNCTION Analyze_some(...)
%%
function [some] = Analyze_some(pt,LAMBDA_some,feather)

%%
Crash_count = 0;

% ------------------------------------------------ Unpack pt data structure   
design = pt.design;
input  = pt.input;

% ------------------------------------------------------------------- Flags
Propeller_flag = input.Propeller_flag; % 0 == turbine, 1 == propeller
Viscous_flag   = input.Viscous_flag;   % 0 == viscous forces off (CD = 0), 
                                       % 1 == viscous forces on

% 0 == no hub, 1 == hub
if isfield(input,'Hub_flag'),    Hub_flag = input.Hub_flag;  
                          else   Hub_flag = 0;  end
                          
% 0 == no duct, 1 == duct
if isfield(input,'Duct_flag'),  Duct_flag = input.Duct_flag;  
                         else   Duct_flag = 0;  end

% 0 == do not display plots, 1 == display plots
if isfield(input,'Plot_flag'),  Plot_flag = input.Plot_flag;  
                         else   Plot_flag = 0;  end 

% 0 == Horseshoe(...,Wrench(...)), 1 == Wake_Horseshoe(...)      
if isfield(input,'Wake_flag'),  Wake_flag = input.Wake_flag;  
                         else   Wake_flag = 0;  end
     
if isfield(input,'Chord_flag'),  Chord_flag = input.Chord_flag;  
                          else   Chord_flag = 0;  end
% -------------------------------------------------------------------------

                         
% '------ Performance inputs ------'                         
Z = input.Z;            % [1 x 1] scalar number of blades

if Propeller_flag == 1
    Js = input.Js;
    L = pi/Js;
else
    L     = input.L;
    Js    = pi/L;
end                         

% '------ Computational inputs ------'
if isfield(input,'ITER'), ITER = input.ITER;  else  ITER = 50;   end
if isfield(input,'Rhv' ), Rhv  = input.Rhv;   else  Rhv  = 0.5;  end

% '------ Design geometry ------'
RC         = design.RC;          % [1 x Mp] control point radii
RV         = design.RV;          % [1 x Mp+1] vortex point radii
DR         = design.DR;          % [1 x Mp] difference in vortex point radii
Rhub_oR    = design.Rhub_oR;     % [1 x 1]

% '------ Design state ------'
G          = design.G';          % [1 x Mp] circulation distribution
VAC        = design.VAC;         % [1 x Mp] 
VTC        = design.VTC;         % [1 x Mp] 
UASTAR     = design.UASTAR;      % [1 x Mp] 
UTSTAR     = design.UTSTAR;      % [1 x Mp] 
if Wake_flag == 1
URSTAR     = design.URSTAR;      % [1 x Mp] 
end
VSTAR      = design.VSTAR;       % [1 x Mp] 
TANBIC     = design.TANBIC;      % [1 x Mp] 
VMIV       = design.VMIV;        % [1 x 1]

% '------ 2D section performance ------'
BetaIC     = atan(TANBIC); % [rad]
ALPHA      = 0*RC;         % [rad], zero at design point

CL         = design.CL;          % [1 x Mp] 
CD         = design.CD;          % [1 x Mp] 
CoD        = design.CoD;         % [1 x Mp] 
t0oc       = design.t0oc;        % [1 x Mp]

Mp  = length(RC);



if isfield(input,'ALPHAstall'), ALPHAstall  = input.ALPHAstall; % [rad]
                          else  ALPHAstall  = 8*pi/180;  end
                          
if length(ALPHAstall) == 1
    ALPHAstall = ALPHAstall * ones(size(RC));
end                          

% --------------------------------------------- dCLdALPHA = 2*pi / (1+2/AR)
if isfield(input,'dCLdALPHA' ), 
    
    dCLdALPHA  = input.dCLdALPHA;
    
    % make dCLdALPHA the same length as RC...
    if      length(dCLdALPHA) == 1
                   dCLdALPHA  = dCLdALPHA * ones(size(RC));        
    elseif (length(dCLdALPHA) == length(input.XR)) && (length(dCLdALPHA) ~= length(RC))
                   dCLdALPHA  = pchip(input.XR,dCLdALPHA,RC);
    end
    
else
    XXR = linspace(Rhub_oR,1,100); % temp radii
    
    if Chord_flag == 0 & isfield(input,'XR') & isfield(input,'XCoD')  
        % Compute  Propeller Aspect Ratio (PAR) from blade outline given in XCoD
        PAR = (1-Rhub_oR)^2 / trapz( XXR , InterpolateChord(input.XR,input.XCoD,XXR) ); 
    else
        PAR = (1-Rhub_oR)^2 / trapz( XXR , InterpolateChord(RC,CoD,XXR)  );     
    end
    
    % Lift curve slope:
    dCLdALPHA = (2*pi / (1 + 2/PAR)) * ones(size(RC));
end
% ------------------------------------------------------------------------- 

% ------------------------------------------------------------------------- 
% '------ Duct parameters ------'
if Duct_flag == 1 
    Rduct_oR    = design.Rduct_oR;    % [1 x 1]
    Cduct_oR    = design.Cduct_oR;    % [1 x 1]
    Xduct_oR    = design.Xduct_oR;    % [1 x 1]
    Gd          = design.Gd;          % [1 x 1]  duct circulation on design
    
    XdRING = design.XdRING;             % [1 x Nd], Nd=12 duct vortex rings    
    VARING = design.VARING;             % [1 x Nd], Nd=12 duct vortex rings    
    UARING = design.UARING;             % [1 x Nd], Nd=12 duct vortex rings
    URRING = design.URRING;             % [1 x Nd], Nd=12 duct vortex rings
    GdRING = design.GdRING;             % [1 x Nd], Nd=12 duct vortex rings
 
    DAHIF_times_TANBIC = design.DAHIFtT; % [Nd,Mp], Duct Horseshoe Influence Functions (influence of propeller on duct)
    DRHIF_times_TANBIC = design.DRHIFtT; % [Nd,Mp], Duct Horseshoe Influence Functions (influence of propeller on duct)
 
    UADIF  = design.UADIF;              % [1 x Mp] influence of duct on propeller
    UADUCT = design.UADUCT;             % [1 x Mp] velocity induced at propeller

    if     isfield(input,'CDd'),      CDd   = input.CDd; else CDd = 0.008; end  % duct drag coefficient
    
    dCLdALPHAd = 2*pi;
    
    % ---------------------------------------------- Flow at duct quarter chord
    DAHIFq_times_TANBIC = design.DAHIFq_times_TANBIC;  % [1,Mp]
    DRHIFq_times_TANBIC = design.DRHIFq_times_TANBIC;  % [1,Mp]
    
    for m = 1:Mp;                               
        DAHIFq(:,m) = DAHIFq_times_TANBIC(:,m) / TANBIC(m);
        DRHIFq(:,m) = DRHIFq_times_TANBIC(:,m) / TANBIC(m);
    end     

    UARINGq = (DAHIFq*G)';  
    URRINGq = (DRHIFq*G)'; 


    VSRINGq = sqrt((VARING+UARINGq)^2+(URRINGq)^2);
    
    CLd     = 4*pi*Gd / (VSRINGq*Cduct_oR);
    
    BetaIDq  = atan( -URRINGq/(VARING+UARINGq) );      % inflow angle
    
    % ------------------------------------------ Record on-design variables 
    CLd0     = CLd;
    BetaIDq0 = BetaIDq;

else
    Rduct_oR    = 1;     
    UADUCT      = 0*RC;     
    CTD         = 0;      
end
% ------------------------------------------------------------------------- 




% ------------------------------------------------------------------------- 
% Find section cross-sectional area for angle of attack correction
if Propeller_flag == 0;
    if isfield(input,'Meanline' ), Meanline  = input.Meanline;   else  Meanline  = 'NACA a=0.8 (modified)';   end
    if isfield(input,'Thickness'), Thickness = input.Thickness;  else  Thickness = 'NACA 65A010 (modified)';  end

    XR = input.XR;
    
    if iscell(Meanline)
    
        if length(Meanline) ~= length(XR)
            Crash_count = -47;
        else
            for j = 1:length(XR)
                [junk1, junk2, junk3, junk4, junk5, junk6, As0(j)] = GeometryFoil2D(Meanline{j},Thickness{j});
                
                As0 = pchip(XR,As0, RC);
            end
        end
    else
        [junk1, junk2, junk3, junk4, junk5, junk6, As0] = GeometryFoil2D(Meanline,Thickness);
    end
end
deltaALPHA1 = 0*RC;
deltaALPHA2 = 0*RC;
deltaALPHA  = deltaALPHA1 + deltaALPHA2;
            
% -------------------------------------------------------------------------              


% ---------------------------------------------- Record on-design variables 
LAMBDA0 = L;
CL0     = CL;
CD0     = CD;
BetaIC0 = BetaIC;    % [rad]

       
% ------------------------- Initialize vortex Horseshoe Influence Functions 
if Wake_flag == 0 
    
    [UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBIC,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);
    
else
    %     disp('Beginning to update {UAHIF,UTHIF,URHIF}'), 
    %     [WX,WY,WZ]          =  Wake_Geometry(Mp,RC,RV,TANBIV,VAC,UASTAR,URSTAR,CTP);
    % %     epsilon = (1-Rhub_oR)/(4*Mp);  % vortex core size in Wake_Influence(...) function
    %     epsilon = min(diff(RV))/100;   % vortex core size in Wake_Influence(...) function
    %    %[UAHIF,UTHIF,URHIF] = Wake_Horseshoe(Mp,Z,RC,SCF,WX,WY,WZ,epsilon);   
    %     [UAHIF,UTHIF,URHIF] = Wake_Horseshoe(Mp,Z,TANBIV,RC,RV,SCF,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR,WX,WY,WZ,epsilon);
    %     disp('Done updating {UAHIF,UTHIF,URHIF}'), 
end

%%
% ---------------------------- Initialize output data structure "some.____"
NL = length(LAMBDA_some);

some.part1  = '------ Off-design states ------';
some.L      = LAMBDA_some(:);  % [NL x 1] tip speed ratio 
some.Js     = pi./some.L;      % [NL x 1] advance coefficient 

some.G      = zeros(NL,Mp);    % [NL x Mp] circulation distribution
some.UASTAR = zeros(NL,Mp);    % [NL x Mp] 
some.UTSTAR = zeros(NL,Mp);    % [NL x Mp] 
if Wake_flag == 1
some.URSTAR = zeros(NL,Mp);    % [NL x Mp] 
end
some.VSTAR  = zeros(NL,Mp);    % [NL x Mp] 
some.TANBC  = zeros(NL,Mp);    % [NL x Mp] 
some.TANBIC = zeros(NL,Mp);    % [NL x Mp]

some.part2  = '------ 2D section performance ------';
some.ALPHA       = zeros(NL,Mp);    % [NL x Mp] 
% some.deltaALPHA1 = zeros(NL,Mp);    % [NL x Mp] 
% some.deltaALPHA2 = zeros(NL,Mp);    % [NL x Mp] 
some.CL          = zeros(NL,Mp);    % [NL x Mp] 
some.CD          = zeros(NL,Mp);    % [NL x Mp] 

some.part3   = '------ Duct parameters ------';
if Duct_flag == 1  % Duct is present 
    Nd = length(design.XdRING);
    some.Gd      = zeros(NL,1);     % [NL x 1]  duct circulation
    some.UARINGq = zeros(NL,1);     % [NL x 1]  UA    at duct quarter chord
    some.URRINGq = zeros(NL,1);     % [NL x 1]  UR    at duct quarter chord 
    some.VSRINGq = zeros(NL,1);     % [NL x 1]  VSTAR at duct quarter chord
    some.BetaIDq = zeros(NL,1);     % [NL x 1]  angle at duct quarter chord 
    some.CLd     = zeros(NL,1);     % [NL x 1]  duct lift coefficient
    some.UARING  = zeros(NL,Nd);    % [NL x Nd], induced velocity at duct
    some.URRING  = zeros(NL,Nd);    % [NL x Nd], induced velocity at duct
    some.UADUCT  = zeros(NL,Mp);    % [NL x Mp], induced velocity at propeller
    some.TAU     = zeros(NL,1);     % [NL x 1]  thrust ratio 
    some.CTD     = zeros(NL,1);     % [NL x 1]  CT for duct
end

some.part4  = '------ Numerical metrics ------';
some.converged  = zeros(NL,1);
some.iteration  = zeros(NL,1);


some.part5  = '------ Performance metrics ------';
some.CT     = zeros(NL,1);     % [NL x 1]
some.CQ     = zeros(NL,1);     % [NL x 1]
some.CP     = zeros(NL,1);     % [NL x 1]
some.KT     = zeros(NL,1);     % [NL x 1]
some.KQ     = zeros(NL,1);     % [NL x 1]
some.EFFY   = zeros(NL,1);     % [NL x 1]
some.EFFYo  = zeros(NL,1);     % open water efficiency
some.Ja     = zeros(NL,1);     % == VMIV*Js
some.VMWV   = zeros(NL,1); 

if Crash_count == -47;
    
    disp('<ERROR>')
    disp('<ERROR> Meanline given as cell array but different length as XR...crash avoided.')
    disp('<ERROR>')    
    
    return
end

% ------------------------------------------------ Implement RepairSpline.m
% To smooth data X(RC):  X_smooth = X*Bsmooth;
Bsmooth = RepairSplineMatrix(RC);
% -------------------------------------------------------------------------

%%


% -------------- For each LAMBDA, find the system state, and compute forces
for i = 1:length(LAMBDA_some)
%     % -------------------------------------------- Initialize the new state 
%     disp(['-------- Determining state:  LAMBDA = ',num2str(LAMBDA_some(i)),' , Js = ',num2str(pi/LAMBDA_some(i))]),
%     disp(' '),
%     disp('Press any key to continue...'),
%     save temp
%     pause,
%     disp(' '),

%% 
    if Crash_count > 3 
        disp(['Crash limit of 3 reached...skipping state:  L = ',num2str(LAMBDA_some(i)),' , Js = ',num2str(pi/LAMBDA_some(i))])
        continue
    end

    converged = 1;
    
    L  =    LAMBDA_some(i);
    Js = pi/LAMBDA_some(i);
    
    TANBC  = (VAC)./(L*RC + VTC);

    
    % Allow twice as many iterations for the first state analyzed
    if i == 1
        ITER = 2*ITER;
    elseif i == 2
        ITER = ITER/2;
    end
         
    % -------------------------------------------------------------------------
    if Plot_flag == 1
        close all,
        % ------------------------------------------------- Uncomment to plot G
        Hgamma = fig;                  set(Hgamma,'Position',[130 680 560 420])
            HHgamma = plot(RC,G,'k.-');   
            set(gca,'XLim',[Rhub_oR 1])
            xlabel('RC [ ]','Fontsize',14,'Fontname','Times')
            ylabel('G [ ]' ,'Fontsize',14,'Fontname','Times')
        % ---------------------------------------------------------------------     

        % ------------------------------ Uncomment to plot UASTAR, UTSTAR, & URSTAR
        Hvel = fig;                      set(Hvel,'Position',[695 680 560 420])
            if Duct_flag == 1
                   HHvel = plot(RC,UASTAR,'b.-',RC,UTSTAR,'r.-',RC,UADUCT,'g.-');
            else
                   HHvel = plot(RC,UASTAR,'b.-',RC,UTSTAR,'r.-');
            end
            axis([Rhub_oR 1 -1 1])
            xlabel('RC [ ]'    ,'Fontsize',14,'Fontname','Times')
            ylabel('UASTAR (blue), UTSTAR (red)','Fontsize',14,'Fontname','Times')
        % --------------------------------------------------------------------- 

        % ---------------------------------------------------------------------
        Hbeta = fig;                    set(Hbeta,'Position',[695 160 560 420])
    %     % ------------------------------------ TANBC & TANBIC
    %         HHbeta = plot(RC,TANBC,'k.-',RC,TANBIC,'r.-'); 
    % 
    %         set(gca,'XLim',[Rhub_oR 1])
    %         xlabel('RC [ ]'    ,'Fontsize',14,'Fontname','Times')
    %         ylabel('TANBIC (red), TANBC (black)','Fontsize',14,'Fontname','Times')
    %     % --------------------------------------------------------------------- 
                % ------------------------------------------------------- ALPHA
                HHbeta = plot(RC,ALPHA*180/pi,'r.'); hold on,            
                plot([0 1],[0 0],'k:')
                plot([0 1],[ALPHAstall(1) ALPHAstall(1)]*180/pi,'k:',[0 1],-[ALPHAstall(1) ALPHAstall(1)]*180/pi,'k:')
                axis([0 1 -20 20])
                xlabel('RC [ ]'    ,'Fontsize',14,'Fontname','Times')
                ylabel('ALPHA, [deg]','Fontsize',14,'Fontname','Times')
                % ----------------------------------------------- percent error
    %             HHbeta = plot([],[],'.');
    %             plot( [Mp,Mp],[-0.1,0.1],'k:',2*[Mp,Mp],[-0.1,0.1],'k:',...
    %                 3*[Mp,Mp],[-0.1,0.1],'k:',4*[Mp,Mp],[-0.1,0.1],'k:',...
    %                 5*[Mp,Mp],[-0.1,0.1],'k:',6*[Mp,Mp],[-0.1,0.1],'k:'),
    %             hold on,
    %             text(10     ,-0.09,'VSTAR'),
    %             text(10+1*Mp,-0.09,'ALPHA'),
    %             text(10+2*Mp,-0.09,'CL'),
    %             text(10+3*Mp,-0.09,'G'),
    %             text(10+4*Mp,-0.09,'UASTAR'),
    %             text(10+5*Mp,-0.09,'UTSTAR'),
    %             plot([0, 6*Mp],[0,0],'k:')
    %             axis([0, 6*Mp, -0.1, 0.1]),
    %             ylabel('percent error, [ ]'    ,'Fontsize',14,'Fontname','Times')            
            % ---------------------------------------------------------------------       


        % ---------------------------------------------------------------------
        Hcod = fig;                      set(Hcod,'Position',[130 160 560 420])
                % ---------------------------------------------- residual error
                HHcod = plot([],[],'.');            
                plot( [Mp,Mp],[-0.1,0.1],'k:',2*[Mp,Mp],[-0.1,0.1],'k:',...
                    3*[Mp,Mp],[-0.1,0.1],'k:',4*[Mp,Mp],[-0.1,0.1],'k:',...
                    5*[Mp,Mp],[-0.1,0.1],'k:',6*[Mp,Mp],[-0.1,0.1],'k:'),
                hold on,
                text(10     ,-0.09,'VSTAR'),
                text(10+1*Mp,-0.09,'ALPHA'),
                text(10+2*Mp,-0.09,'CL'),
                text(10+3*Mp,-0.09,'G'),
                text(10+4*Mp,-0.09,'UASTAR'),
                text(10+5*Mp,-0.09,'UTSTAR'),
                plot([0, 6*Mp],[0,0],'k:')
                axis([0, 6*Mp, -0.1, 0.1]),
                ylabel('residual error, [ ]'    ,'Fontsize',14,'Fontname','Times')
            % ---------------------------------------------------------------------

        pause(0.001), drawnow,
    end  % if Plot_flag == 1
    % -------------------------------------------------------------------------    
    

    
    % =====================================================================
    % ==   FIND STATE OF SYSTEM USING NEWTON SOLVER:
    % ==                   iterate to solve residual equations for unknowns
    N_iter   = 1;                            
    ERROR    = 1;
    ERRORtol = 0.005;
    relax    = 0.9;
%    

            % -------------- Initialize linear system of equations matrices
            R  = zeros(4*Mp,1   );            % R  = vector of residuals
            J  =   eye(4*Mp,4*Mp);            % J  = matrix of derivatives
            DX = zeros(4*Mp,1   );            % DX = vector of change in unknowns 
         %  X  = zeros(4*Mp,1   );            % X  = vector of unknowns              

            % X(     1 :   Mp) =  VSTAR(:);   % X1
            % X(1*Mp+1 : 2*Mp) =  ALPHA(:);   % X2
            % X(2*Mp+1 : 3*Mp) =     CL(:);   % X3
            % X(3*Mp+1 : 4*Mp) =      G(:);   % X4
            % -------------------------------------------------------------------------
            
    while N_iter <= ITER & any(abs(ERROR) > ERRORtol)          % (WHILE LOOP N1)
        % disp(['------- Newton iteration: ',num2str(N_iter)])  % status message
        % disp(' '),

        % ---------------------------------- Store last state of the system
         VSTARlast =  VSTAR;
         alphalast =  ALPHA;
            CLlast =     CL;
             Glast =      G;          

        % ------------- Solve Newton problem to update entire blade at once
        for m = 1:Mp  
            % ------------------------------------------------------ Evaluate residuals
            R(       m) =  VSTAR(m) - sqrt((VAC(m)+UADUCT(m)+UASTAR(m)).^2 + (L*RC(m)+VTC(m)+UTSTAR(m)).^2); % R1

            R(1*Mp + m) =  ALPHA(m) - ( BetaIC0(m) - atan(TANBIC(m)) + feather );                         % R2
            
            R(2*Mp + m) =     CL(m) - CLCD_vs_ALPHA(ALPHA(m)+deltaALPHA(m),ALPHAstall(m),CL0(m),CD0(m),dCLdALPHA(m),Propeller_flag);      % R3
            
            R(3*Mp + m) =      G(m) - (1/(2*pi))*VSTAR(m)*CL(m)*CoD(m);                                % R4        
            % -------------------------------------------------------------------------
                        
            
            % --------------------------------------- Evaluate residual derivatives
            % disp('Evaluating residual derivatives matrix, J...')
            % disp(' '),

            for n = 1:Mp
                J(       m,3*Mp + n) =   (- (VAC(m)+UADUCT(m)+UASTAR(m))/sqrt((VAC(m)+UADUCT(m)+UASTAR(m)).^2 + (L*RC(m)+VTC(m)+UTSTAR(m)).^2) ) * UAHIF(m,n) ...
                                       + (-(L*RC(m)+   VTC(m)+UTSTAR(m))/sqrt((VAC(m)+UADUCT(m)+UASTAR(m)).^2 + (L*RC(m)+VTC(m)+UTSTAR(m)).^2) ) * UTHIF(m,n);

                J(1*Mp + m,3*Mp + n) =   ( (1/(1+(TANBIC(m))^2))*(        1/(L*RC(m)+VTC(m)+UTSTAR(m))) ) * UAHIF(m,n) ... 
                                       + (-(1/(1+(TANBIC(m))^2))*(TANBIC(m)/(L*RC(m)+VTC(m)+UTSTAR(m))) ) * UTHIF(m,n);
            end
            
            J(2*Mp + m,1*Mp + m) = - Find_dCLCDdALPHA(ALPHA(m)+deltaALPHA(m),ALPHAstall(m),CL0(m),CD0(m),dCLdALPHA(m),Propeller_flag);

            J(3*Mp + m,       m) = -(1/(2*pi))*   CL(m)*CoD(m);
            J(3*Mp + m,2*Mp + m) = -(1/(2*pi))*VSTAR(m)*CoD(m);
            % ---------------------------------------------------------------------

        end  
        
        % ------------------------- Update Newton solver vector of unknowns
        DX = linsolve(J,-R);

        if any(isnan(DX))  | ~isreal(DX) | any(abs(DX) > 10^4)
             DX = zeros(4*Mp,1);

             disp('<WARNING>')
             disp('<WARNING> DX == NaN or imaginary... crash avoided...')
             disp('<WARNING>')
             N_iter      = 999;
             converged   = 0;
             Crash_count = Crash_count + 1;
        end         

        VSTAR  = VSTAR + relax*DX( 1:Mp      )';
        ALPHA  = ALPHA + relax*DX((1:Mp)+  Mp)';
        CL     = CL    + relax*DX((1:Mp)+2*Mp)';
        G      = G     + relax*DX((1:Mp)+3*Mp);
        % -----------------------------------------------------------------     

        
        
        % ----------------- Update unkowns not implemented in the Newton solver
        % % Angle of attack correction
        % if Propeller_flag == 0;
        % 
        %     deltaALPHA1 = (Z*As0/pi) .* CoD .* t0oc .* (LAMBDA ./ sqrt((VAC+UADUCT+UASTAR).^2 + (RC*LAMBDA).^2) );
        % 
        %     deltaALPHA2 = (atan((VAC+UADUCT+UASTAR)./(RC*LAMBDA+2*UTSTAR)) - atan((VAC+UADUCT+UASTAR)./(RC*LAMBDA)))/2;
        % 
        %     deltaALPHA = deltaALPHA1 + deltaALPHA2;
        % end    

        
        
        % ----------------------------------------------------------------- 
        % Update induced velocities (influence of propeller on propeller)
        UASTAR = (UAHIF*G)';  
        UTSTAR = (UTHIF*G)';  
        
        TANBIC = (VAC + UADUCT + UASTAR)./(L*RC + VTC + UTSTAR);
                
        % Smooth the inflow angle for numerical stability:
        TANBICsmooth = TANBIC * Bsmooth;
        % ----------------------------------------------------------------- 
        
        
        if Wake_flag == 0 

            [UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBICsmooth,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);  % RELEASE CODE
%               [UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBIC      ,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);  % 12/5/11 test code

        else
            % % disp('Beginning to update {UAHIF,UTHIF,URHIF}'), 
            % 
            % [CT,CQ,CP,KT,KQ,EFFY,TAU,CTP,CTH,CTD,KTP] = ...
            % Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,CD,CoD,G,Z,Js,VMIV,Hub_flag,Rhv,CTD);
            % 
            % [junk1,junk2,URSTAR]= Induced_Velocity(Mp,G,UAHIF,UTHIF,URHIF,Gd,UADIF);        
            % 
            % % [WX,WY,WZ]          =  Wake_Geometry(Mp,RC,RV,TANBIV,VAC,UASTAR,URSTAR,CTP);
            % [WX,WY,WZ] = Wake_Geometry2U(VAC,VTC,UASTAR,UTSTAR,RC,RV,Js);
            % 
            % %[UAHIF,UTHIF,URHIF] = Wake_Horseshoe(Mp,Z,RC,SCF,WX,WY,WZ,epsilon);
            % [UAHIF,UTHIF,URHIF] = Wake_Horseshoe(Mp,Z,TANBIV,RC,RV,SCF,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR,WX,WY,WZ,epsilon);
            % %disp('Done updating {UAHIF,UTHIF,URHIF}'), 
        end   

    % -------------------------------------------------------------------------          
        
        
        
        % ------------------------------------------ Update duct parameters
        if Duct_flag == 1
            % Update the influence of propeller on duct
            for m = 1:Mp;
                DAHIFq(:,m) = DAHIFq_times_TANBIC(:,m) / TANBICsmooth(m);
                DRHIFq(:,m) = DRHIFq_times_TANBIC(:,m) / TANBICsmooth(m);
            end
            
            % Update induced velocities at the duct (influence of propeller on duct)
            UARINGq = (DAHIFq*G)';  
            URRINGq = (DRHIFq*G)'; 
            
            
            VSRINGq = sqrt((VARING+UARINGq)^2+(URRINGq)^2);     % inflow speed

            BetaIDq = atan( -URRINGq/(VARING+UARINGq) );        % inflow angle

            CLd     = dCLdALPHAd * (BetaIDq0 - BetaIDq) + CLd0; % lift coefficient 
            
            Gd      = CLd * VSRINGq*Cduct_oR/(4*pi);            % duct circulation
            
            % Update the induced velocities at the propeller (influence of duct on propeller)
            UADUCT =  UADIF*Gd;
        end
        % -------------------------------------- END update duct parameters
        
        % ------------------------------------------------ End update variables
        % ----------------------------------------------------------------- 
        
        % ------------------------------------------ Evaluate normalized residuals        
        CLtemp = CLCD_vs_ALPHA(ALPHA+deltaALPHA,ALPHAstall,CL0,CD0,dCLdALPHA,Propeller_flag);

        
        ERROR    = [ VSTAR -  sqrt((VAC+UADUCT+UASTAR).^2 + (L*RC+VTC+UTSTAR).^2), ...
                     ALPHA -  (BetaIC0 - atan(TANBIC) + feather), ...
                        CL -  CLtemp, ... 
                        G' -  (1/(2*pi))*VSTAR.*CL.*CoD] ./ ...
                   [max(abs(VSTAR ),1e-4), ...
                    max(abs(ALPHA ),1e-2), ...
                    max(abs(CL    ),1e-4), ...
                    max(abs(G'    ),1e-6)]; 
        % ------------------------------------------ END evaluate normalized residuals
        

        % % ------------------------------------------ Evaluate percent error
        % PERROR    = [ VSTARlast -  VSTAR, ...
        %              alphalast -  ALPHA, ...
        %                 CLlast -     CL, ...
        %                 Glast' -     G'] ./ ...
        %            [max(abs(VSTAR ),1e-4), ...
        %             max(abs(ALPHA ),1e-2), ...
        %             max(abs(CL    ),1e-4), ...
        %             max(abs(G'    ),1e-6)];
        % % ---------------------------------------------- END evaluate percent error


        % --------------------------------------------------- Plot system state
        if Plot_flag == 1
            delete(HHgamma)
            delete(HHvel  )
            delete(HHbeta )
            delete(HHcod  )
            pause(0.0001)
            
            % ---------------------------------------------------- Figure 1
            figure(Hgamma),
                HHgamma = plot(RC,G,'b');
            % ------------------------------------------------------------- 

            % ---------------------------------------------------- Figure 2
            figure(Hvel),               
                if Duct_flag == 1
                       HHvel = plot(RC,UASTAR,'b.-',RC,UTSTAR,'r.-',RC,UADUCT,'g.-');
                else
                       HHvel = plot(RC,UASTAR,'b.-',RC,UTSTAR,'r.-');
                end                
            % --------------------------------------------------------------------- 

            % -------------------------------------------------------- Figure 3
            figure(Hbeta);
                % ------------------------------------------------------ TANBIC
                %  HHbeta = plot(RC,TANBIC0,'r.'); hold on,
                % ------------------------------------------------------- ALPHA
                HHbeta = plot(RC,ALPHA*180/pi,'r'); hold on,        
                % ----------------------------------------------- percent error
                %  HHbeta = plot(1:length(PERROR),PERROR,'.b');
            % -----------------------------------------------------------------


            % ---------------------------------------------------- Figure 4
            figure(Hcod),
                % ------------------------------------------ residual error
                HHcod = plot(1:length(ERROR),ERROR,'b'); 
                
            % ---------------------------------------------------------------------

            pause(0.00001), drawnow, pause(0.0001),
%                  pause,
        end  % --------------------------------------------- END if Plot_flag == 1 


        % -------------------------------------- Prepare for the next iteration
        N_iter  = N_iter + 1;              % iteration in the N loop

%         if N_iter-1 < 10
%             disp(['The max error for iteration  ',num2str(N_iter-1),' is: ',num2str(max(abs(ERROR)))]),  
%            % disp(['Gd = ',num2str(Gd)]),  
%         else
%             disp(['The max error for iteration ' ,num2str(N_iter-1),' is: ',num2str(max(abs(ERROR)))]),  
%            % disp(['Gd = ',num2str(Gd)]),  
%         end
%         % disp(' '),

    end                                                   % (END WHILE LOOP N1)
    
    if N_iter > ITER  &  converged == 1
        disp('WARNING: Newton solver did NOT converge.'),
        disp(' ')
        
        converged   = 0;
        Crash_count = Crash_count + 1;
    end
    
    N_iter  = N_iter - 1;

    % ---------------------------------------------------------------------
    %             Many journeys end here                                  %
    %             But, the secret's told the same:                        %
    %             Life is just a candle, and a dream                      %
    %             Must give it flame                                      %
    % ---------------------------------------------------------------------
    
    % =================================================== END NEWTON SOLVER  
    % =====================================================================
%%
% save temp
% return
    
    % ------------------------------------------------------ Compute forces
    
    % ---------------------------------------------------------------------
    if Duct_flag == 1
        % Update the influence of propeller on duct
        for m = 1:Mp;
            DAHIF(:,m) = DAHIF_times_TANBIC(:,m) / TANBICsmooth(m);
            DRHIF(:,m) = DRHIF_times_TANBIC(:,m) / TANBICsmooth(m);
        end

        % Update induced velocities at the duct (influence of propeller on duct)
        UARING = (DAHIF*G)';  
        URRING = (DRHIF*G)'; 

        % Find duct thrust, (CTDdes==1.0 is unused here)            
        [CTD,junk] = Duct_Thrust(XdRING,Rduct_oR,VARING,UARING,URRING,GdRING,Gd,CDd,1.0);                      
    end
    % ---------------------------------------------------------------------     
    
    
    [CL,CD] = CLCD_vs_ALPHA(ALPHA+deltaALPHA,ALPHAstall,CL0,CD0,dCLdALPHA,Propeller_flag);

    
    if Viscous_flag == 0
        CD = 0*CD;
    end

    
    [CT,CQ,CP,KT,KQ, CTH,TAU, Ja,Jw,VMWV, EFFYo, EFFY,ADEFFY,QF] = ...
          Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,UADUCT,CD,CoD,G,Z,Js,VMIV,Hub_flag,Rhub_oR,Rhv,CTD);
 
 
    
    % ---------------------------------------------------------------------     
    if     (Propeller_flag == 1) & (KT < 0 | KQ < 0)
        Crash_count = Crash_count + 1;
        
    elseif (Propeller_flag == 0) & (CP > 0)
        Crash_count = Crash_count + 1;
    end
    

    % ------------------------- Record output data in structure "some.____"
    some.converged(i) = converged;
    some.iteration(i) = N_iter;
           
         some.G(i,:) = G';       % [NL x Mp] circulation distribution
    some.UASTAR(i,:) = UASTAR;   % [NL x Mp] 
    some.UTSTAR(i,:) = UTSTAR;   % [NL x Mp]
    if Wake_flag == 1
    some.URSTAR(i,:) = URSTAR;   % [NL x Mp]
    end
     some.VSTAR(i,:) = VSTAR;    % [NL x Mp]     
     some.TANBC(i,:) = TANBC;    % [NL x Mp] 
    some.TANBIC(i,:) = TANBIC;   % [NL x Mp] 
    
      some.ALPHA(i,:) = ALPHA;          % [NL x Mp] 
% some.deltaALPHA1(i,:) = deltaALPHA1;    % [NL x Mp] 
% some.deltaALPHA2(i,:) = deltaALPHA2;    % [NL x Mp]      
         some.CL(i,:) = CL;             % [NL x Mp] 
         some.CD(i,:) = CD;             % [NL x Mp] 
        
    if Duct_flag == 1  % Duct is present,  Nd = length(design.Xduct_oR);
        some.Gd(i,:)      = Gd;      % [NL x 1]  duct circulation
        
        some.UARINGq(i,:) = UARINGq; % [NL x 1]  UA    at duct quarter chord
        some.URRINGq(i,:) = URRINGq; % [NL x 1]  UR    at duct quarter chord 
        some.VSRINGq(i,:) = VSRINGq; % [NL x 1]  VSTAR at duct quarter chord
        some.BetaIDq(i,:) = BetaIDq;  % [NL x 1]  angle at duct quarter chord 
        
        some.CLd(i,:)     = CLd;     % [NL x 1]  duct lift coefficient
        
        some.UARING(i,:)  = UARING;  % [NL x Nd], induced velocity at duct
        some.URRING(i,:)  = URRING;  % [NL x Nd], induced velocity at duct
        
        some.UADUCT(i,:)  = UADUCT; % [NL x Mp], induced velocity at propeller
        some.TAU(i,:)     = TAU;     % [NL x 1]  thrust ratio 
        some.CTD(i,:)     = CTD;     % [NL x 1]  CT for duct
    end
     
        some.KT(i) = KT;       % [NL x 1]
        some.KQ(i) = KQ;       % [NL x 1]
        some.CT(i) = CT;       % [NL x 1]
        some.CQ(i) = CQ;       % [NL x 1]
        some.CP(i) = CP;       % [NL x 1]
      some.EFFY(i) = EFFY;     % [NL x 1]

    some.EFFYo(i)  = EFFYo; % open water efficiency
       some.Ja(i)  = Ja; % == VMIV*Js

     some.VMWV(i)  = VMWV;               % [1 x 1]

       
       
% -------------------------------------------------------------------------
    disp(' '),
    disp(['Forces for state:  L = ',num2str(LAMBDA_some(i)),' , Js = ',num2str(pi/LAMBDA_some(i))]),
    disp(' '),
if Propeller_flag == 1
   
    disp(['    Js = ',num2str(Js)]),    
    disp(['    KT = ',num2str(KT)]),   
    disp(['    KQ = ',num2str(KQ)]),   
    disp(['    CT = ',num2str(CT)]),
    if abs(VMIV - 1) > 1e-8,  % i.e. if VMIV is not equal to 1 
    disp(['    Ja = ',num2str(Ja)]),
    end
    disp(['  EFFY = ',num2str(EFFY)]),   
    
else  
    disp(['CT = ' ,num2str(CT)]), 
    disp(['CP = ' ,num2str(CP)]),       
end
    disp(' ')     
% -------------------------------------------------------------------------       
   
    
    if N_iter < 10
        disp(['The max error for iteration  ',num2str(N_iter),' is: ',num2str(max(abs(ERROR)))]),  
    else
        disp(['The max error for iteration ' ,num2str(N_iter),' is: ',num2str(max(abs(ERROR)))]),  
    end
    disp(' '),    
    disp('------------------------------------------------'),
    disp(' ')
   
    
%     if converged == 0
%         disp('WARNING: Newton solver did NOT converge.'),
%         disp(' ')
%         disp('------------------------------------------------'),
%         disp(' ')
%     else    
%         disp('------------------------------------------------'),
%         disp(' ')
%     end    
    
%     disp(' ')
%     disp('Analyze(...) paused.  Press any key to continue... ')
%     pause,
    disp(' ')
    disp(' ')
    
end
    
end % function [...] = Analyze_some(...)                                
