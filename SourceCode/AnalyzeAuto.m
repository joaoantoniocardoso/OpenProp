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
% ======================================================== AnalyzeAuto Function
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
%       delta   	resolution of performance curves (in units of Js or L)
% 
%       Lmin        minimum tip speed ratio to analyze
%
% Outputs:   
%       states      data structure with operating states
%
% -------------------------------------------------------------------------
% Method:
%   Automatically determines Js or L to find entire performance curve.
%
% -------------------------------------------------------------------------

function [s] = AnalyzeAuto(pt,delta,Lmin,CTmin,CTmax)

if nargin == 1
    delta = 0;
    Lmin  = 0;
    CTmin = 0;
    CTmax = Inf;
elseif nargin == 2
    Lmin  = 0;
    CTmin = 0;
    CTmax = Inf;
elseif nargin == 3
    CTmin = 0;
    CTmax = Inf;
end

warning off
disp('------------------------------------------------'),
disp(' ')

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% ------------------------------------------------ Unpack pt data structure
if isfield(pt,'i') & ~isfield(pt,'input' ), input  = pt.i; else input  = pt.input;  end
if isfield(pt,'d') & ~isfield(pt,'design'), design = pt.d; else design = pt.design; end

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

% 0 == Horseshoe(...,Wrench(...)), 1 == Wake_Horseshoe(...)      
if isfield(input,'Wake_flag'),  Wake_flag = input.Wake_flag;  
                         else   Wake_flag = 0;  end
                         
if isfield(input,'Chord_flag'), Chord_flag = input.Chord_flag;  
                         else   Chord_flag = 0;  end
                         
                         
                         
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------- Perform error checking
if Wake_flag == 1

    disp('ERROR: AnalyzeAuto is not set up for generlized the wake model.')
    return
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
                         

          
% ------------------------------------------------------ Performance inputs
Z = input.Z;            % [1 x 1] scalar number of blades

if Propeller_flag == 1
    Js = input.Js;
    L  = pi/Js;
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
VSTAR      = design.VSTAR;       % [1 x Mp] 
TANBIC     = design.TANBIC;      % [1 x Mp] 
VMIV       = design.VMIV;        % [1 x 1]

% '------ 2D section performance ------'
ALPHA      = 0*RC;   % [1 x Mp], [rad]

CL         = design.CL;          % [1 x Mp] 
CD         = design.CD;          % [1 x Mp] 
CoD        = design.CoD;         % [1 x Mp] 

Mp = length(RC);

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
    % Find propeller influence at the duct quarter chord, and
    % put these in pt.design, i.e. where they can be accessed from Analyze_some(...)
    disp(' '), disp('Computing rotor-duct interaction...be patient...'), disp(' '),
    [DAHIFq_times_TANBIC,junk,DRHIFq_times_TANBIC] = Horseshoe_intr_110830(Xduct_oR-Cduct_oR/4,Rduct_oR ,RC,ones(size(RC)),RV,Z,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);    
    
    for m = 1:Mp;                               
        DAHIFq(:,m) = DAHIFq_times_TANBIC(:,m) / TANBIC(m);
        DRHIFq(:,m) = DRHIFq_times_TANBIC(:,m) / TANBIC(m);
    end     

    UARINGq = (DAHIFq*G)';  
    URRINGq = (DRHIFq*G)'; 


    VSRINGq = sqrt((VARING+UARINGq)^2+(URRINGq)^2);
    
    CLd     = 4*pi*Gd / (VSRINGq*Cduct_oR);
    
    BetaIDq  = atan( -URRINGq/(VARING+UARINGq) );      % inflow angle
    
    CTD      = design.CTD;         % [1 x 1]  duct thrust coefficient on design
    
else
    % Variables that never change:
    Rduct_oR            = 1;           
    Cduct_oR            = 1;
    Xduct_oR            = 0;             
    XdRING              = 0;  
    VARING              = 0;    
    GdRING              = 0;
    DAHIF_times_TANBIC  = 0;
    DRHIF_times_TANBIC  = 0;
    UADIF               = 0*RC;
    CDd                 = 0;    
    dCLdALPHAd          = 0;
    DAHIFq_times_TANBIC = 0;
    DRHIFq_times_TANBIC = 0; 
    
    % Variables that determine the operating state (inputs that need to be reset)
    UARINGq             = 0;  
    URRINGq             = 0; 
    VSRINGq             = 0;
    BetaIDq             = 0;    
    CLd                 = 0;
    Gd                  = 0;
    UADUCT              = 0*RC;
    
    % Variables that describe  the operating state (performance outputs only)
    UARING              = 0;
    URRING              = 0;
    CTD                 = 0; 
end
% -------------------------------------------------------------------------

    
% -------------------------------------------------------------------------
% Record on-design variables
TANBIC0 = TANBIC;
CL0     = CL;
CD0     = CD;
    % Duct on-design variables
    UARINGq0             = UARINGq;  
    URRINGq0             = URRINGq; 
    VSRINGq0             = VSRINGq;
    BetaIDq0             = BetaIDq;    
    CLd0                 = CLd;
    Gd0                  = Gd;
    UADUCT0              = UADUCT;
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% ------------------------------- Initialize output data structure "s.____"
% NL == number of off-design states, unknown a-priori

s.part1   = '------ Off-design states ------';
s.L       = [];    % [NL x 1] tip speed ratio 
s.Js      = [];    % [NL x 1] advance coefficient 

s.RC      = RC;    % [ 1 x Mp] radius / propeller radius
s.CoD     = CoD;   % [ 1 x Mp] chord  / propeller diameter

s.G       = [];    % [NL x Mp] circulation distribution
s.UASTAR  = [];    % [NL x Mp] axial      induced velocity distribution
s.UTSTAR  = [];    % [NL x Mp] tangential induced velocity distribution
% if Wake_flag == 1
% s.URSTAR  = [];    % [NL x Mp] radial     induced velocity distribution
% end
s.VSTAR   = [];    % [NL x Mp] total inflow       velocity distribution
s.TANBC   = [];    % [NL x Mp] tangent of free-stream inflow angle
s.TANBIC  = [];    % [NL x Mp] tangent of total inflow angle

s.ALPHA       = [];    % [NL x Mp] ALPHA = alpha - alphaI
% s.deltaALPHA1 = [];    % [NL x Mp] 
% s.deltaALPHA2 = [];    % [NL x Mp]
s.CL          = [];    % [NL x Mp] lift coefficient distribution
s.CD          = [];    % [NL x Mp] drag coefficient distribution

                        
if Duct_flag == 1  % Duct is present 
    s.part3   = '------ Duct parameters ------'; 
    Nd        = length(design.XdRING);
    s.Gd      = [];     % [NL x 1]  duct circulation
    s.UARINGq = [];     % [NL x 1]  UA    at duct quarter chord
    s.URRINGq = [];     % [NL x 1]  UR    at duct quarter chord 
    s.VSRINGq = [];     % [NL x 1]  VSTAR at duct quarter chord
    s.BetaIDq = [];     % [NL x 1]  angle at duct quarter chord 
    s.CLd     = [];     % [NL x 1]  duct lift coefficient
    s.UARING  = [];     % [NL x Nd], induced velocity at duct
    s.URRING  = [];     % [NL x Nd], induced velocity at duct
    s.UADUCT  = [];     % [NL x Mp], induced velocity at propeller
    s.TAU     = [];     % [NL x 1]  thrust ratio 
    s.CTD     = [];     % [NL x 1]  CT for duct

    s.part3   = '------ Performance metrics ------';
else
    s.part4   = '------ Performance metrics ------';
end

s.KT      = [];     % [NL x 1]  thrust coefficient
s.KQ      = [];     % [NL x 1]  torque coefficient
s.CT      = [];     % [NL x 1]  thrust coefficient
s.CQ      = [];     % [NL x 1]  torque coefficient
s.CP      = [];     % [NL x 1]  power  coefficient
s.EFFY    = [];     % [NL x 1]  efficiency
s.VMWV    = [];

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% save temp
% return
%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% ---------------------------------------------- Analyze design state first
Lnew = L;  % design state: (L,G,UASTAR,...), and new state: (Lnew,...)

[L,Js,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,TANBC,ALPHA,CL,Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT,KT,KQ,CT,CQ,CP,EFFY,VMWV,TAU,CTD,converged] =  AnalyzeSingleState(Propeller_flag,Viscous_flag,Hub_flag,Duct_flag,Z,Mp,ITER,Rhv,ALPHAstall,dCLdALPHA,RC,RV,DR,Rhub_oR,CoD,VMIV,TANBIC0,CL0,CD0,L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,ALPHA,CL,Rduct_oR,Cduct_oR,Xduct_oR,XdRING,VARING,GdRING,DAHIF_times_TANBIC,DRHIF_times_TANBIC,UADIF,CDd,dCLdALPHAd,DAHIFq_times_TANBIC,DRHIFq_times_TANBIC, CLd0,BetaIDq0, Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UADUCT, Lnew);

s = AppendState(s,L,G,UASTAR,UTSTAR,VSTAR,TANBC,TANBIC,ALPHA,CL,CD, KT,KQ,CT,CQ,CP,EFFY,VMWV,  Duct_flag,Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT,TAU,CTD);

%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% ------------ Analyze lower LAMBDA states first, then higher LAMBDA states 
%                   (i.e. higher Js states first, then lower Js states)
KeepIterating_flag = 1; 

if  delta ~= 0
    deltaJ = delta;
    deltaL = delta;
else    
    deltaJ = 0.05;
    deltaL = 0.10;
end

while KeepIterating_flag == 1 
       
    % ---------------------------------------------- Increment Js or LAMBDA
    if Propeller_flag == 1
        Js = Js + deltaJ;

        Lnew = pi/Js;

    else
        Lnew = L - deltaL;
        
        Js = pi/Lnew;
    end
    % ---------------------------------------------------------------------
    
    % --------------------------------------------------- Analyze new state
    [L,Js,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,TANBC,ALPHA,CL,Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT,KT,KQ,CT,CQ,CP,EFFY,VMWV,TAU,CTD,converged] =  AnalyzeSingleState(Propeller_flag,Viscous_flag,Hub_flag,Duct_flag,Z,Mp,ITER,Rhv,ALPHAstall,dCLdALPHA,RC,RV,DR,Rhub_oR,CoD,VMIV,TANBIC0,CL0,CD0,L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,ALPHA,CL,Rduct_oR,Cduct_oR,Xduct_oR,XdRING,VARING,GdRING,DAHIF_times_TANBIC,DRHIF_times_TANBIC,UADIF,CDd,dCLdALPHAd,DAHIFq_times_TANBIC,DRHIFq_times_TANBIC, CLd0,BetaIDq0, Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UADUCT, Lnew);
    % ---------------------------------------------------------------------
    
    % --------------------------------- Check if loop should keep iterating
    if Propeller_flag == 1
        KeepIterating_flag = converged & (CQ > 0) & (Lnew >= Lmin) & (CT > CTmin) & (CT < CTmax);
    else
        KeepIterating_flag = converged & (CP < 0) & (Lnew >= Lmin);
    end
    % ---------------------------------------------------------------------
    
    % ----------------- Append new state, or revert to last converged state
    if KeepIterating_flag == 1
        
        % Append the converged state to the states data structure
        s = AppendState(s,L,G,UASTAR,UTSTAR,VSTAR,TANBC,TANBIC,ALPHA,CL,CD, KT,KQ,CT,CQ,CP,EFFY,VMWV,  Duct_flag,Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT,TAU,CTD);
    else
        
        % Reset state variables to last converged values
        NL = length(s.L);

        L       =  s.L(NL);
        Js      = s.Js(NL);
        
        G       =      s.G(NL,:)';
        UASTAR  = s.UASTAR(NL,:);
        UTSTAR  = s.UTSTAR(NL,:);
        VSTAR   =  s.VSTAR(NL,:);
        TANBIC  = s.TANBIC(NL,:);
        ALPHA   =  s.ALPHA(NL,:);
        CL      =     s.CL(NL,:);
        CD      =     s.CD(NL,:);
        
    end
    % ---------------------------------------------------------------------
    
end

%%
% -------------------------------------------------------------------------
% -------------------------------- Resolve the end of the performance curve
KeepTrying_flag = (Lnew >= Lmin); 

while KeepTrying_flag == 1

    % -------------------------------------------------------------------------
    % -------------------- Estimate Js for which KT < 0 (or L for which CP > 0)
    if Propeller_flag == 1
        Js0 = interp1(s.KT,s.Js,0,'linear','extrap');

        deltaJsmall = min([ abs(Js0 - Js)/5,  deltaJ/5]);
        

    else
        L0 = interp1(s.CP,s.L,0,'linear','extrap');

        deltaLsmall = min([ abs(L0 - L)/5,  deltaL/5]);
    end
    
    % -------------------------------------------------------------------------
    KeepIterating_flag = (Lnew >= Lmin); 
    NumberOfNewStates  = 0;


    while KeepIterating_flag == 1 

        % ---------------------------------------------- Increment Js or LAMBDA
        if Propeller_flag == 1
            Js = Js + deltaJsmall;

            Lnew = pi/Js;

        else
            Lnew = L - deltaLsmall;

            Js = pi/Lnew;
        end
        % ---------------------------------------------------------------------

        % --------------------------------------------------- Analyze new state
        [L,Js,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,TANBC,ALPHA,CL,Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT,KT,KQ,CT,CQ,CP,EFFY,VMWV,TAU,CTD,converged] =  AnalyzeSingleState(Propeller_flag,Viscous_flag,Hub_flag,Duct_flag,Z,Mp,ITER,Rhv,ALPHAstall,dCLdALPHA,RC,RV,DR,Rhub_oR,CoD,VMIV,TANBIC0,CL0,CD0,L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,ALPHA,CL,Rduct_oR,Cduct_oR,Xduct_oR,XdRING,VARING,GdRING,DAHIF_times_TANBIC,DRHIF_times_TANBIC,UADIF,CDd,dCLdALPHAd,DAHIFq_times_TANBIC,DRHIFq_times_TANBIC, CLd0,BetaIDq0, Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UADUCT, Lnew);
        % ---------------------------------------------------------------------

        % --------------------------------- Check if loop should keep iterating
        if Propeller_flag == 1
            KeepIterating_flag = converged & (CQ > 0) & (Lnew >= Lmin) & (CT > CTmin) & (CT < CTmax);
        else
            KeepIterating_flag = converged & (CP < 0) & (Lnew >= Lmin);
        end
        % ---------------------------------------------------------------------

        % ----------------- Append new state, or revert to last converged state
        if KeepIterating_flag == 1

            % Append the converged state to the states data structure
            s = AppendState(s,L,G,UASTAR,UTSTAR,VSTAR,TANBC,TANBIC,ALPHA,CL,CD, KT,KQ,CT,CQ,CP,EFFY,VMWV,  Duct_flag,Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT,TAU,CTD);
        
            NumberOfNewStates = NumberOfNewStates + 1;
        else

            % Reset state variables to last converged values
            NL = length(s.L);

            L       =  s.L(NL);
            Js      = s.Js(NL);

            G       =      s.G(NL,:)';
            UASTAR  = s.UASTAR(NL,:);
            UTSTAR  = s.UTSTAR(NL,:);
            VSTAR   =  s.VSTAR(NL,:);
            TANBIC  = s.TANBIC(NL,:);
            ALPHA   =  s.ALPHA(NL,:);
            CL      =     s.CL(NL,:);
            CD      =     s.CD(NL,:);
        end
        % ---------------------------------------------------------------------
        
    end
    
    % --------------------------------- Check if loop should keep iterating
    EFFYtol = 0.1;
      CPtol = 0.1;
    
    if Propeller_flag == 1
        KeepTrying_flag = (s.EFFY(NL) > EFFYtol) & (Lnew >= Lmin) & (NumberOfNewStates > 0);
    else
        KeepTrying_flag = (- s.CP(NL) >   CPtol) & (Lnew >= Lmin) & (NumberOfNewStates > 0);
    end
    % ---------------------------------------------------------------------        
end

%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Revert to design state to analyze the other half of the performance curve
% (L,G,UASTAR,UTSTAR,VSTAR,TANBIC,TANBIV,ALPHA,CL,CD)
if Propeller_flag == 1
    Js = input.Js;
    L = pi/Js;
else
    L     = input.L;
    Js    = pi/L;    
end   
G          = design.G';          % [1 x Mp] 
UASTAR     = design.UASTAR;      % [1 x Mp] 
UTSTAR     = design.UTSTAR;      % [1 x Mp] 
VSTAR      = design.VSTAR;       % [1 x Mp] 
TANBIC     = design.TANBIC;      % [1 x Mp] 
ALPHA      = 0*RC;               % [1 x Mp], [rad]
CL         = design.CL;          % [1 x Mp] 
CD         = design.CD;          % [1 x Mp] 


% '------ Duct parameters ------'
if Duct_flag == 1 
    UARINGq = UARINGq0;  
    URRINGq = URRINGq0; 
    VSRINGq = VSRINGq0;
    BetaIDq = BetaIDq0;    
    CLd     = CLd0;
    Gd      = Gd0;
    UADUCT  = UADUCT0;    
end
    
% ---------------------------------------------- Analyze design state first
Lnew = L;  % design state: (L,G,UASTAR,...), and new state: (Lnew,...)

[L,Js,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,TANBC,ALPHA,CL,Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT,KT,KQ,CT,CQ,CP,EFFY,VMWV,TAU,CTD,converged] =  AnalyzeSingleState(Propeller_flag,Viscous_flag,Hub_flag,Duct_flag,Z,Mp,ITER,Rhv,ALPHAstall,dCLdALPHA,RC,RV,DR,Rhub_oR,CoD,VMIV,TANBIC0,CL0,CD0,L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,ALPHA,CL,Rduct_oR,Cduct_oR,Xduct_oR,XdRING,VARING,GdRING,DAHIF_times_TANBIC,DRHIF_times_TANBIC,UADIF,CDd,dCLdALPHAd,DAHIFq_times_TANBIC,DRHIFq_times_TANBIC, CLd0,BetaIDq0, Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UADUCT, Lnew);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if  delta ~= 0
    deltaJ = delta;
    deltaL = delta;
else    
    deltaJ = 0.05;
    deltaL = 0.5;
end
            
% disp('success')
%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --------------------- Analyze higher LAMBDA states (i.e. lower Js states)
KeepIterating_flag = 1; 

while KeepIterating_flag == 1 
       
    % ---------------------------------------------- Increment Js or LAMBDA
    if Propeller_flag == 1
        Js = Js - deltaJ;

        Lnew = pi/Js;

    else
        Lnew = L + deltaL;
        
        Js = pi/Lnew;
    end
    % ---------------------------------------------------------------------
    
    % --------------------------------------------------- Analyze new state
    [L,Js,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,TANBC,ALPHA,CL,Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT,KT,KQ,CT,CQ,CP,EFFY,VMWV,TAU,CTD,converged] =  AnalyzeSingleState(Propeller_flag,Viscous_flag,Hub_flag,Duct_flag,Z,Mp,ITER,Rhv,ALPHAstall,dCLdALPHA,RC,RV,DR,Rhub_oR,CoD,VMIV,TANBIC0,CL0,CD0,L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,ALPHA,CL,Rduct_oR,Cduct_oR,Xduct_oR,XdRING,VARING,GdRING,DAHIF_times_TANBIC,DRHIF_times_TANBIC,UADIF,CDd,dCLdALPHAd,DAHIFq_times_TANBIC,DRHIFq_times_TANBIC, CLd0,BetaIDq0, Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UADUCT, Lnew);
    % ---------------------------------------------------------------------
    
    % --------------------------------- Check if loop should keep iterating
    if Propeller_flag == 1
        KeepIterating_flag = converged & (CQ > 0) & (CT > CTmin) & (CT < CTmax);
    else
        KeepIterating_flag = converged & (CP < 0);
    end
    % ---------------------------------------------------------------------
    
    % ----------------- Append new state, or revert to last converged state
    if KeepIterating_flag == 1
        
        % Append the converged state to the states data structure
        s = AppendState(s,L,G,UASTAR,UTSTAR,VSTAR,TANBC,TANBIC,ALPHA,CL,CD, KT,KQ,CT,CQ,CP,EFFY,VMWV,  Duct_flag,Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT,TAU,CTD);
    else
        
        % Reset state variables to last converged values
        NL = length(s.L);

        L       =  s.L(NL);
        Js      = s.Js(NL);
        
        G       =      s.G(NL,:)';
        UASTAR  = s.UASTAR(NL,:);
        UTSTAR  = s.UTSTAR(NL,:);
        VSTAR   =  s.VSTAR(NL,:);
        TANBIC  = s.TANBIC(NL,:);
        ALPHA   =  s.ALPHA(NL,:);
        CL      =     s.CL(NL,:);
        CD      =     s.CD(NL,:);
    end
    % ---------------------------------------------------------------------
    
end


%%
% -------------------------------------------------------------------------
% -------------------------------- Resolve the end of the performance curve


% -------------------------------------------------------------------------
% -------------------- Estimate Js for which KT < 0 (or L for which CP > 0)
if Propeller_flag == 1
    Js0 = interp1(s.KT,s.Js,0,'linear','extrap');
    
    deltaJsmall = min([ abs(Js0 - Js)/5,  deltaJ/5]);

else
    L0 = interp1(s.CP,s.L,0,'linear','extrap');

    deltaLsmall = min([ abs(L0 - L)/5,  deltaL/5]);
end

% -------------------------------------------------------------------------
KeepIterating_flag = 1; 

while KeepIterating_flag == 1 
       
    % ---------------------------------------------- Increment Js or LAMBDA
    if Propeller_flag == 1
        Js = Js - deltaJsmall;

        Lnew = pi/Js;

    else
        Lnew = L + deltaLsmall;
        
        Js = pi/Lnew;
    end
    % ---------------------------------------------------------------------
    
    % --------------------------------------------------- Analyze new state
    [L,Js,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,TANBC,ALPHA,CL,Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT,KT,KQ,CT,CQ,CP,EFFY,VMWV,TAU,CTD,converged] =  AnalyzeSingleState(Propeller_flag,Viscous_flag,Hub_flag,Duct_flag,Z,Mp,ITER,Rhv,ALPHAstall,dCLdALPHA,RC,RV,DR,Rhub_oR,CoD,VMIV,TANBIC0,CL0,CD0,L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,ALPHA,CL,Rduct_oR,Cduct_oR,Xduct_oR,XdRING,VARING,GdRING,DAHIF_times_TANBIC,DRHIF_times_TANBIC,UADIF,CDd,dCLdALPHAd,DAHIFq_times_TANBIC,DRHIFq_times_TANBIC, CLd0,BetaIDq0, Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UADUCT, Lnew);
    % ---------------------------------------------------------------------
    
    % --------------------------------- Check if loop should keep iterating
    if Propeller_flag == 1
        KeepIterating_flag = converged & (CQ > 0) & (CT > CTmin) & (CT < CTmax);
    else
        KeepIterating_flag = converged & (CP < 0);
    end
    % ---------------------------------------------------------------------
    
    % ----------------- Append new state, or revert to last converged state
    if KeepIterating_flag == 1
        
        % Append the converged state to the states data structure
        s = AppendState(s,L,G,UASTAR,UTSTAR,VSTAR,TANBC,TANBIC,ALPHA,CL,CD, KT,KQ,CT,CQ,CP,EFFY,VMWV,  Duct_flag,Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT,TAU,CTD);
    else
        
        % No need to reset state variables, since done
     
    end
    % ---------------------------------------------------------------------
    
end


%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------    
% Sort states by asending L

[s.L,order] = sort(s.L);

s.Js  = s.Js(order);

s.G      =      s.G(order,:);
s.UASTAR = s.UASTAR(order,:);
s.UTSTAR = s.UTSTAR(order,:);
s.VSTAR  =  s.VSTAR(order,:);
s.TANBC  =  s.TANBC(order,:);
s.TANBIC = s.TANBIC(order,:);
s.ALPHA  =  s.ALPHA(order,:);
s.CL     =     s.CL(order,:);
s.CD     =     s.CD(order,:);

s.KT     = s.KT(order);
s.KQ     = s.KQ(order);
s.CT     = s.CT(order);
s.CQ     = s.CQ(order);
s.CP     = s.CP(order);
s.EFFY   = s.EFFY(order);
s.VMWV   = s.VMWV(order);



if Duct_flag == 1  % Duct is present 
    s.Gd      = s.Gd(order);
    s.UARINGq = s.UARINGq(order);
    s.URRINGq = s.URRINGq(order);
    s.VSRINGq = s.VSRINGq(order);
    s.BetaIDq = s.BetaIDq(order);
    s.CLd     = s.CLd(order);
    s.UARING  = s.UARING(order,:);
    s.URRING  = s.URRING(order,:);
    s.UADUCT  = s.UADUCT(order,:);
    s.TAU     = s.TAU(order); 
    s.CTD     = s.CTD(order);
end

end
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------- 
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------- 



%%
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------- 
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------- 



%%
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------- 
% -------------------------------------------------------------------------
% ------------------------------------------------------------------------- 
% This function appends a state onto states data structure "s"

function s = AppendState(s,L,G,UASTAR,UTSTAR,VSTAR,TANBC,TANBIC,ALPHA,CL,CD, KT,KQ,CT,CQ,CP,EFFY,VMWV,  Duct_flag,Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT,TAU,CTD)

if nargin == 17
    Duct_flag = 0;
end

% s.part1 '------ Off-design states ------';
s.L       = [s.L;     L];                       % [NL x 1] tip speed ratio 
s.Js      = [s.Js; pi/L];                       % [NL x 1] advance coefficient 

s.G       = [s.G;           G'];                % [NL x Mp] circulation distribution
s.UASTAR  = [s.UASTAR; UASTAR];                 % [NL x Mp] axial      induced velocity distribution
s.UTSTAR  = [s.UTSTAR; UTSTAR];                 % [NL x Mp] tangential induced velocity distribution

s.VSTAR   = [s.VSTAR;  VSTAR];                  % [NL x Mp] total inflow       velocity distribution
s.TANBC   = [s.TANBC;  TANBC];                  % [NL x Mp] tangent of free-stream inflow angle
s.TANBIC  = [s.TANBIC; TANBIC];                 % [NL x Mp] tangent of total inflow angle
s.ALPHA   = [s.ALPHA;  ALPHA];                  % [NL x Mp] ALPHA = alpha - alphaI, RADIANS
% % s.deltaALPHA1 = [];                         % [NL x Mp] 
% % s.deltaALPHA2 = [];                         % [NL x Mp]
s.CL      = [s.CL; CL];                         % [NL x Mp] lift coefficient distribution
s.CD      = [s.CD; CD];                         % [NL x Mp] drag coefficient distribution

% s.part3 '------ Duct parameters ------';                         
if Duct_flag == 1  % Duct is present 
    s.Gd      = [s.Gd; Gd];                     % [NL x 1]  duct circulation
    s.UARINGq = [s.UARINGq; UARINGq];           % [NL x 1]  UA    at duct quarter chord
    s.URRINGq = [s.URRINGq; URRINGq];           % [NL x 1]  UR    at duct quarter chord 
    s.VSRINGq = [s.VSRINGq; VSRINGq];           % [NL x 1]  VSTAR at duct quarter chord
    s.BetaIDq = [s.BetaIDq; BetaIDq];           % [NL x 1]  angle at duct quarter chord 
    s.CLd     = [s.CLd;         CLd];           % [NL x 1]  duct lift coefficient
    
    s.UARING  = [s.UARING;   UARING];           % [NL x Nd], induced velocity at duct
    s.URRING  = [s.URRING;   URRING];           % [NL x Nd], induced velocity at duct
    
    s.UADUCT  = [s.UADUCT;   UADUCT];           % [NL x Mp], induced velocity at propeller
    
    s.TAU     = [s.TAU; TAU];                   % [NL x 1]  thrust ratio 
    s.CTD     = [s.CTD; CTD];                   % [NL x 1]  CT for duct
end

% s.part4 '------ Performance metrics ------';
s.KT      = [s.KT;   KT];     % [NL x 1]  thrust coefficient
s.KQ      = [s.KQ;   KQ];     % [NL x 1]  torque coefficient
s.CT      = [s.CT;   CT];     % [NL x 1]  thrust coefficient
s.CQ      = [s.CQ;   CQ];     % [NL x 1]  torque coefficient
s.CP      = [s.CP;   CP];     % [NL x 1]  power  coefficient
s.EFFY    = [s.EFFY; EFFY];   % [NL x 1]  efficiency
s.VMWV    = [s.VMWV; VMWV];   % [NL x 1]  volumetric mean wake velocity

end

