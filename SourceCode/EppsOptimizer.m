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
% ================================================== EppsOptimizer Function
% Last Modified: 11/2/2011, Brenden Epps
% -------------------------------------------------------------------------
%
% Determines the "optimum" circulation, chord, and thickness distributions
% that satisfy the input operating conditions, using a variational 
% optimization algorithm for the propelller case (Coney, 1989), or a hybrid
% blade element momentum / vortex lattice method for the turbine case
% (Epps and Kimball, 2011).  
%
% Returns: "design" data structure with performance specs such as 
%          thrust coefficient and efficiency, as well as the optimized
%          circulation distribution, chord distribution, etc.
%
% -------------------------------------------------------------------------
%
% This implementation provides:
%   -- Propeller optimization via vortex-lattice theory, solved by...
%           ... linearized system of equations          (EppsOptimizer02.m) "LL-Linear"
%           ... Newton solver                           (EppsOptimizer23.m) "LL-Newton"
%           ... Newton solver with hub drag variation   (EppsOptimizer53.m) 
%   -- Turbine optimization via hybrid vortex lattice / momentum theory, solved by...
%           ... Newton solver                           (EppsOptimizer06.m) "Robust"
%
%   -- Chord length optimization via:
%           ... Specified maximum lift coefficient (default)
%           ... Cavitation mitigation method from (Epps, FAST2011)
%
%   -- Improved wake model for theoretical accuracy and numerical stability.
%   -- Improved duct model implementation for both propeller and turbine cases.
%
% -------------------------------------------------------------------------
          
function [design] = EppsOptimizer(input)
%%
% disp('Begin Epps optimizer'),    
disp(' '),            

% ================================================== Unpack input variables

% --------------------------------------------------------- Geometry inputs
Z = input.Z;           % number of blades

if isfield(input,   'Mp'), Mp    = input.Mp;    else  Mp    = 20;   end
if isfield(input,   'Vs'), Vs    = input.Vs;    else  Vs    = 1;    end % m/s

if     isfield(input,'D'),   R = input.D/2;  
elseif isfield(input,'R'),   R = input.R;
else                         R = 1;
end

if     isfield(input,'Dhub'),   Rhub = input.Dhub/2;
elseif isfield(input,'Rhub'),   Rhub = input.Rhub;
else                            Rhub = 0.2*R;
end

if isfield(input,'Rcirc'), Rcirc = input.Rcirc; else  Rcirc = Rhub; end % m
if isfield(input,'Rroot'), Rroot = input.Rroot; else  Rroot = Rhub; end % m

if Rcirc < Rhub,  disp('ERROR: Rcirc must be >= Rhub.  Setting Rcirc = Rhub.' ), Rcirc = Rhub;  end
if Rroot < Rcirc, disp('ERROR: Rroot must be >= Rcirc. Setting Rroot = Rcirc.'), Rroot = Rcirc; end
% -------------------------------------------------------------------------

% ----------- Blade geometry inputs: XR, XCoD, t0oc0, Xt0oD, XVA, XVT, dXVA
% If propeller geometry or inflow is not given, assume propeller 4118 form 
% with uniform axial inflow, no swirl inflow, and section CD == 0.008.
if isfield(input,'XR'), 
    XR  = input.XR;
    X1  =  ones(size(XR));
    X0  = zeros(size(XR));
    
    if isfield(input,'XCoD'), 
            XCoD  = input.XCoD;  
                      
        if     isfield(input,'Xt0oD'),  Xt0oD = input.Xt0oD;          
        elseif isfield(input,'t0oc0'),  Xt0oD = input.t0oc0 .* XCoD; 
        else                            Xt0oD = 0.1         .* XCoD;      
        end
                                  
    else
        XXR    = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];    % r/R
        XXCoD  = [0.1600 0.1818 0.2024 0.2196 0.2305 0.2311 0.2173 0.1806 0.1387 0.0010];
        Xt0oc0 = [0.2056 0.1551 0.1181 0.0902 0.0694 0.0541 0.0419 0.0332 0.0324 0.0000];
              
        XCoD   = interp1(XXR,XXCoD ,XR,'pchip','extrap');
        t0oc0  = interp1(XXR,Xt0oc0,XR,'pchip','extrap');
             
        if  isfield(input,'Xt0oD'),  
            Xt0oD = input.Xt0oD;
        else
            Xt0oD = t0oc0 .* XCoD;   
        end
    end 
   
    if isfield(input, 'XVA'  ),   XVA = input.XVA;   else   XVA = X1; end
    if isfield(input, 'XVT'  ),   XVT = input.XVT;   else   XVT = X0; end
    if isfield(input,'dXVA'  ),  dXVA = input.dXVA;  else  dXVA = X0; end
else 
    XR    = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];    % r/R
    X1    =  ones(size(XR));
    X0    = zeros(size(XR));
    XCoD  = [0.1600 0.1818 0.2024 0.2196 0.2305 0.2311 0.2173 0.1806 0.1387 0.0010];          % c/D
    t0oc0 = [0.2056 0.1551 0.1181 0.0902 0.0694 0.0541 0.0419 0.0332 0.0324 0.0000];          % t0/c  
    Xt0oD = t0oc0 .* XCoD;              
    XVA   = X1;                                      % Va/Vs
    XVT   = X0;                                      % Vt/Vs
   dXVA   = X0;
end
% -------------------------------------------------------------------------

% ------------------------------------------------ Inflow velocity profiles
if isfield(input,'ri')
    RI  = input.ri/R;

    if isfield(input,'VAI' ),   VAI = input.VAI;   else   VAI =  ones(size(RI)); end  % Va/Vs
    if isfield(input,'VTI' ),   VTI = input.VTI;   else   VTI = zeros(size(RI)); end  % Vt/Vs  
    if isfield(input,'dVAI'),  dVAI = input.dVAI;  else  dVAI = zeros(size(RI)); end  % delta(Va)/Vs       
else
     RI =  XR; %  r/R
    VAI =  XVA;
    VTI =  XVT;
   dVAI = dXVA;    
end

% Smooth the inflow velocity profiles:
VAI = RepairSpline(RI,VAI);
VTI = RepairSpline(RI,VTI);
% -------------------------------------------------------------------------

% ------------------------------------------------ Blade section properties
if isfield(input,'XCD'), XCD = input.XCD; 
    if length(XCD) == 1, XCD =   XCD*X1;  end
else                     XCD = 0.008*X1;    
end 
% -------------------------------------------------------------------------

% ---------------------------------------------------- Computational inputs
if isfield(input,'ITER'), ITER = input.ITER;  else  ITER = 50;   end
if isfield(input,'HUF' ), HUF  = input.HUF;   else  HUF  = 0;    end
if isfield(input,'TUF' ), TUF  = input.TUF;   else  TUF  = 0;    end
if isfield(input,'Rhv' ), Rhv  = input.Rhv;   else  Rhv  = 0.5;  end
% -------------------------------------------------------------------------

% ------------------------------------------------------- Cavitation inputs
if isfield(input,'rho'),  rho = input.rho;  else rho  = 1000;  end % kg/m^3
if isfield(input,'dVs'),  dVs = input.dVs;  else dVs  = 0.30;  end % m/s
if isfield(input,'H'  ),  H   = input.H;    else H    = 3.048; end % m
if isfield(input,'g'  ),  g   = input.g;    else g    = 9.81;  end % m/s^2
if isfield(input,'Patm'), Patm= input.Patm; else Patm = 101325;end % Pa
if isfield(input,'Pv'),   Pv  = input.Pv;   else Pv   = 2500;  end % Pa

SIGMAs  = (Patm + rho*g*H - Pv)/(0.5*rho*Vs^2);
% -------------------------------------------------------------------------

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

% 0 == do not optimize chord lengths, 1 == optimize chord lengths
if isfield(input,'Chord_flag'), Chord_flag = input.Chord_flag;  
                          else  Chord_flag = 0;  end
                          
% Method of chord optimization
if isfield(input,'ChordMethod'), ChordMethod = input.ChordMethod;  
                          else   ChordMethod = 'CLmax';  end
                          

% 0 == do not display plots, 1 == display plots
if isfield(input,'Plot_flag'),  Plot_flag = input.Plot_flag;  
                         else   Plot_flag = 0;  end 

% 0 == Horseshoe(...,Wrench(...)) analytic formulae, 1 == Wake_Horseshoe(...) numerical model possibly including wake roll-up     
if isfield(input,'Wake_flag'),  Wake_flag = input.Wake_flag;  
                         else   Wake_flag = 0;  end
% -------------------------------------------------------------------------
  
% ------------------------------------------------------ Tip speed ratio, L
if Propeller_flag == 1
    Js    = input.Js;
    L     = pi/Js;
else
    L     = input.L;
    Js    = pi/L;
end

D = 2*R;          % [m]
n = Vs/(Js*D);    % [rev/s]
N = 60*n;         % [RPM]
% -------------------------------------------------------------------------

% --------------------------------------------------------------- Duct_flag
if Duct_flag == 1
    if     isfield(input,'TAU'),      TAU   = input.TAU; else TAU = 1; end  % thrust ratio == propeller thrust / total thrust 
    
    if     isfield(input,'Rduct_oR'), Rduct_oR = input.Rduct_oR;  % duct radius
    elseif isfield(input,'Rduct'),    Rduct_oR = input.Rduct /R; 
    else                              Rduct_oR = 1;             
    end
    if     isfield(input,'Cduct_oR'), Cduct_oR = input.Cduct_oR;  % duct chord length
    elseif isfield(input,'Cduct'),    Cduct_oR = input.Cduct /R; 
    else                              Cduct_oR = 1;             
    end
    if     isfield(input,'Xduct_oR'), Xduct_oR = input.Xduct_oR;
    elseif isfield(input,'Xduct'),    Xduct_oR = input.Xduct /R;  % duct axial position downstream
    else                              Xduct_oR = 0;               % i.e. duct mid-chord centered at propeller          
    end
    if     isfield(input,'CDd'),      CDd   = input.CDd; else CDd = 0.008; end  % duct drag coefficient
    
else
    TAU      = 1;  % assume given CTdes is for the propeller alone, not propeller + duct
    Rduct_oR = 1;
end
% -------------------------------------------------------------------------

% ------------------------------------------------------ Thrust requirement
if Propeller_flag == 1
    
    % CTdes == desired total thrust coefficient == CTPdes + CTDdes
    if     isfield(input,'THRUST'),
            CTdes = input.THRUST / (0.5*rho*Vs^2*pi*R^2);
             
    elseif isfield(input,'CTDES'), 
            CTdes = input.CTDES; 
            
    elseif isfield(input,'CT'), 
            CTdes = input.CT; 
            
    elseif isfield(input,'KTDES'),
            CTdes = input.KTDES * (8/pi)/Js^2;
            
    elseif isfield(input,'KT'),
            CTdes = input.KT    * (8/pi)/Js^2;
    else 
            CTdes = 0;
    end
    
    % CQdes == desired torque coefficient
    if     isfield(input,'TORQUE'),
            CQdes = input.TORQUE / (0.5*rho*Vs^2*pi*R^3);
             
    elseif isfield(input,'CQDES'), 
            CQdes = input.CQDES; 
            
    elseif isfield(input,'CQ'), 
            CQdes = input.CQ; 
            
    elseif isfield(input,'KQDES'),
            CQdes = input.KQDES * (16/pi)/Js^2;
            
    elseif isfield(input,'KQ'),
            CQdes = input.KQ    * (16/pi)/Js^2;
    else 
            CQdes = 0;
    end
    
    
    if CTdes > 0 
        TorqueSpec_flag = 0;  % thrust is specified
    else
        TorqueSpec_flag = 1;  % torque is specified
    end

    % ------------------------------- Allocate CTdes between propeller and duct
    CTPdes = CTdes*   TAU;                       % CT desired for the propeller
    CTDdes = CTdes*(1-TAU);                      % CT desired for the duct    
    
else
  
    % CTDdes == desired duct thrust coefficient
    if     isfield(input,'THRUSTduct'),
           CTDdes = input.THRUSTduct / (0.5*rho*Vs^2*pi*R^2); 
 
    elseif isfield(input,'CTD'), 
           CTDdes = input.CTD;       
           
    elseif isfield(input,'KTD'),
           CTDdes = input.KTD   * (8/pi)/Js^2;
    else 
           CTDdes = 0;
    end

    TAU    = 0;  % (not used) assume given CTdes is for the duct alone, not turbine   + duct
    CTPdes = 0;  % (not used)
    
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------- Chord_flag
if Chord_flag == 1
   
    % ---------------------------------------------------------------------
    % Specified maximum allowable CL
    if      isfield(input,'XCLmax')    
        XCLmax = input.XCLmax; if length(XCLmax) == 1, XCLmax = XCLmax*X1;  end
    else
        XCLmax = 0.5 + (0.2-0.5)/(1-XR(1)) * (XR-XR(1));  % XCLmax == 0.5 at root and 0.2 at tip
    end
    % ---------------------------------------------------------------------
    
    % ---------------------------------------------------------------------
    if isfield(input,'CDoCL')
        disp(' '), disp('WARNING: input.CDoCL is not used.  Please give input.XCD and input.XCLmax.'), disp(' '), 
    end
    
    CDoCL = mean(abs(XCD./XCLmax));
    % ---------------------------------------------------------------------
    
    % ---------------------------------------------------------------------  
    % Specified Expanded Area Ratio (EAR)
    if  isfield(input,'EAR'), EARspec = input.EAR; else EARspec = 0; end
    % ---------------------------------------------------------------------    
    
    % ---------------------------------------------------------------------    
    if strcmp(ChordMethod,'FAST2011dCTP') == 1
        % -----------------------------------------------------------------
        % Method: see (Epps et al., FAST'2011)
        % -----------------------------------------------------------------        
        % --------------------------------------------
        % Specify "high-speed" condition
        if isfield(input,'Vh'), Vh = input.Vh;    else  Vh = Vs;   end % m/s
        
        Jh = Js;  % initial guess for advance coefficient for high-speed state
        
        % CTDES == desired total thrust coefficient
        if     isfield(input,'THRUSTh'),
                CTDESh = input.THRUSTh / (0.5*rho*Vh^2*pi*R^2);

        elseif isfield(input,'CTDESh'), 
                CTDESh = input.CTDESh; 

        elseif isfield(input,'CTh'), 
                CTDESh = input.CTh; 
        else
                CTDESh = CTdes;
        end        
        
        
        SIGMAh = (Patm + rho*g*H - Pv)/(0.5*rho*Vh^2);
        % --------------------------------------------       
    end
    % ---------------------------------------------------------------------    
    
else
    XCLmax = X1;
    CDoCL  = 0;
end

if Propeller_flag == 0
    XCLmax = -abs(XCLmax);
end
% -------------------------------------------------------------------------

% ------------------------------------------------------------ Viscous_flag 
if Viscous_flag == 0
    XCD   = X0;
    CDoCL = 0;
    CDd   = 0;
end
% -------------------------------------------------------------------------

% ------------- Geometry inputs for (Coney, 1989) chord optimization method
if isfield(input, 'Meanline'), Meanline  = input.Meanline;   else  Meanline  = 'NACA a=0.8 (modified)';       end
if isfield(input,'Thickness'), Thickness = input.Thickness;  else  Thickness = 'NACA 65A010'; end

if iscell(Meanline)

    if length(Meanline) ~= length(XR)
        disp('<ERROR>')
        disp('<ERROR> Meanline given as cell array but different length as XR...crash.')
        disp('<ERROR>') 
        return
    else
        Xf0octilde = 0*XR; % memory allocation
         XCLItilde = 0*XR; % memory allocation
         
        for j = 1:length(XR)
            [Xf0octilde(j), XCLItilde(j), junk1, junk2, junk3, junk4, junk5] = GeometryFoil2D(Meanline{j},Thickness{j});
        end
    end
else
    [f0octilde, CLItilde, junk1, junk2, junk3, junk4, junk5] = GeometryFoil2D(Meanline,Thickness);
    
    Xf0octilde = f0octilde * ones(size(XR));
     XCLItilde =  CLItilde * ones(size(XR));
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Check that turbine optimization is only attempted for uniform inflow
if (Propeller_flag == 0) && any(VAI < 0.99)
    disp('')
    disp('WARNING: OpenProp v3 turbine optimization is only valid for uniform inflow (Va/Vs == 1).')
    disp('')
end
% -------------------------------------------------------------------------

% =========================================================================
% =========================================================================
% =========================================================================



% =========================================================================
% =========================================================================
% ======================================================= Initialize design
Rhub_oR  = Rhub/R;    % [ ], hub                  radius / propeller radius
Rroot_oR = Rroot/R;   % [ ], blade root           radius / propeller radius
Rcirc_oR = Rcirc/R;   % [ ], circular section max radius / propeller radius

% ------------- Compute the Volumetric Mean Inflow Velocity, eqn 163, p.138
XRtemp   = linspace(Rhub_oR,1,100);         % (temp) radius / propeller radius
XVAtemp  = pchip(RI,VAI,XRtemp);            % (temp) axial inflow velocity /ship velocity
VMIV     = 2*trapz(XRtemp,XRtemp.*XVAtemp)/(1^2-Rhub_oR^2);  % [ ], VMIV/ship velocity

% ------------------------------------ Compute vortex & control point radii
[RC,RV,DR] = LLPanelRadii(Mp,Rhub_oR,Hub_flag,Duct_flag);


% ------------ Interpolate Va, Vt, Cd, and c/D at vortices & control points
 VAC = pchip(RI, VAI  ,RC);   %        axial inflow vel.  / Vs at ctrl pts
dVAC = pchip(RI,dVAI  ,RC);   % delta( axial inflow vel.) / Vs due to wake, for chord optimization
 VTC = pchip(RI, VTI  ,RC);   % tangential inflow vel. / ship vel. at ctrl pts
CD   = pchip(XR,XCD   ,RC);   % section drag coefficient           at ctrl pts
t0oD = pchip(XR,Xt0oD ,RC);   % section thickness / propeller dia. at ctrl pts
CLmax= pchip(XR,XCLmax,RC);   % maximum allowable lift coefficient at ctrl pts

f0octilde = pchip(XR,Xf0octilde,RC);
 CLItilde = pchip(XR, XCLItilde,RC);

if (abs(XR(end)-1) < 1e-4) && (XCoD(end) <= 0.01)  % if XR == 1 and XCoD == 0
    
    CoD  = InterpolateChord(XR,XCoD,RC);   % section chord / propeller diameter at ctrl pts
else
    CoD  =            pchip(XR,XCoD,RC);   % section chord / propeller diameter at ctrl pts
end

t0oc = t0oD./CoD; % input t0/c values




% ------------------------- Initial estimates based on Actuator Disc Theory
if Propeller_flag == 1                            
    
    G	   = zeros(Mp,1);			                %  G = Gamma / 2*pi*R*Vs

    UASTAR = 0*RC + 0.5*(sqrt(1+CTPdes)-1);         % Kerwin & Hadler (2010), eqn (4.26)
    UTSTAR = 0*RC;
    
    if (VMIV < 0.05)                       % assume bollard pull propeller design
        UASTAR = 0*RC + 0.5*sqrt(CTPdes);  % (Epps, 2013)
    end                    
else                                                %  element is a turbine
    
    [CPBetz,BetzRC,BetzG,BetzUA,BetzUT,BetzTAN, BetzGRC] = Turbine_ADS_Theory(L,Z,CDoCL,RC);
    
    G      = BetzGRC';                              % G = Gamma / 2*pi*R*Vs
           
    UASTAR = pchip(BetzRC,BetzUA,RC);
    UTSTAR = pchip(BetzRC,BetzUT,RC);
end

    TANBC  = (VAC         )./(L*RC + VTC);
    TANBIC = (VAC + UASTAR)./(L*RC + VTC + UTSTAR);

    VSTAR  = sqrt((VAC+UASTAR).^2 + (L*RC+VTC+UTSTAR).^2);        % V* / Vs

    dVdG   = zeros(Mp,Mp);					   	 % (2piR) * d(V*) /d(Gamma)
% -------------------------------------------------------------------------


% ------------ If optimizing chord, set chord very small on first iteration
if Chord_flag == 1
    if Propeller_flag == 1
        CoD = 0.01*ones(size(RC));
    else
        CoD = 2*pi*G'./(VSTAR.*CLmax);
    end 
end


% ------------------------- Initialize vortex Horseshoe Influence Functions 
if Wake_flag == 0 
    [UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBIC,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);
    
else
    URSTAR = zeros(1,Mp);                                         % URSTAR / Vs
    
    disp('ERROR: Wake_Geometry and Wake_Horseshoe are no longer supported...please re-implement these...')
    return
        % disp('Beginning to update {UAHIF,UTHIF,URHIF}'), 
        % [WX,WY,WZ]          =  Wake_Geometry(Mp,RC,RV,TANBIV,VAC,UASTAR,URSTAR,CTPdes);
        % epsilon = min(diff(RV))/100;   % vortex core size in Wake_Influence(...) function
        % %[UAHIF,UTHIF,URHIF] = Wake_Horseshoe(Mp,Z,RC,SCF,WX,WY,WZ,epsilon);
        % [UAHIF,UTHIF,URHIF] = Wake_Horseshoe(Mp,Z,TANBIV,RC,RV,SCF,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR,WX,WY,WZ,epsilon);
        % disp('Done updating {UAHIF,UTHIF,URHIF}'), 
end
% ------------------------------------------------------------------------- 

% -------------------------------------------------------------------------
%------------------------------------------------ Initialize duct variables
if Duct_flag == 1  % duct is present
    %     XdRING            [1,Nd], x/R location of each vortex ring downstream of propeller
    %     VARING            [1,Nd], (axial free-stream velocity at duct)/Vs
    %
    %     GdRING            [1,Nd],  fraction of total non-dimensional duct circulation, sucht that sum(GdRING) = 1 
    %                                (i.e. non-dimensional circulation per unit Gd)
    %     Gd                         total non-dimensional circulation about the duct, Gd == Gamma_d/(2*pi*R*Vs),
    %                                such that the non-dimensional circulation of ring n is Gd*GdRING(n).
    %     Gamma_d [m^2/s]            total     dimensional circulation about the duct [m^2/s]
    %
    %     CDd                       section drag coefficient for the duct
    %     CTDdes                    desired duct thrust coefficient
    %     CTD                       duct thrust coeff (viscous drag included) with total duct circulation of Gd 
    %
    % Influence of propeller on duct:
    %
    %   (DAHIF ,DRHIF)     [Nd,Mp], (axial/radial horseshoe influence functions of prop on duct)*(2*pi*R)
    %   (UARING,URRING)    [1,Nd],  (axial/radial velocity induced at duct (XdRING,Rduct_oR) by prop) / Vs
    %
    %
    % Influence of duct on propeller:
    %
    %   UARING      [1,Mp]  (axial velocity induced on PROP by duct)/Vs
    %   UADIF       [1,Mp]   axial velocity induced on PROP by duct per unit Gd, 
    %                        i.e. non-dimensional Duct Influence Function
    %                         == 2*pi*R * (duct influence function with unit dimensional duct circulation)
    %    
    [XdRING,GdRING,UADIF] = Duct_Influence(Rduct_oR,Cduct_oR,Xduct_oR,RC); % Duct geometry, and influence of duct on propeller

    % ---------------------------- (axial inflow velocity / Vs) at Rduct_oR
    VARING = pchip(RI,VAI,Rduct_oR);  
    
    % ------------------------------------------------ Initial guess for Gd
    Gd     = 0;                  % total duct circulation / (2*pi*R*Vs)
    UADUCT = UADIF * Gd;         % axial velocity induced at RC by the duct
    
    % ---------------------------------------------------------------------
    % Initialize Duct Horseshoe Influence Functions (influence of propeller on duct)
    %
    % DAHIF(n,m) = influence of m-th horseshoe vortex shed from propeller (Mp panels) 
    %                    on the n-th control point of the duct            (Nd rings) 
    %
    disp(' '), disp('Computing rotor-duct interaction...be patient...'), disp(' '),
    [DAHIF_times_TANBIC, DTHIF, DRHIF_times_TANBIC] = Horseshoe_intr_110830(XdRING,Rduct_oR, RC,ones(size(TANBIC)),RV,Z,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR); 

    DAHIF = 0* DAHIF_times_TANBIC; % memory allocation
    DRHIF = 0* DRHIF_times_TANBIC; % memory allocation
    
    for m = 1:Mp;                               
        DAHIF(:,m) = DAHIF_times_TANBIC(:,m) / TANBIC(m);
        DRHIF(:,m) = DRHIF_times_TANBIC(:,m) / TANBIC(m);
    end  
    
else
    Gd       = 1;  % if set Gd==0, then Gd_res would be Inf  
    UADUCT   = 0*RC;
    CTD      = 0;  
end
% -------------------------------------------------------------------------


% ------------------------------------------------ Implement RepairSpline.m
% To smooth data X(RC):  X_smooth = X*Bsmooth;
Bsmooth = RepairSplineMatrix(RC);
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------- Configure plot windows
if Plot_flag == 1
    % close all,
    set(0,'DefaultAxesFontSize',12)    

    % Rainbow Color matrix
    CLR = [     1       0       0;      ... % (1) Red
                1       0.5     0;      ... % (2) Orange
                0.8     0.8     0;      ... % (3) Burnt Yellow  
                0       0.8     0;      ... % (4) Green
                0       0       1;      ... % (5) Blue
                0.75    0       0.75;   ... % (6) Purple
                0.3     0.3     0.3];   ... % (7) Gray   

    color_count = 1;    
    
    % ------------------------------------------------- Uncomment to plot G
    Hgamma = fig;                  set(Hgamma,'Position',[130 680 560 420])
        HHgamma = plot(RC,G,'k.-');   
        plot([Rhub_oR,1],[0 0],'k')
        if Propeller_flag == 0, plot(RC,BetzGRC,'r--');  end 
        set(gca,'XLim',[Rhub_oR 1])
        xlabel('r/R','Fontsize',14,'Fontname','Times')
        ylabel('G' ,'Fontsize',14,'Fontname','Times')
    % ---------------------------------------------------------------------     

    % ------------------------------ Uncomment to plot UASTAR, UTSTAR, & URSTAR
    Hvel = fig;                      set(Hvel,'Position',[695 680 560 420])
        if Duct_flag == 1
               HHvel = plot(RC,UASTAR,'b.-',RC,UTSTAR,'r.-',RC,UADUCT,'g.-');
        else
               HHvel = plot(RC,UASTAR,'b.-',RC,UTSTAR,'r.-');
        end
        plot([Rhub_oR,1],[0 0],'k')
        if Propeller_flag == 0, plot(BetzRC,BetzUA,'b--',BetzRC,BetzUT,'r--');  end
        axis([Rhub_oR 1 -1 1])
        xlabel('r/R'    ,'Fontsize',14,'Fontname','Times')
        ylabel('UASTAR (blue), UTSTAR (red)','Fontsize',14,'Fontname','Times')
    % --------------------------------------------------------------------- 

    % ------------------------------------ Uncomment to plot TANBC & TANBIC
    Hbeta = fig;                    set(Hbeta,'Position',[695 160 560 420])
                 plot(RC,TANBC,'k.-');
        HHbeta = plot(RC,TANBIC,'r.-'); 
        plot([Rhub_oR,1],[0 0],'k')
        if Propeller_flag == 0, plot(BetzRC,BetzTAN,'r--');  end
        set(gca,'XLim',[Rhub_oR 1])
        xlabel('r/R'    ,'Fontsize',14,'Fontname','Times')
        ylabel('TANBIC (red), TANBC (black)','Fontsize',14,'Fontname','Times')
    % ---------------------------------------------------------------------    

    % ----------------------------------------------- Uncomment to plot CoD
    Hcod = fig;                      set(Hcod,'Position',[130 160 560 420])
        XXR   = linspace(Rhub_oR,1,100);
                % -------------------------------------------------------------------------        
                % ------------------------------------------------------- Interpolate chord
                if Chord_flag == 0  % if chord optimization WAS NOT performed, 
                                    % then interpolate chord and thickness the same as in EppsOptimizer.m

                    if  (abs(XR(end)-1) < 1e-4) && (XCoD(end) <= 0.01)  % if XR == 1 and XCoD == 0

                        XXCoD  = InterpolateChord(XR,XCoD, XXR);   % section chord / propeller diameter at ctrl pts
                    else
                        XXCoD  =            pchip(XR,XCoD, XXR);   % section chord / propeller diameter at ctrl pts
                    end
                else % chord optimization WAS performed, so we need to be careful to iterpolate in such a way as to  
                     % preserve a   zero-chord-length tip (as in the case of a free-tip propeller) or 
                     % preserve a finite-chord-length tip (as in the case of a ducted   propeller).

                    if (Duct_flag == 0) || (Rduct_oR > 1.001)   

                        XXCoD  = InterpolateChord(RC,CoD ,XXR);  % yields zero   chord at the tip       
                    else                                                
                        XXCoD  =            pchip(RC,CoD ,XXR);  % yields finite chord at the tip   
                    end     
                end
                % -------------------------------------------------------------------------        
                % -------------------------------------------------------------------------        
        HHcod = plot(XXR,-XXCoD,'k-',XXR,XXCoD,'k-',RC,-CoD,'b.',RC,CoD,'b.');
        set(gca,'XLim',[Rhub_oR 1]),  % axis equal, 
        xlabel('r/R'    ,'Fontsize',14,'Fontname','Times')
        ylabel('expanded blade  (c/R)','Fontsize',14,'Fontname','Times')
    % ---------------------------------------------------------------------

    pause(0.001), drawnow,
else
    Hvel = 0;
end  % if Plot_flag == 1
% -------------------------------------------------------------------------


%%
% =========================================================================
% =========================================================================
disp('--------- Begin circulation optimization'),
disp(' '),

% ------------------------- Propeller optimization method (choose one only)
if isfield(input,'EppsOptimizer02_flag'), EppsOptimizer02_flag = input.EppsOptimizer02_flag; else EppsOptimizer02_flag = 1; end  % Propeller: LL-Linear
if isfield(input,'EppsOptimizer23_flag'), EppsOptimizer23_flag = input.EppsOptimizer23_flag; else EppsOptimizer23_flag = 0; end  % Propeller: LL-Newton (standard    hub drag model)
if isfield(input,'EppsOptimizer53_flag'), EppsOptimizer53_flag = input.EppsOptimizer53_flag; else EppsOptimizer53_flag = 0; end  % Propeller: LL-Newton (variational hub drag model)
if (EppsOptimizer23_flag == 1 || EppsOptimizer53_flag == 1), EppsOptimizer02_flag = 0; end

% ------------------------------ Initialize optimization routine parameters
LM      = -1;                       % Lagrange Multiplier (propeller case)
LM_last = LM;                       % last value of the Lagrange Multiplier
 G_last = 0*G;                      % last (previous) value of G
Gd_last = 0;                        % last (previous) value of Gd
 G_iter = 1;                        % iteration in the G loop
 G_res  = 1;                        % residual G  between iterations
Gd_res  = 0;                        % residual Gd between iterations
 C_res  = 0;                        % residual chord length during optimization  
relax   = 0.9;                      % Newton solver relaxation parameter
 G_TOL  = 1e-4;                     % tolerance on the circulation for convergence of the optimizer

if (Chord_flag == 1)  &&  ( (strcmp(ChordMethod,'FAST2011dCTP') == 1) || (strcmp(ChordMethod,'FAST2011dVAC') == 1) || (strcmp(ChordMethod,'Brizzolara2007') == 1) )
    ITER = ITER*3;  % allow three circulation iterations per chord iteration
end

%% 
while G_iter <= ITER && any([G_res;Gd_res] > G_TOL)        % (WHILE LOOP G1)
%     disp(['----- Gamma iteration: ',num2str(G_iter)])   % status message

    % =@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=
    % UPDATE: G, UASTAR, UTSTAR, TANBIC, (duct stuff), LM
    % --------------------------------------------------------------------- 
    if Propeller_flag == 1 % optimize a propeller
        
        
        % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-        
        % "LL-Linear" (EppsOptimizer02.m, which is similar to the v2.4.1 optimizer)       
        % ----------------------------------------------------------------- 
        if EppsOptimizer02_flag == 1
            % ---------------------- Set up simultaneous equations for G and LM
            A  = zeros(Mp+1,Mp+1);    % A matrix for linear system of equations
            B  = zeros(Mp+1,1);       % B matrix for linear system of equations

            for i = 1:Mp                           % for each equation for G(i)
                for m = 1:Mp                       % for each vortex panel, m        
                    A(i,m) =  UAHIF(m,i)*RC(m)*DR(m)    ...                
                            + UAHIF(i,m)*RC(i)*DR(i)    ...
                            + LM_last*UTHIF(m,i)*DR(m)  ...
                            + LM_last*UTHIF(i,m)*DR(i); 
                end   

                % 11/16/2011 BEPPS: This was the implementation in EppsOptimizer02.m, but the inclusion
                %       of these drag terms has negligable effect on the converged circulation distribution.
                %
                % B(i)  = -(VAC(i) + UADUCT(i))*RC(i)*DR(i); ...
                %         -(1/(2*pi))*sum(CD .* dVdG(:,i)'.* CoD .* (L*RC + VTC + UTSTAR).*RC.*DR) ...
                %         -(1/(2*pi))*sum(CD .* VSTAR     .* CoD .* UTHIF(:,i)'          .*RC.*DR);  
                %
                % A(i,Mp+1) =  (L*RC(i) + VTC(i))*DR(i);   ...
                %            - (1/(2*pi))*sum(CD .* dVdG(:,i)' .* CoD .* (VAC + UADUCT + UASTAR).*DR) ...
                %            - (1/(2*pi))*sum(CD .* VSTAR      .* CoD .*        UAHIF(:,i)'   .*DR);   
                %
                % ------
                B(i)      = -(VAC(i) + UADUCT(i))*RC(i)*DR(i);

                A(i,Mp+1) =  (L*RC(i) + VTC(i))*DR(i);
                % ------
            end

            % The (Mp+1) equation is either the thrust constraint or torque constraint
            % Note: A(Mp+1,Mp+1) = 0            
            if TorqueSpec_flag == 0  % thrust is specified, CT_prop_inviscid/(4*Z) == CTPdes/(4*Z) + CT_prop_viscous/(4*Z) + CT_hub/(4*Z)
                
                for m = 1:Mp                          
                    A(Mp+1,m) = (L*RC(m) + VTC(m) + UTSTAR(m))*DR(m);  
                end

                B(Mp+1) =  CTPdes/(4*Z)  +  (1/(2*pi))*sum(CD.*VSTAR.*CoD.*(VAC + UADUCT + UASTAR).*DR);

                if Hub_flag == 1
                    B(Mp+1) = B(Mp+1) + (Z/8)*(log(1/Rhv)+3)*(G_last(1)^2);
                end
                
            else % torque is specified, CQ_prop_inviscid/(4*Z) == CQdes/(4*Z) - CQ_prop_viscous/(4*Z)
                
                for m = 1:Mp
                    A(Mp+1,m) = (VAC(m) + UADUCT(m) + UASTAR(m))*RC(m)*DR(m);
                end

                B(Mp+1) =  CQdes/(4*Z)  -  (1/(2*pi))*sum(CD.*VSTAR.*CoD.*(L*RC + VTC + UTSTAR).*RC.*DR);
            end
                
            % -------------------------------- Solve linear system of equations
            GLM = linsolve(A,B);             
            G   = GLM(1:Mp);                     % G is the first Mp entries
            LM  = GLM(Mp+1);                     % LM is the last entry


            % ----------------------------------------------------------------- 
            % Update induced velocities (influence of propeller on propeller)
            UASTAR = (UAHIF*G)';  
            UTSTAR = (UTHIF*G)';  

            TANBIC = (VAC + UADUCT + UASTAR)./(L*RC + VTC + UTSTAR);

            % Smooth the inflow angle for numerical stability:
            TANBICsmooth = TANBIC * Bsmooth;
            % -----------------------------------------------------------------

            if any(isnan(GLM))  || ~isreal(GLM) || max(G-G_last) > 10
                G      = 0*RC';
                UASTAR = 0*RC;
                UTSTAR = 0*RC;
                TANBIC = TANBC;

                disp(' ')
                disp('<WARNING>')
                disp('<WARNING> GLM == NaN or imaginary... crash avoided...')
                disp('<WARNING>')
                disp(' ')
                disp('Switching numerical method to Newton solver...')
                disp(' ')

                EppsOptimizer02_flag = 0;
                EppsOptimizer23_flag = 1;
            end   
            % -----------------------------------------------------------------      

            % ----------------------------------------------------------------- 
            if Duct_flag == 1
                % Update the influence of propeller on duct
                for m = 1:Mp;
                    DAHIF(:,m) = DAHIF_times_TANBIC(:,m) / TANBICsmooth(m);
                    DRHIF(:,m) = DRHIF_times_TANBIC(:,m) / TANBICsmooth(m);
                end

                % Update induced velocities at the duct (influence of propeller on duct)
                UARING = (DAHIF*G)';  
                URRING = (DRHIF*G)'; 
                
                % If propeller torque specified, then need to update desired duct thrust based on current prop thrust and specified thrust ratio, TAU
                if TorqueSpec_flag == 1
                                        
                    CTP = 4*Z*sum(  G'.*(L*RC + VTC + UTSTAR).*DR  -  (1/(2*pi)).*VSTAR.*CoD.*CD.*(VAC + UADUCT + UASTAR).*DR  );
              
                    if Hub_flag == 1
                        CTP = CTP - 0.5*(log(1/Rhv)+3)*(Z*G(1))^2;   
                    end

                    CTDdes = ( (1-TAU)/TAU ) * CTP;
                end

                % Update duct circulation such that CTD == CTDdes            
                [junk,Gd] = Duct_Thrust(XdRING,Rduct_oR,VARING,UARING,URRING,GdRING,Gd,CDd,CTDdes);

                % Update the induced velocities at the propeller (influence of duct on propeller)
                UADUCT =  UADIF*Gd;                        
            end
            % ----------------------------------------------------------------- 
        
        end % if EppsOptimizer02_flag == 1
        % ----------------------------------------------------------------- 
        % END: "LL-Linear"
        % -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

       
        % -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
        % "LL-Newton" (EppsOptimizer23.m)
        % -------------------------------------------------------------------------
        if (EppsOptimizer23_flag == 1 || EppsOptimizer53_flag == 1)
            % ------------------------------------------ Execute Newton solver   
            RNS = zeros(4*Mp+1,     1);     % Residual vector
            JNS = zeros(4*Mp+1,4*Mp+1);     % Jacobian

            % Evaluate induced velocities
            UASTARtemp = (UAHIF*G)';  
            UTSTARtemp = (UTHIF*G)';

            % ---------------------------------------------- Evaluate residuals        
            for i = 1:Mp

                % 11/16/2011 BEPPS: This was the implementation in EppsOptimizer23.m, but the inclusion
                %       of these drag terms has negligable effect ( O(10^-4) ) on the converged circulation distribution.
                %            
                % RNS(i) =       (VAC(i) + UADUCT(i) + UASTAR(i))     *RC(i)*DR(i)  ...
                %          + sum( UAHIF(:,i)'.*G'.*RC  .*DR   ) ...
                %          + sum( (1/(2*pi))*CD.*CoD .* dVdG(:,i)' .* (L*RC + VTC + UTSTAR) .* RC .* DR ) ...
                %          + sum( (1/(2*pi))*CD.*CoD .*  VSTAR     .* (      UTHIF(:,i)'      ) .* RC .* DR ) ...
                %          ...    
                %          + LM * ( ...
                %                        (L*RC(i) + VTC(i) + UTSTAR(i)) * DR(i)  ...
                %                  + sum( UTHIF(:,i)' .* G'                .* DR   ) ...
                %                  - sum( (1/(2*pi))*CD.*CoD .* dVdG(:,i)' .* (  VAC + UADUCT + UASTAR     ) .* DR ) ...
                %                  - sum( (1/(2*pi))*CD.*CoD .*  VSTAR     .* (                 UAHIF(:,i)') .* DR ) ...
                %                 );
                %
                % ------
                RNS(i) =  (VAC(i) + UADUCT(i) + UASTAR(i)) *RC(i) *DR(i)  ...
                         +        sum( UAHIF(:,i)' .* G'  .*RC   .*DR   ) ...
                         ...
                         + LM * ( sum( UTHIF(:,i)' .* G'        .* DR   ) ...
                                 +(L*RC(i) + VTC(i) + UTSTAR(i)) * DR(i)  ...
                                );            
                % ------

                RNS(i+  Mp) = UASTAR(i) - UASTARtemp(i);
                RNS(i+2*Mp) = UTSTAR(i) - UTSTARtemp(i);
                RNS(i+3*Mp) = TANBIC(i) - (VAC(i) + UADUCT(i) + UASTAR(i))/(L*RC(i) + VTC(i) + UTSTAR(i)); 
            end    

            % The (Mp+1) equation is either the thrust constraint or torque constraint
            %
            if TorqueSpec_flag == 0  % thrust is specified
                
                    RNS(1+4*Mp) = sum((L*RC + VTC + UTSTAR).*G'.*DR - (1/(2*pi))*CD.*CoD.*VSTAR.*(VAC + UADUCT + UASTAR).*DR)  - CTPdes/(4*Z);                

                if     Hub_flag == 1 && EppsOptimizer23_flag == 1

                    RNS(1+4*Mp) = RNS(1+4*Mp) - (Z/8)*(log(1/Rhv)+3)*(G_last(1)^2);  % EppsOptimizer23.m hub drag treatment  (do NOT include hub drag in variational optimization)

                elseif Hub_flag == 1 && EppsOptimizer53_flag == 1

                    RNS(1)      = RNS(1)      - (Z/4)*(log(1/Rhv)+3)*(G(1)*LM);      % EppsOptimizer53.m hub drag treatment  (include hub drag in variational optimization) 

                    RNS(1+4*Mp) = RNS(1+4*Mp) - (Z/8)*(log(1/Rhv)+3)*(G(1)^2);       % EppsOptimizer53.m hub drag treatment  (include hub drag in variational optimization)   
                end
                
            else  % torque is specified
                
                    RNS(1+4*Mp) = sum( (VAC + UADUCT + UASTAR).*G'.*RC.*DR + (1/(2*pi))*CD.*CoD.*VSTAR.*(L*RC + VTC + UTSTAR).*RC.*DR )  - CQdes/(4*Z);
            end

            % ----------------------------------------------- Evaluate Jacobian        
            for i = 1:Mp 
                JNS(i     ,(1:Mp)     ) = UAHIF(:,i)' .* RC .* DR  +  LM * UTHIF(:,i)' .* DR;

                JNS(i     , i    +  Mp) = JNS(i,i+  Mp) + RC(i)*DR(i);
                JNS(i     , i    +2*Mp) = JNS(i,i+2*Mp) + LM   *DR(i);

                % 11/16/2011 BEPPS: This was the implementation in EppsOptimizer23.m, but the inclusion
                %       of these drag terms has negligable effect ( O(10^-4) ) on the converged circulation distribution.
                %
                % JNS(i     ,(1:Mp)+  Mp) = - LM * (1/(2*pi))*CD.*CoD.*dVdG(:,i)'    .*DR;  
                % JNS(i     ,(1:Mp)+2*Mp) =        (1/(2*pi))*CD.*CoD.*dVdG(:,i)'.*RC.*DR;  
                %
                % JNS(i ,1+4*Mp)      =       (L*RC(i) + VTC(i) + UTSTAR(i)) * DR(i)  ...
                %                       + sum( UTHIF(:,i)' .* G'                .* DR   ) ...
                %                       - sum( (1/(2*pi))*CD.*CoD .* dVdG(:,i)' .* (  VAC + UADUCT + UASTAR     ) .* DR ) ...
                %                       - sum( (1/(2*pi))*CD.*CoD .*  VSTAR     .* (        UAHIF(:,i)') .* DR );
                %
                % ------
                JNS(i ,1+4*Mp)      =       (L*RC(i) + VTC(i) + UTSTAR(i)) * DR(i)  ...
                                      + sum( UTHIF(:,i)' .* G'            .* DR   );
                % ------


                JNS(i+  Mp,1:Mp  ) = - UAHIF(i,1:Mp);
                JNS(i+  Mp,i+  Mp) = 1;

                JNS(i+2*Mp,1:Mp  ) = - UTHIF(i,1:Mp);
                JNS(i+2*Mp,i+2*Mp) = 1;

                JNS(i+3*Mp,i+  Mp) = -1                          /(L*RC(i) + VTC(i) + UTSTAR(i));
                JNS(i+3*Mp,i+2*Mp) = (VAC(i)+UADUCT(i)+UASTAR(i))/(L*RC(i) + VTC(i) + UTSTAR(i))^2;
                JNS(i+3*Mp,i+3*Mp) = 1;

                
                % The (Mp+1) equation is either the thrust constraint or torque constraint
                %
                if TorqueSpec_flag == 0  % thrust is specified
                
                    JNS(1+4*Mp,i     ) = (L*RC(i) + VTC(i) + UTSTAR(i))*DR(i);
                    JNS(1+4*Mp,i+  Mp) = - (1/(2*pi))*CD(i)*CoD(i)*VSTAR(i)*DR(i);    
                    JNS(1+4*Mp,i+2*Mp) = G(i)*DR(i);


                    if Hub_flag == 1 && EppsOptimizer53_flag == 1     
                        JNS(1,1)      = JNS(1,1)      - (Z/4)*(log(1/Rhv)+3)*LM;   % EppsOptimizer53.m hub drag treatment     
                        JNS(1,1+4*Mp) = JNS(1,1+4*Mp) - (Z/4)*(log(1/Rhv)+3)*G(1); % EppsOptimizer53.m hub drag treatment 
                        JNS(1+4*Mp,1) = JNS(1+4*Mp,1) - (Z/4)*(log(1/Rhv)+3)*G(1); % EppsOptimizer53.m hub drag treatment
                    end            

                else  % torque is specified

                    JNS(1+4*Mp,i     ) = (VAC(i) + UADUCT(i) + UASTAR(i))*RC(i)*DR(i);
                    JNS(1+4*Mp,i+  Mp) =                            G(i) *RC(i)*DR(i);
                    JNS(1+4*Mp,i+2*Mp) = (1/(2*pi))*CD(i)*CoD(i)*VSTAR(i)*RC(i)*DR(i);             
                end
                
            end

            % ----------------------------- Update Newton solver vector of unknowns
            DX     = linsolve(JNS,-RNS);
            G      = G      + relax*DX( 1:Mp      ) ;
            UASTAR = UASTAR + relax*DX((1:Mp)+  Mp)';
            UTSTAR = UTSTAR + relax*DX((1:Mp)+2*Mp)';
            TANBIC = TANBIC + relax*DX((1:Mp)+3*Mp)';
            LM     = LM     + relax*DX(     1+4*Mp) ;

            % Smooth the inflow angle for numerical stability:
            TANBICsmooth = TANBIC * Bsmooth;
            % -----------------------------------------------------------------  

            if any(isnan(DX))  || ~isreal(DX)  || any(DX > 999)
                 G      = 0*RC';
                 UASTAR = 0*RC;
                 UTSTAR = 0*RC;
                 TANBIC = TANBC;

                 disp(' ')
                 disp('<WARNING>')
                 disp('<WARNING> DX == NaN or imaginary... crash avoided...')
                 disp('<WARNING>')
                 disp(' ')
                 G_iter = 999;
            end
            % ----------------------------------------------------------------- 

            % -----------------------------------------------------------------
            if Duct_flag == 1
                % Update the influence of propeller on duct
                for m = 1:Mp;
                    DAHIF(:,m) = DAHIF_times_TANBIC(:,m) / TANBICsmooth(m);
                    DRHIF(:,m) = DRHIF_times_TANBIC(:,m) / TANBICsmooth(m);
                end

                % Update induced velocities at the duct (influence of propeller on duct)
                UARING = (DAHIF*G)';  
                URRING = (DRHIF*G)'; 

                
                % If propeller torque specified, then need to update desired duct thrust based on current prop thrust and specified thrust ratio, TAU
                if TorqueSpec_flag == 1
                                        
                    CTP = 4*Z*sum(  G'.*(L*RC + VTC + UTSTAR).*DR  -  (1/(2*pi)).*VSTAR.*CoD.*CD.*(VAC + UADUCT + UASTAR).*DR  );
              
                    if Hub_flag == 1
                        CTP = CTP - 0.5*(log(1/Rhv)+3)*(Z*G(1))^2;   
                    end

                    CTDdes = ( (1-TAU)/TAU ) * CTP;
                end
                
                % Update duct circulation such that CTD == CTDdes            
                [junk,Gd] = Duct_Thrust(XdRING,Rduct_oR,VARING,UARING,URRING,GdRING,Gd,CDd,CTDdes);

                % Update the induced velocities at the propeller (influence of duct on propeller)
                UADUCT =  UADIF*Gd;                        
            end
            % ----------------------------------------------------------------- 
        
        end % if (EppsOptimizer23_flag == 1 || EppsOptimizer53_flag == 1)
        % -----------------------------------------------------------------
        % END: "LL-Newton"
        % -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


    else % optimize a turbine
  
        % ----------------------------------------------------------------- 
        %               To seek the sacred river Alph                     %
        %               To walk the caves of ice                          %
        %               To break my fast on honeydew                      %
        %               And drink the milk of Paradise...                 %      
        % -----------------------------------------------------------------  
       
 

        % -!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-       
        % -----------------------------------------------------------------        
        RNS = zeros(4*Mp,1);     % Residual vector
        JNS = zeros(4*Mp,4*Mp);  % Jacobian

        % Evaluate induced velocities
        UASTARtemp = (UAHIF*G)';  
        UTSTARtemp = (UTHIF*G)';

        % ---------------------------------------------- Evaluate residuals        
        for i = 1:Mp
             
            if (Chord_flag == 1) && (strcmp(ChordMethod,'CLmax') == 1)
                % ==========
                % 11/16/2011 BEPPS: This new drag treatment is consistent with actuator disk theory.
                %                   The drag term here assumes d(CL)/d(G) == 0, which is true in chord optimization with CL == CLmax.
                % 11/17/2011 BEPPS: However, in the no-chord-optimization case, d(CL)/d(G) is not zero, so the drag term is zero to the leading order.  
                %                   Furthermore, this code as is crashes for the no-chord-optimization case, so use the formulation below instead.
                RNS(i) =   ( VAC(i) + UADUCT(i) + 2*UASTAR(i)) * ( VAC(i) + UADUCT(i) +   UASTAR(i))             ...    
                         -                                       (L*RC(i) +    VTC(i) + 2*UTSTAR(i)) * UTSTAR(i) ...
                         + ( VAC(i) + UADUCT(i) + 2*UASTAR(i)) * (L*RC(i) +    VTC(i) + 2*UTSTAR(i)) * CD(i)/CLmax(i);
                % ==========
            else
                % 11/16/2011 BEPPS: This was the "Robust method" implementation in EppsOptimizer06.m, but the inclusion
                %                   of these drag terms has negligable effect ( O(10^-4) ) on the converged circulation distribution.
                %            
                % Note, the first line is equivalent to:
                %
                %             (VAC(i) + UADUCT(i))^2 +  (3*(VAC(i) + UADUCT(i))+ 2*UASTAR(i))*UASTAR(i) ...       
                %
                % RNS(i) =   ( VAC(i) + UADUCT(i) + 2*UASTAR(i)) * (VAC(i) + UADUCT(i) + UASTAR(i))  ...    
                %          - (L*RC(i) +    VTC(i) + 2*UTSTAR(i)) * UTSTAR(i) ... 
                %          + ( VAC(i) + UADUCT(i) + 2*UASTAR(i)) * (1/(2*pi))*CD(i)*CoD(i)*dVdG(i,i)*(L*RC(i) + VTC(i) + UTSTAR(i)) ...
                %          + ( VAC(i) + UADUCT(i) + 2*UASTAR(i)) * (1/(2*pi))*CD(i)*CoD(i)*VSTAR(i) * UTHIF(i,i);  
                % ------  
                % Inviscid terms:    
                RNS(i) =   ( VAC(i) + UADUCT(i) + 2*UASTAR(i)) * (VAC(i) + UADUCT(i) + UASTAR(i))  ...    
                         - (L*RC(i) +    VTC(i) + 2*UTSTAR(i)) * UTSTAR(i);                    
                % ------   
            end

            RNS(i+  Mp) = UASTAR(i) - UASTARtemp(i);
            RNS(i+2*Mp) = UTSTAR(i) - UTSTARtemp(i);
            RNS(i+3*Mp) = TANBIC(i) - (VAC(i) + UADUCT(i) + UASTAR(i))/(L*RC(i) + VTC(i) + UTSTAR(i));            

        end

        % ----------------------------------------------- Evaluate Jacobian        
        for i = 1:Mp 

            if (Chord_flag == 1) && (strcmp(ChordMethod,'CLmax') == 1)
                % ==========
                % 11/16/2011 BEPPS: New drag treatment...   
                JNS(i     ,i+  Mp) =  3*(VAC(i) + UADUCT(i)) + 4*UASTAR(i)   + 2 * (L*RC(i) +    VTC(i) + 2*UTSTAR(i)) * CD(i)/CLmax(i);                              

                JNS(i     ,i+2*Mp) = - (L*RC(i) +    VTC(i)  + 4*UTSTAR(i))  + 2 *  (VAC(i) + UADUCT(i) + 2*UASTAR(i)) * CD(i)/CLmax(i);               
                % ==========
            else                  
                % 11/16/2011 BEPPS: This was the implementation in EppsOptimizer06.m...
                %                  
                % JNS(i     ,i+  Mp) =    3*(VAC(i) + UADUCT(i)) + 4*UASTAR(i) ...
                %                      + (1/pi)*CD(i)*CoD(i)*dVdG(i,i)*(L*RC(i) + VTC(i) + UTSTAR(i)) ...
                %                      + (1/pi)*CD(i)*CoD(i)*VSTAR(i) * UTHIF(i,i); 
                % 
                % JNS(i     ,i+2*Mp) = - (L*RC(i) +    VTC(i) + 4*UTSTAR(i)) ... 
                %                      + ( VAC(i) + UADUCT(i) + 2*UASTAR(i))*(1/(2*pi))*CD(i)*CoD(i)*dVdG(i,i);                 
                % ------ 
                % Inviscid terms:           
                JNS(i     ,i+  Mp) =  3*(VAC(i) + UADUCT(i)) + 4*UASTAR(i);                              

                JNS(i     ,i+2*Mp) = - (L*RC(i) +    VTC(i)  + 4*UTSTAR(i));               
                % ------
            end

            JNS(i+  Mp,1:Mp  ) = - UAHIF(i,1:Mp);
            JNS(i+  Mp,i+  Mp) = 1;
            JNS(i+2*Mp,1:Mp  ) = - UTHIF(i,1:Mp);
            JNS(i+2*Mp,i+2*Mp) = 1;

            JNS(i+3*Mp,i+  Mp) = -1       /(L*RC(i) + VTC(i) + UTSTAR(i));
            JNS(i+3*Mp,i+2*Mp) = TANBIC(i)/(L*RC(i) + VTC(i) + UTSTAR(i));
            JNS(i+3*Mp,i+3*Mp) = 1;
        end

        % ----------------------------- Update Newton solver vector of unknowns
        DX     = linsolve(JNS,-RNS);

        G      = G      + relax*DX( 1:Mp      ) ;
        UASTAR = UASTAR + relax*DX((1:Mp)+  Mp)';
        UTSTAR = UTSTAR + relax*DX((1:Mp)+2*Mp)';
        TANBIC = TANBIC + relax*DX((1:Mp)+3*Mp)'; 

        % Smooth the inflow angle for numerical stability:
        TANBICsmooth = TANBIC * Bsmooth;
        % -----------------------------------------------------------------

        % -----------------------------------------------------------------
        if any(isnan(DX))  || ~isreal(DX)
             G      = 0*RC';
             UASTAR = 0*RC;
             UTSTAR = 0*RC;
             TANBIC = TANBC;

             disp(' ')
             disp('<WARNING>')
             disp('<WARNING> DX == NaN or imaginary... crash avoided...')
             disp('<WARNING>')
             disp(' ')
             G_iter = 999;
        end
        % ----------------------------------------------------------------- 
        
        % ----------------------------------------------------------------- 
        if Duct_flag == 1            
            % Update the influence of propeller on duct
            for m = 1:Mp;
                DAHIF(:,m) = DAHIF_times_TANBIC(:,m) / TANBICsmooth(m);
                DRHIF(:,m) = DRHIF_times_TANBIC(:,m) / TANBICsmooth(m);
            end
            
            % Update induced velocities at the duct (influence of propeller on duct)
            UARING = (DAHIF*G)';  
            URRING = (DRHIF*G)'; 
        
            % Update duct circulation such that CTD == CTDdes            
            [junk,GdNEW] = Duct_Thrust(XdRING,Rduct_oR,VARING,UARING,URRING,GdRING,Gd,CDd,CTDdes);
            
            Gd = 0.5*GdNEW + (1-0.5)*Gd;
            
            % Update the induced velocities at the propeller (influence of duct on propeller)
            UADUCT =  UADIF*Gd;                        
        end
        % ----------------------------------------------------------------- 
        
        % -----------------------------------------------------------------
        % END: "Robust" method (EppsOptimizer06.m)
        % -!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!- 

    
        % % -----------------------------------------------------------------        
        % % "SIMPLE (INCORRECT) OPTIMIZER"        
        % % -----------------------------------------------------------------
        % % ---------------------- Set up simultaneous equations for G and LM
        % A  = zeros(Mp,Mp);         % A matrix for linear system of equations
        % B  = zeros(Mp,1);          % B matrix for linear system of equations
        % 
        % for i = 1:Mp                           % for each equation for G(i)
        %     for m = 1:Mp                       % for each vortex panel, m        
        %         A(i,m) =  UAHIF(m,i)*RC(m)*DR(m)    ...                
        %                 + UAHIF(i,m)*RC(i)*DR(i); 
        %     end   
        % 
        %     B(i)  = -VAC(i)*RC(i)*DR(i) ...
        %             -(1/(2*pi))*sum(CD.*dVdG(:,i)'.*CoD.*(L*RC + VTC + UTSTAR).*RC.*DR) ...
        %             -(1/(2*pi))*sum(CD.*VSTAR.*CoD.*UTHIF(:,i)'.*RC.*DR);             
        % end
        % 
        % % -------------------------------- Solve linear system of equations
        % GLM = linsolve(A,B);             
        % G   = GLM(1:Mp);                     % G is the first Mp entries
        % 
        % % Update induced velocities (influence of propeller on propeller)
        % UASTAR = (UAHIF*G)';  UTSTAR = (UTHIF*G)';  
        % 
        % % -------------------------  Repair outlier {UASTAR, UTSTAR, URSTAR} values
        % UASTAR = RepairSpline(RC,UASTAR,'UASTAR',Plot_flag*Hvel);
        % UTSTAR = RepairSpline(RC,UTSTAR,'UTSTAR',Plot_flag*Hvel);
        %
        % % -------------------------------------------- Update TANBIC
        % TANBIC = (VAC + UASTAR)./(L*RC + VTC + UTSTAR);
        %          
        % % -----------------------------------------------------------------        
        % % END "SIMPLE (INCORRECT) OPTIMIZER"        
        % % -----------------------------------------------------------------           
    end 
    % ---------------------------------------------------------------------
    % END UPDATE: G, UASTAR, UTSTAR, TANBIC, (duct stuff), LM
    % =@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=@=
    


    % ------------------------------------ Update VSTAR and its derivatives
    VSTAR  = sqrt((VAC+UADUCT+UASTAR).^2 + (L*RC+VTC+UTSTAR).^2);       
    
    % 11/16/2011 BEPPS: Inclusion of these drag terms has negligable effect ( O(10^-4) ) on the converged circulation distribution.
    % 
    % betaIC = atand(TANBIC);                                     % [deg] 
    % 
    % 
    % DUADUT = -(L*RC + VTC + 2*UTSTAR)./(VAC + UADUCT + 2*UASTAR);            % d(UASTAR)/d(UTSTAR)
    % 
    % for i = 1:Mp
    %     dVdG(:,i) = sind(betaIC').*DUADUT(:).*UTHIF(:,i) + cosd(betaIC').*UTHIF(:,i);
    % end            
    % ---------------------------------------------------------------------
    
    % --------------------- Update the vortex Horseshoe Influence Functions
    if Wake_flag == 0 

        % 11/17/2011 BEPPS: If you use TANBICsmooth here, then the circulation and 
        %                   induced velocities in the turbine case will be wavy!
        [UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBIC,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);
        
    else
        % Update radial induced velocity
        URSTAR = (URHIF*G)';

        disp('ERROR: Wake_Geometry and Wake_Horseshoe are no longer supported...please re-implement these...')
            % disp('Beginning to update {UAHIF,UTHIF,URHIF}'), 
            % [WX,WY,WZ]          =  Wake_Geometry(Mp,RC,RV,TANBIV,VAC+UADUCT,UASTAR,URSTAR,CTPdes);
            % %[UAHIF,UTHIF,URHIF] = Wake_Horseshoe(Mp,Z,RC,SCF,WX,WY,WZ,epsilon);
            % [UAHIF,UTHIF,URHIF] = Wake_Horseshoe(Mp,Z,TANBIV,RC,RV,SCF,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR,WX,WY,WZ,epsilon);
            % disp('Done updating {UAHIF,UTHIF,URHIF}'), 
    end        
    % ---------------------------------------------------------------------    
    

    % ------------------------------------------- Update chord distribution
    if (Chord_flag == 1) && (G_iter <= ITER)
        
        if strcmp(ChordMethod,'CLmax') == 1
            % -----------------------------------------------------------------  
            % Method: Choose chord length based on CLmax
            % -----------------------------------------------------------------        
            CoD = 2*pi*G'./(VSTAR.*CLmax); % scale CoD to keep CL == CLmax
            % ----------------------------------------------------------------- 

        elseif strcmp(ChordMethod,'ConeyPLL') == 1
            % -----------------------------------------------------------------
            % Method: (Coney, 1989) cavitation method -- ASSUMES GIVEN THICKNESS DISTRIBUTION t0oD      
            % -----------------------------------------------------------------   
            SIGMA = SIGMAs./VSTAR.^2;    % local cavitation number

            f0oD = (2*pi*G'./ VSTAR) .* f0octilde ./ CLItilde;

            CoD  = (8.09*f0oD+3.033*t0oD)./(2*SIGMA) + sqrt((8.09*f0oD+3.033*t0oD).^2 + 4* SIGMA .* (26.67*f0oD.^2 + 10*f0oD.*t0oD) )./(2*SIGMA);
            % ----------------------------------------------------------------- 
            
            
        elseif strcmp(ChordMethod,'FAST2011dCTP') == 1
            % -----------------------------------------------------------------  
            % Method: see (Epps et al., FAST'2011)
            % -----------------------------------------------------------------
            % Every third iteration, update the chord and thickness.
            if mod(G_iter,3) == 0
            
               [CoD, t0oc, C_res] = Chord_FAST2011_dCTP(SIGMAh,CTDESh, Jh,CoD,t0oc,...      
                                     Propeller_flag,Viscous_flag,Hub_flag,Duct_flag, ...
                                     Z,Mp,ITER,Rhv,RC,RV,DR,Rhub_oR,VMIV,CTD,...
                                     L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,CD,...
                                     R,rho,Vs,N);            
            end
            
        elseif strcmp(ChordMethod,'FAST2011dVAC') == 1
            % -----------------------------------------------------------------  
            % Method: see (Epps et al., FAST'2011)
            % -----------------------------------------------------------------
            % Every third iteration, update the chord and thickness.
            if mod(G_iter,3) == 0
                
               [CoD, t0oc, C_res] = Chord_FAST2011_dVAC(dVAC,SIGMAs,L,RC,VAC,VTC,UADUCT,UASTAR,UTSTAR,G,CoD,t0oc,...
                                           DR,CD,Z,Js,VMIV,Hub_flag,Rhub_oR,Rhv,CTD,Mp,R,rho,Vs,N);
            end
            
        elseif strcmp(ChordMethod,'Brizzolara2007') == 1
            % -----------------------------------------------------------------  
            % Method: see (Brizzolara et al., 2007)
            % -----------------------------------------------------------------
            % Every third iteration, update the chord and thickness.
            if mod(G_iter,3) == 0
                
%                 save temp
%                 return
                
               [CT,CQ,CP,KT,KQ, CTH,TAU, Ja,Jw,VMWV, EFFYo, EFFY,ADEFFY,QF] = ...
                    Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,UADUCT,CD,CoD,G,Z,Js,VMIV,Hub_flag,Rhub_oR,Rhv,CTD);

               [CoD, t0oc, C_res] = Chord_Brizzolara(rho,n,D,Vs,H,Mp,Z,KT,KQ, RC,G,VSTAR,TANBIC,CoD,t0oc);
            end
        end
        
        % -----------------------------------------------------------------
        % Scale CoD to give specified Expanded Area Ratio (EAR == EARspec)
        if EARspec > 0
           EAR = (2*Z/pi) * trapz(   linspace(Rhub_oR,1,100), interp1(RC,CoD, linspace(Rhub_oR,1,100), 'spline','extrap')   );  
 
           CoD = (EARspec/EAR) * CoD; 
        end
        % -----------------------------------------------------------------
        
        % -----------------------------------------------------------------
        if all(CoD == 0) || any(isnan(CoD))
             G      = 0*RC';
             UASTAR = 0*RC;
             UTSTAR = 0*RC;
             TANBIC = TANBC;
             CoD    = 0*RC;
             
             disp(' ')
             disp('<WARNING>')
             disp('<WARNING> CoD == NaN or zero... crash avoided...')
             disp('<WARNING>')
             disp(' ')
             G_iter = 999;
        end
        % -----------------------------------------------------------------            
    end
    % ---------------------------------------------------------------------
        

    
    % ---------------------------------------------------------------------
    if Plot_flag == 1
        % Cycle through the colors    
        color_index = mod(color_count-1,size(CLR,1))+1;

        color_count = color_count + 1;
                
        % ---------------------------------------------------------------------    
        figure( Hgamma),
        delete(HHgamma),
               HHgamma = plot(RC,G,'.-','Linewidth',2,'Color',CLR(color_index,:));
        % ---------------------------------------------------------------------     

        % --------------------------------------------------------------------- 
        figure( Hvel),
        delete(HHvel),
        if Duct_flag == 1
               HHvel = plot(RC,UASTAR,'b.-',RC,UTSTAR,'r.-',RC,UADUCT,'g.-');
        else
               HHvel = plot(RC,UASTAR,'b.-',RC,UTSTAR,'r.-');
        end
        % --------------------------------------------------------------------- 
        
        % --------------------------------------------------------------------- 
        figure( Hbeta),
        delete(HHbeta),
               HHbeta = plot(RC,TANBIC,'.-','Linewidth',2,'Color',CLR(color_index,:));
        % ---------------------------------------------------------------------         
        
        % ---------------------------------------------------------------------
        if Chord_flag == 1  % if chord optimization WAS NOT performed, then no need to replot the blade
                            % else
                            % chord optimization WAS performed, so we need to be careful to iterpolate in such a way as to  
                            % preserve a   zero-chord-length tip (as in the case of a free-tip propeller) or 
                            % preserve a finite-chord-length tip (as in the case of a ducted   propeller).
            figure( Hcod),
            delete(HHcod),
                % ------------------------------------------------------------------------- Interpolate chord
                    if (Duct_flag == 0) || (Rduct_oR > 1.001)   

                        XXCoD  = InterpolateChord(RC,CoD ,XXR);  % yields zero   chord at the tip       
                    else                                                
                        XXCoD  =            pchip(RC,CoD ,XXR);  % yields finite chord at the tip   
                    end     
                % -------------------------------------------------------------------------          
            HHcod = plot(XXR,-XXCoD,'k-',XXR,XXCoD,'k-',RC,-CoD,'b.',RC,CoD,'b.');
        end
        % ---------------------------------------------------------------------

        pause(0.0001),   drawnow,
    end  % if Plot_flag == 1     
    % ---------------------------------------------------------------------
    

    % ---------------------------------------------------------------------
    % ---------------------------------- Prepare for the next iteration
    G_res   = abs((G - G_last)./G);    % residual G
    G_last  = G;                       % the last value of G
    Gd_res  = abs((Gd - Gd_last)/Gd);  % residual Gd
    Gd_last = Gd;                      % the last value of Gd
    LM_last = LM;                      % last value of the Lagrange Multiplier    
    
    if G_iter < 10
        disp(['The max  G_res for iteration  ',num2str(G_iter),' is: ',num2str(max(G_res))]),  
    else
        disp(['The max  G_res for iteration ',num2str(G_iter),' is: ',num2str(max(G_res))]),  
    end    
    
%     if Duct_flag == 1
%         disp(['The     Gd_res for iteration ',num2str(G_iter),' is: ',num2str(Gd_res)]),        
%     end
    % ---------------------------------------------------------------------
    
    % ---------------------------------------------------------------------
    if (Chord_flag == 1)  &&  ( (strcmp(ChordMethod,'FAST2011dCTP') == 1) || (strcmp(ChordMethod,'FAST2011dVAC') == 1) || (strcmp(ChordMethod,'Brizzolara2007') == 1) )
        if mod(G_iter,3) ~= 0
            G_res = 1;
        else
            G_res = max([G_res;C_res/10]);
        end
    end
    % ---------------------------------------------------------------------
    
    % ---------------------------------------------------------------------
    G_iter  = G_iter + 1;              % iteration in the G loop
    % ---------------------------------------------------------------------   
end % while G_iter <= ITER && any([G_res;Gd_res] > 1e-4)  % (WHILE LOOP G1)
% -------------------------------------------------------------------------

if G_iter > ITER
    disp(' '),
    disp('<WARNING> While loop G1 did NOT converge.'),
    Converge_flag = 0;
else
    disp(' '),
    disp('Design optimization complete'),
    Converge_flag = 1;
end


% =========================================================================
% =========================================================================
%%

% =========================================================================
% 11/17/2011 BEPPS: Note that the OpenProp 3.2.0 version of Unload_Blade.m 
%                   does not do a very good job of scaling the circulation
%                   such that CT == CTdes.  This code should be improved.
%
% If required, unload the hub and tip, then rescale the circulation
% distribution to get the desired value of the thrust coefficient.
if Hub_flag && (HUF > 0 || TUF > 0)                       % (IF STATEMENT U1)

    if Duct_flag == 0
        DAHIF_times_TANBIC  = 0;
        DRHIF_times_TANBIC 	= 0;
        XdRING              = 0;
        Rduct_oR            = 1;
        VARING              = 0;
        GdRING              = 0;
        Gd                  = 0;
        CDd                 = 0;
        CTDdes              = 0;
    end
    
    [G,UASTAR,UTSTAR,TANBIC,UARING,URRING,Gd,UADUCT] = ...
                                    Unload_Blade(HUF,TUF,RC,Rhub_oR, G,  VAC,VTC, TANBIC,RV,DR,L,Mp,Z, ...
                                                 Hub_flag,ITER,Bsmooth,CD,CoD,Js,VMIV,Rhv,CTPdes,...
                                                 Duct_flag,UADUCT,XdRING,Rduct_oR,VARING,GdRING,UADIF,Gd,CDd,CTDdes,...
                                                 DAHIF_times_TANBIC,DRHIF_times_TANBIC,...
                                                 Plot_flag,Hgamma,HHgamma,Hvel,HHvel,Hbeta,HHbeta);    
    
end
% =========================================================================

% if Rroot > Rhub 
%     % --------------------------------------------------------------------- 
%     % Develop code for turbine case with circular blade sections near the root.
%     % ---------------------------------------------------------------------
% end


% --------------------- Compute thrust & torque coefficients and efficiency
% Compute the actual CT for the duct with the current circulation Gd
if Duct_flag == 1
    [CTD,junk] = Duct_Thrust(XdRING,Rduct_oR,VARING,UARING,URRING,GdRING,Gd,CDd,CTDdes);
end
 
[CT,CQ,CP,KT,KQ, CTH,TAU, Ja,Jw,VMWV, EFFYo, EFFY,ADEFFY,QF,QFo,QFw] = ...
	Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,UADUCT,CD,CoD,G,Z,Js,VMIV,Hub_flag,Rhub_oR,Rhv,CTD);


% -------------------------------------------------------------------------
CL   = 2*pi*G'./(VSTAR.*CoD);     % lift coefficient

% Expanded Area Ratio
EAR = (2*Z/pi) * trapz(   linspace(Rhub_oR,1,100), interp1(RC,CoD, linspace(Rhub_oR,1,100), 'spline','extrap')   );  

% -------------------------------------------------------------------------


% ------------------------------------------- Update chord distribution
if (Chord_flag == 1) && ( (strcmp(ChordMethod,'FAST2011dCTP') == 1) || (strcmp(ChordMethod,'FAST2011dVAC') == 1) || (strcmp(ChordMethod,'Brizzolara2007') == 1) )
    % Use with chord optimization methods that optimize t0oc
    t0oD = t0oc.*CoD;
else   
    % Use when no chord optimization or with chord optimization methods 
    % that do not optimize t0oc (e.g. CLmax or Coney use given t0oD)
    t0oc = t0oD ./ CoD;    
end
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
    disp(' '),
    disp('Forces after circulation optimization:')
    
if Propeller_flag == 1
   
    disp(['    Js = ',num2str(Js)]),    
    disp(['    KT = ',num2str(KT)]),   
    disp(['    KQ = ',num2str(KQ)]),   
    disp(['    CT = ',num2str(CT)]),
    if abs(VMIV - 1) > 1e-8,  % i.e. if VMIV is not equal to 1 
    disp(['    Ja = ',num2str(Ja)]),
    end
    disp(['  EFFY = ',num2str(EFFY)]),
    disp(['ADEFFY = ',num2str(ADEFFY)]),    
    disp(['    QF = ',num2str(QF)]),    
    
else
    disp(['L      =  ',num2str(L)]),    
    disp(['CP     = ' ,num2str(CP)]),  
    disp(['CPBetz = ' ,num2str(CPBetz)]),  
    disp(['QF     =  ',num2str(CP/CPBetz)]),      
end
    disp(' ')     
    disp(' ')     
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
% --------------------------- Package results in "design._____" data struct
design.part1      = '------ Section properties, size (1,Mp) ------';
design.RC         = RC;                 % [1 x Mp] control point radii
design.DR         = DR;                 % [1 x Mp] difference in vortex point radii
design.G          = G';                 % [1 x Mp] circulation distribution
design.VAC        = VAC;                % [1 x Mp] 
design.VTC        = VTC;                % [1 x Mp] 
design.UASTAR     = UASTAR;             % [1 x Mp] 
design.UTSTAR     = UTSTAR;             % [1 x Mp] 
if Wake_flag == 1
design.URSTAR     = URSTAR;             % [1 x Mp] 
end
design.VSTAR      = VSTAR;              % [1 x Mp]  
design.TANBC      = TANBC;              % [1 x Mp] 
design.TANBIC     = TANBIC;             % [1 x Mp] 
design.CL         = CL;                 % [1 x Mp] 
design.CD         = CD;                 % [1 x Mp] 
design.CoD        = CoD;                % [1 x Mp]
design.t0oc       = t0oc;               % [1 x Mp] 
design.t0oD       = t0oD;               % [1 x Mp] 

design.part2      = '------ Other properties  ------';
design.converged  = Converge_flag;
design.iteration  = G_iter;
design.RV         = RV;                 % [1 x Mp+1] vortex point radii
design.Rhub_oR    = Rhub_oR;            % [1 x 1]
if Rroot_oR > Rhub_oR, design.Rroot_oR   = Rroot_oR; end
if Rcirc_oR > Rhub_oR, design.Rcirc_oR   = Rcirc_oR; end
design.EAR        = EAR;
design.LM         = LM;                 % [1 x 1]
design.VMIV       = VMIV;               % [1 x 1]
design.VMWV       = VMWV;               % [1 x 1]
design.SIGMAs     = SIGMAs;


if Duct_flag == 1  
    design.part3      = '------ Duct parameters ------';
    design.Rduct_oR   = Rduct_oR;           % [1 x 1]
    design.Cduct_oR   = Cduct_oR;           % [1 x 1]
    design.Xduct_oR   = Xduct_oR;           % [1 x 1]
    design.Gd         = Gd;                 % [1 x 1], duct circulation   
    design.VARING     = VARING;             % [1 x  1]
    
    design.XdRING     = XdRING;             % [1 x Nd], Nd=12 duct vortex rings 
    design.UARING     = UARING;             % [1 x Nd], Nd=12 duct vortex rings
    design.URRING     = URRING;             % [1 x Nd], Nd=12 duct vortex rings
    design.GdRING     = GdRING;             % [1 x Nd], Nd=12 duct vortex rings

    design.DAHIFtT    = DAHIF_times_TANBIC; % [Nd,Mp], Duct Horseshoe Influence Functions (influence of propeller on duct)
    design.DRHIFtT    = DRHIF_times_TANBIC; % [Nd,Mp], Duct Horseshoe Influence Functions (influence of propeller on duct)
    
    

    design.UADIF      = UADIF;              % [1 x Mp] 
    design.UADUCT     = UADUCT;             % [1 x Mp] 
    design.CTPdes     = CTPdes;             % [1 x 1], desired propeller CT
    design.CTDdes     = CTDdes;             % [1 x 1], desired duct      CT

    design.TAU        = TAU;                  % [1 x 1]
    design.CTD        = CTD;                  % [1 x 1]
    design.part4      = '------ Performance metrics ------';
else
    design.part3      = '------ Performance metrics ------';
end


if Propeller_flag == 1
    design.L      = L;
    design.Js     = Js;
    design.KT     = KT;                   % [1 x 1]
    design.KQ     = KQ;                   % [1 x 1]
    design.CT     = CT;                   % [1 x 1]
    design.CQ     = CQ;                   % [1 x 1]
    design.CP     = CP;                   % [1 x 1]
    design.CTH    = CTH;                  % [1 x 1]
    
    if abs(VMIV - 1) > 1e-8,  % i.e. if VMIV is not equal to 1 
    design.EFFYo  = EFFYo;                 % [1 x 1]
    design.Ja     = Ja;
    end
    design.EFFY   = EFFY;                 % [1 x 1],  EFFY = EFFYa by convention (Kerwin)
    design.ADEFFY = ADEFFY;               % [1 x 1], prop.   actuator disk
    design.QF     = QF;                   % [1 x 1]
    
    if (VMIV < 0.05)   % assume design for bollard pull
        design.QFo    = QFo;                  % [1 x 1]
        design.QFw    = QFw;                  % [1 x 1]
    end
else    
    design.L      = L;
    design.Js     = Js;
%     design.KT     = KT;                   % [1 x 1]
%     design.KQ     = KQ;                   % [1 x 1] 
    design.CT     = CT;                   % [1 x 1]
    design.CQ     = CQ;                   % [1 x 1]
    design.CP     = CP;                   % [1 x 1]
    design.CPBetz = CPBetz;               % [1 x 1], turbine actuator disk
    design.QF     = CP/CPBetz;
end
% -------------------------------------------------------------------------
  
% ============================================== END EppsOptimizer Function
% =========================================================================