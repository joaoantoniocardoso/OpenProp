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
% ================================================= LerbsOptimizer Function
% Last Modified: 11/2/2011, Brenden Epps
%
% The LerbsOptimizer function determines the circulation distribution that
% satisfies the input operating conditions and the Lerbs criterion. Returns 
% performance specs, such as thrust coefficient and efficiency, as well as
% the circulation distribution, ect.
%
% -------------------------------------------------------------------------
% Note: This implementation does not include the duct.
% -------------------------------------------------------------------------

function design = LerbsOptimizer(input)

% -------------------------------------------------------------------------
if isfield(input,'Propeller_flag')
    if input.Propeller_flag == 0
        disp('Sorry, LerbsOptimizer is only set up for propellers...please update the code...')
        return
    end
end

if isfield(input,'Duct_flag')
       Duct_flag = input.Duct_flag;
    if Duct_flag == 1
        disp('Sorry, LerbsOptimizer does not model the duct...please update the code...')
        return
    end
else 
    Duct_flag = 0;
end
% -------------------------------------------------------------------------

%% ========================================================================            
% -------------------------------------------------- Unpack input variables
% '------ Performance inputs ------'
Z       = input.Z;           % [1 x 1] number of blades
Js      = input.Js;          % [1 x 1]
L       = pi/Js;           % [1 x 1]
   
% -------------------------------------------------------------------------
if isfield(input,  'Mp'), Mp   = input.Mp;    else  Mp   = 20;  end
if isfield(input,  'Vs'), Vs   = input.Vs;    else  Vs   = 1;   end % m/s
if isfield(input, 'rho'), rho  = input.rho;   else  rho  = 1000;end % kg/m^3

if     isfield(input,'D'),   R = input.D/2;
elseif isfield(input,'R'),   R = input.R;
else                         R = 1;
end

if     isfield(input,'Dhub'),   Rhub = input.Dhub/2;
elseif isfield(input,'Rhub'),   Rhub = input.Rhub;
else                            Rhub = 0.2*R;
end

    % CTdes == desired total thrust coefficient   
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
    
    

% -------------------------------------------------------------------------
% If propeller geometry or inflow is not given, assume propeller 4118 form 
% with uniform axial inflow, no swirl inflow, and section CD == 0.008.
% Note, the real propeller 4118 has CoD == 0.001 at the tip.
if isfield(input,'XR'), 
    XR  = input.XR;
    X1  =  ones(size(XR));
    X0  = zeros(size(XR));
    
    if isfield(input,'XCoD'), 
            XCoD  = input.XCoD;  
                                  
    else
        XXR    = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];    % r/R
        XXCoD  = [0.1600 0.1818 0.2024 0.2196 0.2305 0.2311 0.2173 0.1806 0.1387 0.0100];
             
        XCoD   = interp1(XXR,XXCoD ,XR,'pchip','extrap');

    end 
   
    if isfield(input,'XVA'  ),   XVA = input.XVA;   else   XVA = X1; end
    if isfield(input,'XVT'  ),   XVT = input.XVT;   else   XVT = X0; end
else 
    XR    = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];    % r/R
    X1    =  ones(size(XR));
    X0    = zeros(size(XR));
    XCoD  = [0.1600 0.1818 0.2024 0.2196 0.2305 0.2311 0.2173 0.1806 0.1387 0.0100];          % c/D
    XVA   =       X1;                                      % Va/Vs
    XVT   =       X0;                                      % Vt/Vs
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Blade section properties
if isfield(input,'XCD'), XCD = input.XCD; 
    if length(XCD) == 1, XCD =   XCD*X1;  end
else                     XCD = 0.008*X1;    
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Inflow velocity profiles
if isfield(input,'ri')
    RI  = input.ri/R;

    if isfield(input,'VAI'),   VAI = input.VAI;   else   VAI =  ones(size(ri)); end  % Va/Vs
    if isfield(input,'VTI'),   VTI = input.VTI;   else   VTI = zeros(size(ri)); end  % Vt/Vs       
else
     RI = XR; %  r/R
    VAI = XVA;
    VTI = XVT;  
end

% Smooth the inflow velocity profiles:
VAI = RepairSpline(RI,VAI);
VTI = RepairSpline(RI,VTI);
% -------------------------------------------------------------------------

  
% '------ Computational inputs ------'
Viscous_flag = input.Viscous_flag;  % [1 x 1]   0 == viscous forces off (CD = 0), 1 == viscous forces on
Hub_flag     = input.Hub_flag;      % [1 x 1]   0 == no hub, 1 == hub


% 0 == do not optimize chord lengths, 1 == optimize chord lengths
if isfield(input,'Chord_flag'), Chord_flag = input.Chord_flag;  
                          else  Chord_flag = 0;  end

% -------------------------------------------------------------- Chord_flag
if Chord_flag == 1
        
    % Specified Expanded Area Ratio (EAR)
    if  isfield(input,'EAR'), EARspec = input.EAR; else EARspec = 0; end
    
    
    % Specified maximum allowable CL
    if  isfield(input,'XCLmax')    
        XCLmax = input.XCLmax; if length(XCLmax) == 1, XCLmax = XCLmax*X1;  end
    else
        XCLmax = 0.5 + (0.2-0.5)/(1-XR(1)) * (XR-XR(1));  % XCLmax == 0.5 at root and 0.2 at tip
    end
else
    XCLmax = X1;
end
% -------------------------------------------------------------------------

% ------------------------------------------------------------ Viscous_flag 
if Viscous_flag == 0
    XCD   = X0;
end
% -------------------------------------------------------------------------


% ---------------------------------------------------- Computational inputs
if isfield(input,'ITER'), ITER = input.ITER;  else  ITER = 20;   end
if isfield(input,'HUF' ), HUF  = input.HUF;   else  HUF  = 0;    end
if isfield(input,'TUF' ), TUF  = input.TUF;   else  TUF  = 0;    end
if isfield(input,'Rhv' ), Rhv  = input.Rhv;   else  Rhv  = 0.5;  end

%------------------- Initialize duct related variables needed in functions
% Note: This implementation does not include the duct.
CTD      = 0;     % CT for duct
Rduct_oR = 1;
UADUCT   = 0;
% =========================================================================


% =========================================================================

% ------ Compute the Volumetric Mean Inflow Velocity, Kerwin eqn 163, p.138
Rhub_oR  = Rhub/R;               % [ ], hub       radius / propeller radius
XRtemp   = linspace(Rhub_oR,1,100);                           % (temp) radius / propeller radius
VAItemp  = interp1(RI,VAI,XRtemp,'pchip','extrap');            % (temp) axial inflow velocity /ship velocity
VMIV     = 2*trapz(XRtemp,XRtemp.*VAItemp)/(1-Rhub_oR^2);  % [ ], VMIV/ship velocity

% -- Compute cosine spaced vortex & control pt. radii, Kerwin eqn 255, p179
% RV = radius of vortex  points / propeller radius [ ]
% RC = radius of control points / propeller radius [ ]
DEL  = pi/(2*Mp);
Rdif = 0.5*(1 - Rhub_oR);

for m = 1:Mp+1
    RV(m) = Rhub_oR + Rdif*(1-cos(2*(m-1)*DEL));
end

for n = 1:Mp
    RC(n) = Rhub_oR + Rdif*(1-cos((2*n-1)*DEL));
end

DR = diff(RV);

% ------------ Interpolate Va, Vt, Cd, and c/D at vortices & control points
VAC  = pchip(RI,VAI   ,RC);   % axial      inflow vel. / ship vel. at ctrl pts
VTC  = pchip(RI,VTI   ,RC);   % tangential inflow vel. / ship vel. at ctrl pts
CD   = pchip(XR,XCD   ,RC);   % section drag coefficient           at ctrl pts
CLmax= pchip(XR,XCLmax,RC);   % maximum allowable lift coefficient at ctrl pts

if (abs(XR(end)-1) < 1e-4) && (XCoD(end) <= 0.01)  % if XR == 1 and XCoD == 0
    
    CoD  = InterpolateChord(XR,XCoD,RC);   % section chord / propeller diameter at ctrl pts
else
    CoD  =            pchip(XR,XCoD,RC);   % section chord / propeller diameter at ctrl pts
end


% --- Do first estimation of tanBi based on 90% of actuator disk efficiency
EDISK  = 1.8/(1+sqrt(1+CTdes/VMIV^2));          % efficiency estimate

TANBC  = VAC./(pi.*RC./Js + VTC);               % tan(Beta) at control pts.

VBAC   = VTC.*TANBC./VAC;

TANBXC = TANBC.*sqrt(VMIV./(VAC-VBAC))/EDISK;   % estimation of tan(BetaI)

% -------------------------- Unload hub and tip as specified by HUF and TUF
Rmean = (Rhub_oR+1)/2;

for i=1:Mp
    if RC(i)<Rmean
        HRF=HUF;
    else
        HRF=TUF;
    end
    DTANB     = HRF*(TANBXC(i)-TANBC(i))*((RC(i)-Rmean)/(Rhub_oR-Rmean))^2;
    TANBXC(i) = TANBXC(i)-DTANB;
end
% --------------------------

% ------- Iterate to scale tanBi to get desired value of thrust coefficient
CT_iter   = 1;                          % iteration in the CT loop
CT_res    = 1;                          % residual CT
CT_last2  = 0;                          % the CT prior to the last CT
CT_last1  = 0;                          % the last value of CT
TMF_last2 = 0;                          % the TMF prior to the last TMF
TMF_last1 = 0;                          % the last value of TMF
B         = zeros(Mp,1);
A         = zeros(Mp,Mp);

while CT_iter <= ITER & any(CT_res > 1e-5)               % (WHILE LOOP MF1)
    % disp(['LerbsOptimizer iteration...',num2str(CT_iter)])
    
    % ------------------------------ Scale TANBI by a multiplication factor
    if CT_iter == 1
        TMF = 1;

    elseif CT_iter == 2
        TMF = 1+(CTdes-CT)/(5*CTdes);

    elseif CT_iter > 2
        TMF = TMF_last1 + (TMF_last1-TMF_last2)*(CTdes-CT_last1)/...
                          ( CT_last1- CT_last2);
    end

    TANBIC = TMF.*TANBXC;              % TMF = TANBIC Multiplication Factor

    % ------ Compute the vortex Horseshoe Influence Functions, Kerwin p.179
    [UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBIC,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);
    
    % --------- Solve simutaneous equations for circulation strengths G(i)
    %           A = (2*pi*R*Vs)*(LHS of eqn 258), B = RHS of eqn 258, p.181
    %           NOTE: There is an error in the LHS of Kerwin, eqn 258.  It 
    %                 should have a leading 1/Vs to agree with eqn 254.

    for n = 1:Mp                            % for each control point, n
        B(n) = VAC(n)*((TANBIC(n)/TANBC(n))-1);

        for m = 1:Mp                            % for each vortex  panel, m
            A(n,m) = UAHIF(n,m)-UTHIF(n,m)*TANBIC(n);
        end
    end

    G = linsolve(A,B);                     % G = Gamma / (2*pi*R*Vs) = [ ]

    
    % -- Compute induced velocities at control points. Kerwin eqn 254, p179
    UASTAR = (UAHIF*G)';  
    UTSTAR = (UTHIF*G)';
    
    VSTAR  = sqrt((VAC+UADUCT+UASTAR).^2 + (L*RC+VTC+UTSTAR).^2); 
    
    % ------------------------------------------- Update chord distribution
    if Chord_flag == 1
        % -----------------------------------------------------------------  
        % (1) Choose chord length based on CLmax, or
        % -----------------------------------------------------------------        
        CoD = 2*pi*abs(G')./(VSTAR.*CLmax); % scale CoD to keep CL == CLmax
        % -----------------------------------------------------------------   
        
        % -----------------------------------------------------------------  
        % (2) IMPLEMENT FAST'2011 code here, or Coney method, or...
        % ----------------------------------------------------------------- 
        
        
        % -----------------------------------------------------------------
        % Scale CoD to give specified Expanded Area Ratio (EAR == EARspec)
        if EARspec > 0
           EAR = (2*Z/pi) * trapz(   linspace(Rhub_oR,1,100), interp1(RC,CoD, linspace(Rhub_oR,1,100), 'spline','extrap')   );  
 
           CoD = (EARspec/EAR) * CoD; 
        end
        % -----------------------------------------------------------------
    end
    % ---------------------------------------------------------------------     
     
     
    % ----------------- Compute thrust & torque coefficients and efficiency
    [CT,CQ,CP,KT,KQ, CTH,TAU, Ja,Jw,VMWV, EFFYo, EFFY,ADEFFY,QF] = ...
          Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,UADUCT,CD,CoD,G,Z,Js,VMIV,Hub_flag,Rhub_oR,Rhv,CTD);                               

    % -------------------------------------- Prepare for the next iteration
    CT_iter   = CT_iter + 1;           % iteration in the CT loop
    CT_res    = abs(CT - CTdes);       % residual CT
    CT_last2  = CT_last1;              % the CT prior to the last CT
    CT_last1  = CT;                    % the last value of CT
    TMF_last2 = TMF_last1;             % the TMF prior to the last TMF
    TMF_last1 = TMF;                   % the last value of TMF

    
%      disp(['LerbsOptimizer iteration: ',num2str(CT_iter-1),',  max percent error: ',num2str(max(CT_res))]),
%     disp(' '),    
    
    if CT_iter > ITER
        disp('WARNING: Lerbs optimization did NOT converge.'),
    end
end                                                  % (END WHILE LOOP MF1)
% ------------------------------------------------------------------------- 



VSTAR  = sqrt((VAC+UASTAR).^2 + (pi*RC/Js+VTC+UTSTAR).^2);        % V* / Vs
CL     = 2*pi*G'./(VSTAR.*CoD);  % [ ] lift coefficient

% Expanded Area Ratio
EAR = (2*Z/pi) * trapz(   linspace(Rhub_oR,1,100), interp1(RC,CoD, linspace(Rhub_oR,1,100), 'spline','extrap')   );  


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
    disp(' '),
    disp('Forces after circulation optimization:')

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
    disp(' ')     
    disp(' ')      

% -------------------------------------------------------------------------
% --------------------------- Package results in "design._____" data struct
design.part1      = '------ Section properties, size (1,Mp) ------';
design.RC         = RC;                 % [1 x Mp] control point radii
design.G          = G';                 % [1 x Mp] circulation distribution
design.VAC        = VAC;                % [1 x Mp] 
design.VTC        = VTC;                % [1 x Mp] 
design.UASTAR     = UASTAR;             % [1 x Mp] 
design.UTSTAR     = UTSTAR;             % [1 x Mp] 
design.VSTAR      = VSTAR;              % [1 x Mp]  
design.TANBC      = TANBC;              % [1 x Mp] 
design.TANBIC     = TANBIC;             % [1 x Mp] 
design.CL         = CL;                 % [1 x Mp] 
design.CD         = CD;                 % [1 x Mp] 
design.CoD        = CoD;                % [1 x Mp]

design.part2      = '------ Other properties  ------';
design.RV         = RV;
design.Rhub_oR    = Rhub_oR;            % [1 x 1]
design.EAR        = EAR;
design.VMIV       = VMIV;               % [1 x 1]
design.VMWV       = VMWV;               % [1 x 1]

design.part3      = '------ Performance metrics ------';
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
design.EFFY   = EFFY;                 % [1 x 1]
design.ADEFFY = ADEFFY;               % [1 x 1], prop.   actuator disk
design.QF     = QF;                   % [1 x 1]

% ============================================= END LerbsOptimizer Function
% =========================================================================