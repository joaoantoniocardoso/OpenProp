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
% ================================================ LerbsParametric Function
%
% The LerbsOptimizer function determines the  "optimum" circulation distribution;
% it satisfies the parinput operating conditions and the Lerbs criterion.  The
% LerbsOptimizer function returns performance specs, such as thrust coefficient and
% efficiency, as well as the circulation distribution, ect.
%
% -------------------------------------------------------------------------

function paroutput = LerbsParametric(parinput)
           
% -------------------------------------------------------------------------
if isfield(parinput,'Propeller_flag')
    if parinput.Propeller_flag == 0
        disp('Sorry, LerbsParametric is only set up for propellers...please update the code...')
        return
    end
end

if isfield(parinput,'Duct_flag')
       Duct_flag = parinput.Duct_flag;
    if Duct_flag == 1
        disp('Sorry, LerbsParametric does not model the duct...please update the code...')
        return
    end
else
    Duct_flag = 0;
end
% -------------------------------------------------------------------------

%% ------------------------------------------------- Unpack input variables
% '------ Performance inputs ------'
Zall    = parinput.Z;           % [numberZ x 1], [ ]   number of blades
Nall    = parinput.N;           % [numberN x 1], [RPM] rotation rate
Dall    = parinput.D;           % [numberD x 1], [m]   propeller diameter

Vs      = parinput.Vs;          % [1 x 1]
THRUST  = parinput.THRUST;       % [1 x 1]
% '------ Geometry inputs ------'
Mp      = parinput.Mp;          % [1 x 1]
Rhub    = parinput.Rhub;        % [1 x 1]


% -------------------------------------------------------------------------
% If propeller geometry or inflow is not given, assume propeller 4118 form 
% with uniform axial inflow, no swirl inflow, and section CD == 0.008.
% Note, the real propeller 4118 has CoD == 0.001 at the tip.
if isfield(parinput,'XR'), 
    XR  =  parinput.XR;
    X1  =  ones(size(XR));
    X0  = zeros(size(XR));
    
    if isfield(parinput,'XCoD'), 
            XCoD  = parinput.XCoD;  
                                  
    else
        XXR    = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];    % r/R
        XXCoD  = [0.1600 0.1818 0.2024 0.2196 0.2305 0.2311 0.2173 0.1806 0.1387 0.0100];
             
        XCoD   = interp1(XXR,XXCoD ,XR,'pchip','extrap');

    end 
   
    if isfield(parinput,'XVA'  ),   XVA = parinput.XVA;   else   XVA = X1; end
    if isfield(parinput,'XVT'  ),   XVT = parinput.XVT;   else   XVT = X0; end
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
if isfield(parinput,'XCD'), XCD = parinput.XCD; 
    if length(XCD) == 1,    XCD =   XCD*X1;  end
else                        XCD = 0.008*X1;    
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Inflow velocity profiles
if isfield(parinput,'ri')    
    ri  = parinput.ri;

    if isfield(parinput,'VAI'),   VAI = parinput.VAI;   else   VAI =  ones(size(ri)); end  % Va/Vs
    if isfield(parinput,'VTI'),   VTI = parinput.VTI;   else   VTI = zeros(size(ri)); end  % Vt/Vs       
else    
    ri  = XR*max(Dall/2);     % set bigger than any possible radii needed 
    VAI =  ones(size(XR));
    VTI = zeros(size(XR)); 
end 

% Smooth the inflow velocity profiles:
VAI = RepairSpline(ri,VAI);
VTI = RepairSpline(ri,VTI);
% -------------------------------------------------------------------------



% '------ Computational inputs ------'
Viscous_flag = parinput.Viscous_flag; % 0 == viscous forces off (CD = 0), 1 == viscous forces on
Hub_flag     = parinput.Hub_flag;     % 0 == no hub, 1 == hub

% 0 == do not optimize chord lengths, 1 == optimize chord lengths
if isfield(parinput,'Chord_flag'), Chord_flag = parinput.Chord_flag;  
                             else  Chord_flag = 0;  end

% -------------------------------------------------------------- Chord_flag
if Chord_flag == 1
        
    % Specified Expanded Area Ratio (EAR)
    if  isfield(parinput,'EAR'), EARspec = parinput.EAR; else EARspec = 0; end
    
    
    % Specified maximum allowable CL
    if   isfield(parinput,'XCLmax')    
        XCLmax = parinput.XCLmax; if length(XCLmax) == 1, XCLmax = XCLmax*X1;  end
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
if isfield(parinput,'ITER'), ITER = parinput.ITER;  else  ITER = 20;   end
if isfield(parinput,'HUF' ), HUF  = parinput.HUF;   else  HUF  = 0;    end
if isfield(parinput,'TUF' ), TUF  = parinput.TUF;   else  TUF  = 0;    end
if isfield(parinput,'Rhv' ), Rhv  = parinput.Rhv;   else  Rhv  = 0.5;  end

if isfield(parinput,'rho' ), rho  = parinput.rho;   else  Rhv  = 1000;  end


%------------------- Initialize duct related variables needed in functions
% Note: This implementation does not include the duct.
CTD      = 0;     % CT for duct
Rduct_oR = 1;
UADUCT   = 0;
% -------------------------------------------------------------------------
  
% =========================================================================

% -------------------------------------------------------------------------
%         From first to last                                              %
%         The peak is never passed                                        %
%         Something always fires the light that gets in your eyes         %
%         One moment's high, and glory rolls on by                        %
%         Like a streak of lightning                                      %
%         That flashes and fades in the summer sky                        %
% -------------------------------------------------------------------------


% =========================================================================
% ----------------------------------- Initialize output performance metrics
numberZ = length(Zall);
numberN = length(Nall);
numberD = length(Dall);

JS     = zeros(numberZ,numberN,numberD);
KT     = zeros(numberZ,numberN,numberD);
KQ     = zeros(numberZ,numberN,numberD);
CT     = zeros(numberZ,numberN,numberD);
CQ     = zeros(numberZ,numberN,numberD);
CP     = zeros(numberZ,numberN,numberD);
CTH    = zeros(numberZ,numberN,numberD);
VMIV   = zeros(numberZ,numberN,numberD);
VMWV   = zeros(numberZ,numberN,numberD);
  EFFYo= zeros(numberZ,numberN,numberD);
  EFFY = zeros(numberZ,numberN,numberD);
ADEFFY = zeros(numberZ,numberN,numberD);
QF     = zeros(numberZ,numberN,numberD);
EAR    = zeros(numberZ,numberN,numberD);


propNo = 1;                       % propeller number
Nprops = numberZ*numberN*numberD; % total number of props to analyze


for iZ = 1:numberZ
    for iN = 1:numberN
        for iD = 1:numberD
            if     propNo < 10
                disp(['Optimizing propeller   ',num2str(propNo),' of ',num2str(Nprops)])
            elseif propNo < 100
                disp(['Optimizing propeller  ' ,num2str(propNo),' of ',num2str(Nprops)])
            else
                disp(['Optimizing propeller '  ,num2str(propNo),' of ',num2str(Nprops)])
            end
                
            Z = Zall(iZ);
            N = Nall(iN);
            D = Dall(iD);
            
            R  = D/2;                                    
            Js = Vs/((N/60)*D);
            L  = pi/Js;
            
            JS(iZ,iN,iD) = Js;
            
            CTDES = THRUST/(0.5*rho*Vs^2*pi*R^2); % CT thrust coefficient required          
            
            % Scale the inflow velocity profile radii
            RI = ri/R;           


            % ------ Compute the Volumetric Mean Inflow Velocity, Kerwin eqn 163, p.138
            Rhub_oR  = Rhub/R;               % [ ], hub       radius / propeller radius
            XRtemp   = linspace(Rhub_oR,1,100);                           % (temp) radius / propeller radius
            VAItemp  = interp1(RI,VAI,XRtemp,'pchip','extrap');            % (temp) axial inflow velocity /ship velocity
            VMIV(iZ,iN,iD) = 2*trapz(XRtemp,XRtemp.*VAItemp)/(1-Rhub_oR^2);  % [ ], VMIV/ship velocity

            
            % -- Compute cosine spaced vortex & control pt. radii, Kerwin eqn 255, p179
            % RV = radius of vortex  points / propeller radius [ ]
            % RC = radius of control points / propeller radius [ ]
            DEL  = pi/(2*Mp);
            Rdif = 0.5*(1 - Rhub_oR);

            RC = zeros(1,Mp  );
            RV = zeros(1,Mp+1);
            
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

            CoD  = InterpolateChord(XR,XCoD,RC);   % section chord / propeller diameter at ctrl pts



            % --- Do first estimation of tanBi based on 90% of actuator disk efficiency
            EDISK  = 1.8/(1+sqrt(1+CTDES/VMIV(iZ,iN,iD)^2));          % efficiency estimate

            TANBC  = VAC./(pi.*RC./Js + VTC);               % tan(Beta) at control pts.

            VBAC   = VTC.*TANBC./VAC;

            TANBXC = TANBC.*sqrt(VMIV(iZ,iN,iD)./(VAC-VBAC))/EDISK;   % estimation of tan(BetaI)

            % -------------------------- Unload hub and tip as specified by HUF and TUF
            if HUF > 0 || TUF > 0
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
            end

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
                    TMF = 1+(CTDES-CT(iZ,iN,iD))/(5*CTDES);

                elseif CT_iter > 2
                    TMF = TMF_last1 + (TMF_last1-TMF_last2)*(CTDES-CT_last1)/...
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
                

                % --------------------------------------------------------- 
                if any(isnan(G))  | ~isreal(G) 
                    G      = 0*RC';
                    UASTAR = 0*RC;
                    UTSTAR = 0*RC;
                    TANBIC = TANBC;
                    
                    CT_iter = 999;
                    
                    disp('<WARNING>')
                    disp('<WARNING> G == NaN or imaginary... crash avoided...')
                    disp('<WARNING>')
                    disp(' ')
                end
                % --------------------------------------------------------- 

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
                [CT(iZ,iN,iD),CQ(iZ,iN,iD),CP(iZ,iN,iD),KT(iZ,iN,iD),KQ(iZ,iN,iD), CTH(iZ,iN,iD),TAU, Ja(iZ,iN,iD),Jw,VMWV(iZ,iN,iD), EFFYo(iZ,iN,iD), EFFY(iZ,iN,iD),ADEFFY(iZ,iN,iD),QF(iZ,iN,iD)] = ...
                    Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,UADUCT,CD,CoD,G,Z,Js,VMIV(iZ,iN,iD),Hub_flag,Rhub_oR,Rhv,CTD);
      
                % -------------------------------------- Prepare for the next iteration
                CT_iter   = CT_iter + 1;           % iteration in the CT loop
                CT_res    = abs(CT(iZ,iN,iD) - CTDES);       % residual CT
                CT_last2  = CT_last1;              % the CT prior to the last CT
                CT_last1  = CT(iZ,iN,iD);                    % the last value of CT
                TMF_last2 = TMF_last1;             % the TMF prior to the last TMF
                TMF_last1 = TMF;                   % the last value of TMF


            %     disp(['LerbsOptimizer iteration: ',num2str(CT_iter-1),',  max percent error: ',num2str(max(CT_res))]),
            %     disp(' '),    

                if CT_iter > ITER & CT_iter < 999
                    disp('WARNING: Lerbs optimization did NOT converge.'),
                end
            end                                                  % (END WHILE LOOP MF1)
            
            % Expanded Area Ratio
            EAR(iZ,iN,iD) = (2*Z/pi) * trapz(   linspace(Rhub_oR,1,100), InterpolateChord(RC,CoD, linspace(Rhub_oR,1,100) )   );  

            
            propNo = propNo + 1;
        end
    end
end
% -------------------------------------- Pack up output performance metrics
paroutput.part1 = '------   Free parameters   ------';
paroutput.Z     = Zall;
paroutput.N     = Nall;
paroutput.D     = Dall;
paroutput.part2 = '------ Performance metrics (iZ,iN,iD) ------';
paroutput.Js    = JS;
paroutput.KT    = KT;
paroutput.KQ    = KQ;
paroutput.CT    = CT;
paroutput.CQ    = CQ;
paroutput.CP    = CP;
paroutput.CTH   = CTH;
paroutput.EFFYo = EFFYo;              %   open-water   efficiency

paroutput.VMIV  = VMIV;
paroutput.Ja    = Ja;
paroutput.EFFY  = EFFY;
paroutput.ADEFFY= ADEFFY;
paroutput.QF    = QF;

paroutput.VMWV  = VMWV;
paroutput.EAR   = EAR;
% ============================================= END LerbsOptimizer Function
% =========================================================================