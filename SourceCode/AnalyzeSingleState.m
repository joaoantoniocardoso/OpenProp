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
% ============================================= AnalyzeSingleState Function
%
% This function computes the state of a given propeller/turbine
% at a single given off-design tip speed ratio, Lnew.
%
% -------------------------------------------------------------------------
%
%   Inputs:
%       Propeller_flag,Viscous_flag,Hub_flag,Duct_flag                     flags
%       Z,Mp,ITER,Rhv,ALPHAstall,dCLdALPHA,RC,RV,DR,Rhub_oR,CoD,VMIV       propeller properties independent of state
%       TANBIC0,CL0,CD0,                                                   on-design state
%       L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,ALPHA,CL,                   state variables, corresponding to tip speed ratio of L
%
%       Rduct_oR,Cduct_oR,Xduct_oR,XdRING,VARING,GdRING,                   duct properties independent of state
%       DAHIF_times_TANBIC,DRHIF_times_TANBIC,UADIF,CDd,                   duct properties independent of state
%       dCLdALPHAd,DAHIFq_times_TANBIC,DRHIFq_times_TANBIC,                duct properties independent of state
%       CLd0,BetaIDq0,                                                     duct on-design variables
%       Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UADUCT                      duct state variables
%
%       Lnew                                                               tip speed ratio for new state
%
%   NOTE: angle of attack must be given in RADIANS!
%     ALPHA  	[1 x Mp], [rad]
%
%
%   Outputs:
%       L,Js,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,TANBC,ALPHA,CL           state variables, corresponding to tip speed ratio of Lnew
%       Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT        duct state variables
%       KT,KQ,CT,CQ,CP,EFFY,VMWV,TAU,CTD,                                  performance characteristics
%       converged                                                          1 == Newton solver converged, 0 == Newton solver did not converge
%
%   One-line function call: (overwrites current state (L == L) with new state (L == LAMBDA) )
%
% [L,Js,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,TANBC,ALPHA,CL, Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT, KT,KQ,CT,CQ,CP,EFFY,VMWV,TAU,CTD,converged] =  AnalyzeSingleState(Propeller_flag,Viscous_flag,Hub_flag,Duct_flag, Z,Mp,ITER,Rhv,ALPHAstall,dCLdALPHA,RC,RV,DR,Rhub_oR,CoD,VMIV, TANBIC0,CL0,CD0, L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,ALPHA,CL, Rduct_oR,Cduct_oR,Xduct_oR,XdRING,VARING,GdRING, DAHIF_times_TANBIC,DRHIF_times_TANBIC,UADIF,CDd, dCLdALPHAd,DAHIFq_times_TANBIC,DRHIFq_times_TANBIC, CLd0,BetaIDq0, Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UADUCT, Lnew);
%
% -------------------------------------------------------------------------

function [L,Js,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,TANBC,ALPHA,CL,                   ...   
          Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT,                              ...
          KT,KQ,CT,CQ,CP,EFFY,VMWV,TAU,CTD,converged] =                 ...
                                                                                      ...
      AnalyzeSingleState(Propeller_flag,Viscous_flag,Hub_flag,Duct_flag,              ...
                         Z,Mp,ITER,Rhv,ALPHAstall,dCLdALPHA,RC,RV,DR,Rhub_oR,CoD,VMIV,...
                         TANBIC0,CL0,CD0,                                             ...
                         L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,ALPHA,CL,             ...
                                                                                      ... 
                         Rduct_oR,Cduct_oR,Xduct_oR,XdRING,VARING,GdRING,             ...
                         DAHIF_times_TANBIC,DRHIF_times_TANBIC,UADIF,CDd,             ...
                         dCLdALPHAd,DAHIFq_times_TANBIC,DRHIFq_times_TANBIC,          ...
                         CLd0,BetaIDq0,          ...
                         Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UADUCT, ...
                         Lnew)
    
Wake_flag = 0;

% -------------------------------------------------------------------------   
% -------------------------------------------------------------------------  
% Set the tip-speed ratio to the given new value, Lnew
L  = Lnew;
Js = pi/L;
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------   

% ------------------------------------------------ Implement RepairSpline.m
% To smooth data X(RC):  X_smooth = X*Bsmooth;
Bsmooth = RepairSplineMatrix(RC);
% -------------------------------------------------------------------------

%%

TANBC = (VAC)./(L*RC + VTC); % inflow angle 

% ---------------------------------------------------- inflow angle 
BetaIC0 = atan(TANBIC0);    % [rad]
BetaIC  = atan(TANBIC );    % [rad]

       
% ------------------------- Initialize vortex Horseshoe Influence Functions 
[UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBIC,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);
  
%% 
    feather     = 0;
    deltaALPHA  = 0*RC;
    Crash_count = 0;
    
    
    % =====================================================================
    %
    % COPY AND PASTE OF THE CODE IN Analyze.m NEWTON SOLVER
    %
    % =====================================================================
    % ==   FIND STATE OF SYSTEM USING NEWTON SOLVER:
    % ==                   iterate to solve residual equations for unknowns
    N_iter   = 1;                            
    ERROR    = 1;
    ERRORtol = 0.005;
    relax    = 0.9;
    converged = 1;  % 1 == state converged, 0 == state did not converge

    while N_iter <= ITER & any(abs(ERROR) > ERRORtol)          % (WHILE LOOP N1)
        % disp(['------- Newton iteration: ',num2str(N_iter)])  % status message
        % disp(' '),

        % ---------------------------------- Store last state of the system
         VSTARlast =  VSTAR;
         alphalast =  ALPHA;
            CLlast =     CL;
             Glast =      G;

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
            

        % ------------- Solve Newton problem to update entire blade at once
        for m = 1:Mp  
            % ------------------------------------------------------ Evaluate residuals
            R(       m) =  VSTAR(m) - sqrt((VAC(m)+UADUCT(m)+UASTAR(m)).^2 + (L*RC(m)+VTC(m)+UTSTAR(m)).^2); % R1

            R(1*Mp + m) =  ALPHA(m) - (BetaIC0(m) - atan(TANBIC(m)) + feather);                         % R2
            
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

        end   % ------------------------------------ END for each blade section

        
        % ------------------------- Update Newton solver vector of unknowns
        DX     = linsolve(J,-R);

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


        % ----------------------------------------------------------------- 
        % Update induced velocities (influence of propeller on propeller)
        UASTAR = (UAHIF*G)';  
        UTSTAR = (UTHIF*G)';  
        
        TANBIC = (VAC + UADUCT + UASTAR)./(L*RC + VTC + UTSTAR);
                
        % Smooth the inflow angle for numerical stability:
        TANBICsmooth = TANBIC * Bsmooth;
        % ----------------------------------------------------------------- 
        
        
        if Wake_flag == 0 

            [UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBICsmooth,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);

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
        
        
        % ------------------------------------------------ End update variables
        
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
        

        % ------------------------------------------ Evaluate percent error
        PERROR    = [ VSTARlast -  VSTAR, ...
                     alphalast -  ALPHA, ...
                        CLlast -     CL, ...
                        Glast' -     G'] ./ ...
                   [max(abs(VSTAR ),1e-4), ...
                    max(abs(ALPHA ),1e-2), ...
                    max(abs(CL    ),1e-4), ...
                    max(abs(G'    ),1e-6)];
        % ---------------------------------------------- END evaluate percent error


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
    
    % =================================================== END NEWTON SOLVER  
    % =====================================================================    
%%
    
    
    % ------------------------------------------------------ Compute forces
    
    % ---------------------------------------------------------------------
    if Duct_flag == 1
        DAHIF = 0*DAHIF_times_TANBIC;
        DRHIF = 0*DRHIF_times_TANBIC;
        
        % Update the influence of propeller on duct
        for m = 1:Mp;
            DAHIF(:,m) = DAHIF_times_TANBIC(:,m) / TANBIC(m);
            DRHIF(:,m) = DRHIF_times_TANBIC(:,m) / TANBIC(m);
        end

        % Update induced velocities at the duct (influence of propeller on duct)
        UARING = (DAHIF*G)';  
        URRING = (DRHIF*G)'; 

        % Find duct thrust, (CTDdes==1.0 is unused here)            
        [CTD,junk] = Duct_Thrust(XdRING,Rduct_oR,VARING,UARING,URRING,GdRING,Gd,CDd,1.0); 
    else 
        UARING = 0*RC;
        URRING = 0*RC;
        CTD = 0;
    end
    % ---------------------------------------------------------------------        
    
    [CL,CD] = CLCD_vs_ALPHA(ALPHA,ALPHAstall,CL0,CD0,dCLdALPHA,Propeller_flag);
    
    if Viscous_flag == 0
        CD = 0*CD;
    end

    [CT,CQ,CP,KT,KQ, CTH,TAU, Ja,Jw,VMWV, EFFYo, EFFY,ADEFFY,QF] = ...
              Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,UADUCT,CD,CoD,G,Z,Js,VMIV,Hub_flag,Rhub_oR,Rhv,CTD);

    
    % -------------------------------------------------------------------------
    disp(' '),
    disp(['Forces for state:  L = ',num2str(L),' , Js = ',num2str(Js)]),
    
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

    
        if N_iter-1 < 10
            disp(['The max error for iteration  ',num2str(N_iter-1),' is: ',num2str(max(abs(ERROR)))]),  
        else
            disp(['The max error for iteration ' ,num2str(N_iter-1),' is: ',num2str(max(abs(ERROR)))]),  
        end
        disp(' '),    
    
    if converged == 0
        disp('WARNING: Newton solver did NOT converge.'),
        disp(' ')
        disp('------------------------------------------------'),
        disp(' ')
    else    
        disp('------------------------------------------------'),
        disp(' ')
    end
        
end                             
