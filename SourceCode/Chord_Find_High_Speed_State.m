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
% =========================================== Chord_Find_High_Speed_State.m
%  Created:  2/14/2011 Oscar Viquez 
% Modified: 11/18/2011 Brenden Epps
% -------------------------------------------------------------------------

function [Gh, VSTARh, Jh, UASTARh,UTSTARh,CDh] = Chord_Find_High_Speed_State(CTreq, J0,...
                         ...       
                         Propeller_flag,Viscous_flag,Hub_flag,Duct_flag,              ...
                         Z,Mp,ITER,Rhv,ALPHAstall,dCLdALPHA,RC,RV,DR,Rhub_oR,CoD,VMIV,...
                         TANBIC0,CL0,CD0,                                             ...
                         L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,ALPHA,CL);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
if Duct_flag == 1
    disp('ERROR: Find_High_Speed_State.m is not configured for the ducted case...please update the code...')
    return
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

    % Duct on-design variables
    BetaIDq0             = BetaIDq;    
    CLd0                 = CLd;    
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
              
                     
                     
                     
                            
% -------------------------------------------------------------------------
% One-line function call: (overwrites current state (L == L) with new state (L == Lnew) )
Lnew = pi/J0;
                            
[L,Js,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,TANBC,ALPHA,CL, Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT, KT,KQ,CT,CQ,CP,EFFY,VMWV,TAU,CTD,converged] =  AnalyzeSingleState(Propeller_flag,Viscous_flag,Hub_flag,Duct_flag, Z,Mp,ITER,Rhv,ALPHAstall,dCLdALPHA,RC,RV,DR,Rhub_oR,CoD,VMIV, TANBIC0,CL0,CD0, L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,ALPHA,CL, Rduct_oR,Cduct_oR,Xduct_oR,XdRING,VARING,GdRING, DAHIF_times_TANBIC,DRHIF_times_TANBIC,UADIF,CDd, dCLdALPHAd,DAHIFq_times_TANBIC,DRHIFq_times_TANBIC, CLd0,BetaIDq0, Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UADUCT, Lnew);

KT0 = KT;

if converged == 0
    % Output variables
    Gh      = 0*G;
    VSTARh  = 0*VSTAR;
    Jh      = 0*Js;
    UASTARh = 0*UASTAR;
    UTSTARh = 0*UTSTAR;
    CDh     = 0*G;
    return
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% One-line function call: (overwrites current state (L == L) with new state (L == Lnew) )
J1 = J0 - 0.1;

Lnew = pi/J1;
                            
[L,Js,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,TANBC,ALPHA,CL, Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT, KT,KQ,CT,CQ,CP,EFFY,VMWV,TAU,CTD,converged] =  AnalyzeSingleState(Propeller_flag,Viscous_flag,Hub_flag,Duct_flag, Z,Mp,ITER,Rhv,ALPHAstall,dCLdALPHA,RC,RV,DR,Rhub_oR,CoD,VMIV, TANBIC0,CL0,CD0, L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,ALPHA,CL, Rduct_oR,Cduct_oR,Xduct_oR,XdRING,VARING,GdRING, DAHIF_times_TANBIC,DRHIF_times_TANBIC,UADIF,CDd, dCLdALPHAd,DAHIFq_times_TANBIC,DRHIFq_times_TANBIC, CLd0,BetaIDq0, Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UADUCT, Lnew);

KT1   = KT;
CTNEW = CT;
JsNEW = J1;

if converged == 0
    % Output variables
    Gh      = 0*G;
    VSTARh  = 0*VSTAR;
    Jh      = 0*Js;
    UASTARh = 0*UASTAR;
    UTSTARh = 0*UTSTAR;
    CDh     = 0*G;
    return
end
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
dKTdJs = (KT1-KT0)/(J1-J0);  % finite difference approximation
JsOLD  = 0.5*( J0 +  J1);
KTOLD  = 0.5*(KT0 + KT1);
% -------------------------------------------------------------------------


%%
% Iterate
while abs(CTNEW - CTreq) > 1e-3
    % disp('-----')

    JsNEW = ( dKTdJs + sqrt(dKTdJs^2 - 4*(CTreq*pi/8)*(dKTdJs*JsOLD-KTOLD)) ) / (2 * CTreq*pi/8);
   
    
    % ---------------------------------------------------------------------
    Lnew = pi/JsNEW;
    
    [L,Js,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,TANBC,ALPHA,CL, Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UARING,URRING,UADUCT, KT,KQ,CT,CQ,CP,EFFY,VMWV,TAU,CTD,converged] =  AnalyzeSingleState(Propeller_flag,Viscous_flag,Hub_flag,Duct_flag, Z,Mp,ITER,Rhv,ALPHAstall,dCLdALPHA,RC,RV,DR,Rhub_oR,CoD,VMIV, TANBIC0,CL0,CD0, L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,ALPHA,CL, Rduct_oR,Cduct_oR,Xduct_oR,XdRING,VARING,GdRING, DAHIF_times_TANBIC,DRHIF_times_TANBIC,UADIF,CDd, dCLdALPHAd,DAHIFq_times_TANBIC,DRHIFq_times_TANBIC, CLd0,BetaIDq0, Gd,UARINGq,URRINGq,VSRINGq,BetaIDq,CLd,UADUCT, Lnew);

    KTNEW = KT;   
    CTNEW = CT;
    
    if converged == 0
        % Output variables
        Gh      = 0*G;
        VSTARh  = 0*VSTAR;
        Jh      = 0*Js;
        UASTARh = 0*UASTAR;
        UTSTARh = 0*UTSTAR;
        CDh     = 0*G;
        return
    end
    % ---------------------------------------------------------------------
    
    dKTdJs = (KTNEW-KTOLD)/(JsNEW-JsOLD);
    JsOLD  = 0.5*(JsNEW+JsOLD);
    KTOLD  = 0.5*(KTNEW+KTOLD);

    disp(['CT req = ',num2str(CTreq), '   CT new = ',num2str(CTNEW)])
%     disp(['CT new = ',num2str(CTNEW)])
%     disp(['Js new = ',num2str(JsNEW)])

    disp(' ')
    disp('------------------------------------------------'),
    disp(' ')    
end


% Output variables
Gh      = G;
VSTARh  = VSTAR;
Jh      = JsNEW;
UASTARh = UASTAR;
UTSTARh = UTSTAR;

[CLh,CDh] = CLCD_vs_ALPHA(ALPHA,ALPHAstall,CL0,CD0,dCLdALPHA,Propeller_flag);



end