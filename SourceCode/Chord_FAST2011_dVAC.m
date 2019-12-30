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

% -------------------------------------------------------------------------
% 11/18/2011 Brenden Epps
%
% This function optimizes blade chord and thickness using the method of
% (Epps et al., FAST 2011), reconfigured to optimize given the design point
% and axial inflow deviation dVAC, as if the blade were in the ship's wake.
%
% This implementation enforces the ABS thickness rule.
%
% -------------------------------------------------------------------------


function [CoD, t0oc, C_res] = Chord_FAST2011_dVAC(dVAC,SIGMAs,L,RC,VAC,VTC,UADUCT,UASTAR,UTSTAR,G,CoD,t0oc,...
                                           DR,CD,Z,Js,VMIV,Hub_flag,Rhub_oR,Rhv,CTD,Mp,R,rho,Vs,N)

                                       
% -------------------------------------------------------------------------
CoD_last = CoD;    % last value of CoD
% -------------------------------------------------------------------------

disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Performing FAST2011 chord optimization...')
disp(' ')                                       
                                       
% -------------------------------------------------------------------------
% Parameters:
NACAa       = 0.8;               % 'NACA a=0.8' meanline
RHOLE       = 0.640279919559199; % == 0.00639/(2*0.04995)^2 == leading edge radius / chord, for t0oc==1 (RLE = RHOLE*t0oc^2)
CAVmargin   = 1.0;               % allowable SIGMAs = CAVmargin*SIGMAs
CoDmax      = 0.65;              % maximum allowable c/D
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
dALPHA = atan( (VAC + UADUCT + UASTAR)./(L*RC + VTC + UTSTAR) )   -   atan( (VAC + UADUCT + UASTAR - dVAC)./(L*RC + VTC + UTSTAR) );  % [rad] inflow angle variation    

  VSTAR  = sqrt( (VAC + UADUCT + UASTAR       ).^2 + (L*RC + VTC + UTSTAR).^2 );
% VSTARe = sqrt( (VAC + UADUCT + UASTAR - dVAC).^2 + (L*RC + VTC + UTSTAR).^2 );

  SIGMA  = SIGMAs./VSTAR.^2;  % local cavitation number
% SIGMA  = SIGMAs./VSTARe.^2; % local cavitation number
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% CASE 2: Find optimum using Newton solver
g1   = 4/pi;
g2   = pi*G'./((1+NACAa)*VSTAR);
g3sq = (2/RHOLE) .* (dALPHA).^2;

for i = 1:Mp
    [CoD_2(i),t0oc(i)] = Chord_ctSolver_deltaALPHA(CAVmargin*SIGMA(i), g1, g2(i), g3sq(i) , CoD(i), t0oc(i));
end

% -----------------------------------------------------------------
% CASE 1: Solve for CoD required by CASE 1 for given t0oc calculated using CASE 2.
%
CoD_1 = g2 ./ (-dALPHA + sqrt(dALPHA.^2 + (CAVmargin*SIGMA+1)) - (1 + g1*t0oc) );
% -----------------------------------------------------------------


% -----------------------------------------------------------------
% Set CoD to maximum of CASE 1 or CASE 2
%
CoD  = max(CoD_1,CoD_2);

% Enforce maximum allowable c/D 
CoD = min(CoD,CoDmax);


t0oD = t0oc .* CoD;
% -----------------------------------------------------------------
% -----------------------------------------------------------------
% -----------------------------------------------------------------


% -----------------------------------------------------------------
% Augment CoD and t0oD using ABS_Blade method (where "rated speed" is the endurance speed, Vs)
% -----------------------------------------------------------------
% alphaItilde = 1.54; % [deg]
%    CLItilde = 1;

CL       = 2*pi*G'./(VSTAR.*CoD);           % lift coefficient

alphaI   = 1.54*CL;                         % == alphaItilde.*CL/CLItilde;        

TANBIC   = (VAC + UADUCT + UASTAR)./(L*RC + VTC + UTSTAR);

TANtheta = tand( atand(TANBIC) + alphaI );  % tan(pitch angle)

% CP == power coefficient
[junk1,junk2,CP] = Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,UADUCT,CD,CoD,G,Z,Js,VMIV,Hub_flag,Rhub_oR,Rhv,CTD);


[CoD, t0oD] = Chord_ABS_Blade(RC,t0oD,t0oc,TANtheta,  CP,R,rho,Vs,N,Z);          
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------        
C_res = max( abs(CoD-CoD_last) ); % residual CoD
% -------------------------------------------------------------------------
      
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp(['The max  C_res is: ',num2str(C_res)]),
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Continuing circulation optimization:')
disp(' ')


end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%%
% -------------------------------------------------------------------------
% Find optimum c and tau using Newton solver -- given deltaALPHA
% -------------------------------------------------------------------------
function [c,tau] = Chord_ctSolver_deltaALPHA(SIGMA, g1, g2, g3sq, c, tau)

if	c < 0.2
    c = 0.2;
end

if  tau < 0.02
    tau = 0.02;
elseif tau > 0.1;
    tau = 0.1;
end

relax   = 0.9;

J = zeros(2,2);
R =  ones(2,1);

iter = 0;

while (abs(R(1)) > 1e-6 || abs(R(2)) > 1e-6) && iter < 1000
    % given: SIGMA, g1, g2, g3sq, 

    % Residual equations:
    R1      =      (1 + g1*tau + g2/c)^2 +   g3sq/(tau^2) - 1 - SIGMA; % => 0

    R2      = 2*g1*(1 + g1*tau + g2/c)   - 2*g3sq/(tau^3);             % => 0

    R       = [R1;R2];

    % Jacobian == J = [dR1dc dR1dt; dR2dc dR2dt]

    dR1dc   = 2 * (1 + g1*tau + g2/c) * (-g2 / c^2);

    dR1dt   = R2;

    dR2dc   = -2 * g1 * g2 / c^2;

    dR2dt   =  2 * g1^2            + 6 * g3sq / tau^4;



    J       = [dR1dc dR1dt; ...
               dR2dc dR2dt];

    % dx == delta[c; tau]

    dx      = linsolve(J, -R);

    % enforce a maximum step size
    dx(1) = sign(dx(1)) * min([abs(dx(1)),0.05]);
    dx(2) = sign(dx(2)) * min([abs(dx(2)),0.02]);

    % update the state variables
    if       (c  + relax * dx(1)) > 0
        c  =  c  + relax * dx(1);
    else
        c  =  c/2;
    end

    if         (tau  + relax * dx(2)) > 0
        tau  =  tau  + relax * dx(2);
    else
        tau  = tau/2;
    end

    iter = iter + 1;
end


end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------




%%
% -----------------------------------------------------------------

%             % =============================================================
%             % Plot contours of (-CPmax) for each blade section -- Case 1
%             % =============================================================
%             tt = [0.001:0.001:0.4];
%             cc = [0.001:0.001:1];
%             [t,c] = ndgrid(tt,cc);
% 
%             close all,
%             
%             for i = 1:Mp
% 
% %                 % Case 2:
% %                 v = (1 + g1*t + g2(i)./c).^2 + g3sq(i)./(t.^2) - 1;
% 
%                 % Case 1:
%                 v = (1 + g1*t + g2(i)./c).^2 + 2*dALPHA(i)*(1 + g1*t + g2(i)./c) - 1;
% 
%                 
%                 fig; axis([0 1 0 0.4])
%                 
%                 contourf(c,t,v,[-0.1:0.1:2])
%                 H1 = contour(c,t,v,SIGMA(i)*[1 1],'color','r','linewidth',2); 
% 
%                 % Case 1 vs case 2 cross-over point:  Case 2 if CoD is greater than case2_if_CoD_greater_than_this
%                 case2_if_CoD_greater_than_this = g2(i)./( (dALPHA(i)/RHOLE)./tt.^2 - (1 + g1*tt) );
%                 case2_if_CoD_greater_than_this(case2_if_CoD_greater_than_this < 0) = NaN;
%           plot( case2_if_CoD_greater_than_this ,  tt,'m--','linewidth',2)
%                 
% %                 % Case 2 constraint: optimum CoD vs tau satisfies this constraint:  2*g1*(1 + g1*tau + g2/c) - 2*g3sq/(tau^3) = 0
% %                 CoD_case2 = g2(i)./((2*g3sq(i)./tt.^3)/(2*g1) - (1 + g1*tt));
% %                 CoD_case2(CoD_case2 < 0) = NaN;
% %           plot( CoD_case2 ,  tt,'y--','linewidth',2)
% 
%                
% %                 % Case 1: optimum CoD vs tau (for which -CPmax == CAVmargin*SIGMA)
% %                 CoD_case1 = g2(i) ./ (-dALPHA(i) + sqrt(dALPHA(i).^2 + (CAVmargin*SIGMA(i)+1)) - (1 + g1*tt) );
% %                 CoD_case1(CoD_case1 < 0) = NaN;
% %           plot( CoD_case1 , tt,'r--','linewidth',2)
% 
%           
% 
%               % plot(CoD81e(i),t0oc(i),'yx','markersize',18)
% %                 plot(CoD82e(i),t0oc(i),'g.','markersize',18)
% 
%                 xlabel('c/D'), ylabel('t0/c'),% colorbar
% 
%                 
%                 
%                 % Formatting                                
%                 set(gca,'XTick',[0:0.2:1])
%                 set(gca,'YTick',[0:0.05:0.4])
%             
%                 set(gca,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'})
%                 set(gca,'YTickLabel',{'0','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40'})
%             
%                 caxis([0 2]);  % h = colorbar; 
%             
%                 axis([0 1 0 0.4])
%                 box on
%             
%                 set(gca,'FontName','Times','FontSize',18)
%                 xlabel('c / D','FontName','Times','FontSize',18)
%                 ylabel('t0 / c','FontName','Times','FontSize',18)
%             
% %                 if i < 10
% %                     saveas(gcf,['fig_CPmax_case1_0',num2str(i)],'epsc')
% %                 else
% %                     saveas(gcf,['fig_CPmax_case1_',num2str(i)],'epsc')
% %                 end
%                 
%             end
%             %%
%             % =============================================================
%             % Plot contours of (-CPmax) for each blade section -- Case 2
%             % =============================================================
%             tt = [0.001:0.001:0.4];
%             cc = [0.001:0.001:1];
%             [t,c] = ndgrid(tt,cc);
% 
%             close all,
%             
%             for i = 1:Mp
% 
%                 % Case 2:
%                 v = (1 + g1*t + g2(i)./c).^2 + g3sq(i)./(t.^2) - 1;
% 
% %                 % Case 1:
% %                 v = (1 + g1*t + g2(i)./c).^2 + 2*dALPHA(i)*(1 + g1*t + g2(i)./c) - 1;
% 
%                 
%                 fig; 
%                 contourf(c,t,v,[-0.1:0.1:2])
%                 
%                 H1 = contour(c,t,v,SIGMA(i)*[1 1],'color','r','linewidth',2); 
% 
%                 
%                 % Case 1 vs case 2 cross-over point:  Case 2 if CoD is greater than case2_if_CoD_greater_than_this
%                 case2_if_CoD_greater_than_this = g2(i)./( (dALPHA(i)/RHOLE)./tt.^2 - (1 + g1*tt) );
%                 case2_if_CoD_greater_than_this(case2_if_CoD_greater_than_this < 0) = NaN;
%           plot( case2_if_CoD_greater_than_this ,  tt,'m--','linewidth',2)
%                 
%                 % Case 2 constraint: optimum CoD vs tau satisfies this constraint:  2*g1*(1 + g1*tau + g2/c) - 2*g3sq/(tau^3) = 0
%                 CoD_case2 = g2(i)./((2*g3sq(i)./tt.^3)/(2*g1) - (1 + g1*tt));
%                 CoD_case2(CoD_case2 < 0) = NaN;
%           plot( CoD_case2 ,  tt,'y--','linewidth',2)
% 
%                
%                 % Case 1: optimum CoD vs tau (for which -CPmax == CAVmargin*SIGMA)
%                 CoD_case1 = g2(i) ./ (-dALPHA(i) + sqrt(dALPHA(i).^2 + (CAVmargin*SIGMA(i)+1)) - (1 + g1*tt) );
%                 CoD_case1(CoD_case1 < 0) = NaN;
%           plot( CoD_case1 , tt,'r--','linewidth',2)
% 
%           
% 
%               % plot(CoD81e(i),t0oc(i),'yx','markersize',18)
%                 plot(CoD82e(i),t0oc(i),'g.','markersize',18)
% 
%                 xlabel('c/D'), ylabel('t0/c'),% colorbar
% 
%                 
%                 
%                 % Formatting                                
%                 set(gca,'XTick',[0:0.2:1])
%                 set(gca,'YTick',[0:0.05:0.4])
%             
%                 set(gca,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'})
%                 set(gca,'YTickLabel',{'0','0.05','0.10','0.15','0.20','0.25','0.30','0.35','0.40'})
%             
%                 caxis([0 2]);  % h = colorbar; 
%             
%                 axis([0 1 0 0.4])
%                 box on
%             
%                 set(gca,'FontName','Times','FontSize',18)
%                 xlabel('c / D','FontName','Times','FontSize',18)
%                 ylabel('t0 / c','FontName','Times','FontSize',18)
%             
% %                 if i < 10
% %                     saveas(gcf,['fig_CPmax_case2_0',num2str(i)],'epsc')
% %                 else
% %                     saveas(gcf,['fig_CPmax_case2_',num2str(i)],'epsc')
% %                 end
%                 
%             end            
%             % =============================================================
%             % =============================================================
            
%%            
