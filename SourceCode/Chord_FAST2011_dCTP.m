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
% (Epps et al., FAST 2011), given the design point "endurance speed" 
% and a "high speed" state.
%
% This implementation enforces the ABS thickness rule.
%
% -------------------------------------------------------------------------


function [CoD, t0oc, C_res] = Chord_FAST2011_dCTP(SIGMAh,CTDESh, Jh,CoD,t0oc,...      
                         Propeller_flag,Viscous_flag,Hub_flag,Duct_flag, ...
                         Z,Mp,ITER,Rhv,RC,RV,DR,Rhub_oR,VMIV,CTD,...
                         L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,CD,...
                         R,rho,Vs,N);
                     
if Duct_flag == 1
    disp('ERROR: Chord_FAST2011_dCTP.m is not configured for the ducted case...please update the code...')
    return
else 
    UADUCT = 0*RC;
end
   
% -------------------------------------------------------------------------
CoD_last = CoD;    % last value of CoD
% -------------------------------------------------------------------------

disp(' ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Performing FAST2011 chord optimization...')
disp(' ')

% -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
% Find off-design state (Gh, VSTARh, Jh) - h == high-speed state
% -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
%
% Compute design state variables
ALPHA   = 0*RC;    
CL      = 2*pi*G'./(VSTAR.*CoD); % [ ] lift coefficient
Js      = pi/L;

TANBIC0 = TANBIC;
CL0     = CL;
CD0     = CD;


% Off-design analysis parameters
ALPHAstall  = 8*pi/180 * ones(size(RC));

% -------------------------------------------------------------------------
% Propeller Aspect Ratio (PAR)
PAR = (1-Rhub_oR)^2 / trapz( linspace(Rhub_oR,1,100) , InterpolateChord(RC,CoD,linspace(Rhub_oR,1,100))  );

% Lift curve slope:
dCLdALPHA = (2*pi / (1 + 2/PAR)) * ones(size(RC));
% -------------------------------------------------------------------------

% Compute off-design state                     
[Gh, VSTARh, Jh, UASTARh,UTSTARh,CDh] = Chord_Find_High_Speed_State(CTDESh, Jh,...
                         ...       
                         Propeller_flag,Viscous_flag,Hub_flag,Duct_flag,              ...
                         Z,Mp,ITER,Rhv,ALPHAstall,dCLdALPHA,RC,RV,DR,Rhub_oR,CoD,VMIV,...
                         TANBIC0,CL0,CD0,                                             ...
                         L,G,VAC,VTC,UASTAR,UTSTAR,VSTAR,TANBIC,ALPHA,CL);                     
       
% If Find_High_Speed_State crashed...     
if Jh == 0
    CoD   = 0*CoD;
    t0oc  = 0*t0oc;
    C_res = 0;    
    
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('Aborting FAST2011 chord optimization... (c/D = 0, t0/c = 0)'),
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(' ')
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('Continuing circulation optimization:')
    disp(' ')
    
    return
end
% -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
% -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

% -!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
% UPDATE CHORD and THICKNESS DISTRIBUTIONS -- CoD and t0oc
% -!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!- 

% -------------------------------------------------------------------------
% Parameters:
NACAa       = 0.8;               % 'NACA a=0.8' meanline
RHOLE       = 0.640279919559199; % == 0.00639/(2*0.04995)^2 == leading edge radius / chord, for t0oc==1 (RLE = RHOLE*t0oc^2)
CAVmargin   = 1.0;               % allowable SIGMAs = CAVmargin*SIGMAs
CoDmax      = 0.65;              % maximum allowable c/D
% -------------------------------------------------------------------------

SIGMA     = SIGMAh./VSTARh.^2;      % local cavitation number

% -------------------------------------------------------------------------
% CASE 2: Find optimum using Newton solver
g1   = 4/pi;
g2   = pi*G'./((1+NACAa)*VSTAR);
g3   =              (Gh'./VSTARh - G'./VSTAR);
g3sq = (2/RHOLE) .* (Gh'./VSTARh - G'./VSTAR).^2;  % not exactly g3^2

for i = 1:Mp
    [CoD_2(i),t0oc(i)] = Chord_ctSolver(CAVmargin*SIGMA(i), g1, g2(i), g3sq(i) , CoD(i), t0oc(i));
end


% -----------------------------------------------------------------
% CASE 1:
a1   = 4/pi;
a2   = pi*G'./((1+NACAa)*VSTAR);
a3   = (Gh'./VSTARh - G'./VSTAR);

% Solve for CoD required by CASE 1 for given t0oc calculated using CASE 2 or 3.
%          (1+CAVmargin*SIGMA-(1+a1*t0oc).^2)
%         -(2*(1+a1*t0oc).*(a2+a3))
%         -(a2.^2+2*a2.*a3)
%
CoD_1 = ( (2*(1+a1*t0oc).*(a2+a3)) + sqrt( (2*(1+a1*t0oc).*(a2+a3)).^2 + 4*(1+CAVmargin*SIGMA-(1+a1*t0oc).^2).*(a2.^2+2*a2.*a3) ) )./( 2*(1+CAVmargin*SIGMA-(1+a1*t0oc).^2) );


% -----------------------------------------------------------------


% Set CoD to maximum of CASE 1 or CASE 2
CoD = max(CoD_1,CoD_2);

% Enforce maximum allowable c/D 
CoD = min(CoD,CoDmax);

t0oD = t0oc.*CoD;
% -----------------------------------------------------------------
% -----------------------------------------------------------------
% -----------------------------------------------------------------

% save temp
% disp('paused...saved temp')
% pause,

% -----------------------------------------------------------------
% Augment CoD and t0oD using ABS_Blade method (where "rated speed" is the endurance speed, Vs)
% -----------------------------------------------------------------
% alphaItilde = 1.54; % [deg]
%    CLItilde = 1;

CL       = 2*pi*G'./(VSTAR.*CoD);           % lift coefficient

alphaI   = 1.54*CL;                         % == alphaItilde.*CL/CLItilde;        

TANtheta = tand( atand(TANBIC) + alphaI );  % tan(pitch angle)


% CP == power coefficient
[junk1,junk2,CP] = Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,UADUCT,CD,CoD,G,Z,Js,VMIV,Hub_flag,Rhub_oR,Rhv,CTD);


% fig(1234); plot(RC,CoD,'k')
% fig(1235); plot(RC,t0oD,'k')

[CoD, t0oD] = Chord_ABS_Blade(RC,t0oD,t0oc,TANtheta,  CP,R,rho,Vs,N,Z);  


% fig(1234); plot(RC,CoD,'r')
% fig(1235); plot(RC,t0oD,'r')


% -------------------------------------------------------------------------
% 12/16/2011 "Design Study 27" Given new thickness distribution, find chord distribution following [-CP]max = SIGMA contour
CoD_2 = (g1.*t0oD + g2) ./ ( sqrt(1 + CAVmargin*SIGMA - g3sq./t0oD.^2) - 1 );

CoD_1 = (g1.*t0oD + g2 + g3)./(CAVmargin*SIGMA) + sqrt( ((g1.*t0oD + g2 + g3)./(CAVmargin*SIGMA)).^2 + (2*g3.*(g1.*t0oD + g2) + (g1.*t0oD + g2).^2)./(CAVmargin*SIGMA) );

t0oc_1 = t0oD./CoD_1;
t0oc_2 = t0oD./CoD_2;


% Set CoD to maximum of CASE 1 or CASE 2
CoD = max(CoD_1,CoD_2);

% Enforce maximum allowable c/D 
CoD = min(CoD,CoDmax);

t0oc = t0oD./CoD;

% fig(1234); 
%     plot(RC,CoD_1,'m')
%     plot(RC,CoD_2,'c')
%     plot(RC,CoD,'g')
% 
% % for i = 1:Mp
% %     fig(i);
% %     plot(CoD_1(i),t0oc_1(i),'.m','markersize',20)
% %     plot(CoD_2(i),t0oc_2(i),'.c','markersize',20)
% % end

% -------------------------------------------------------------------------


% disp('paused...')
% pause,
% close(1234),
% close(1235),
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


        % % -----------------------------------------------------------------
        % % Augment CoD and t0oD using ABS_Blade method (where "rated speed" is the max speed, Vh)
        % % -----------------------------------------------------------------
        % % alphaItilde = 1.54; % [deg]
        % %    CLItilde = 1;
        % 
        % CL       = 2*pi*G'./(VSTAR.*CoD);           % lift coefficient
        % 
        % alphaI   = 1.54*CL;                         % == alphaItilde.*CL/CLItilde;        
        % 
        % TANtheta = tand( atand(TANBIC) + alphaI );  % tan(pitch angle)
        % 
        % % CP == power coefficient
        % [junk1,junk2,CPh] = Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,UADUCT,CDh,CoD,Gh,Z,Js,VMIV,Hub_flag,Rhub_oR,Rhv,CTD);
        % 
        % 
        % [CoD, t0oD] = Chord_ABS_Blade(RC,t0oD,t0oc,TANtheta,  CPh,R,rho,Vs,N,Z);          
        % % -------------------------------------------------------------------------
        % % -------------------------------------------------------------------------

        


        
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

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% 



%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [c,tau] = Chord_ctSolver(SIGMA, g1, g2, g3sq, c, tau)

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

while abs(R(1)) > 1e-6 || abs(R(2)) > 1e-6
    % given: SIGMA, g1, g2, g3sq, 

    % Residual equations:
    R1      =      (1 + g1*tau + g2/c)^2 +   g3sq/(c^2*tau^2) - 1 - SIGMA; % => 0

    R2      = 2*g1*(1 + g1*tau + g2/c)   - 2*g3sq/(c^2*tau^3);             % => 0

    R       = [R1;R2];

    % J = Jacobian = [dR1dc dR1dt; dR2dc dR2dt]

    dR1dc   = 2 * (1 + g1 * tau + g2 / c) * (-g2 / c^2) - 2 * g3sq / (c^3 * tau^2);

    dR1dt   = 2 * (1 + g1 * tau + g2 / c) * ( g1      ) - 2 * g3sq / (c^2 * tau^3);

    dR2dc   = -2 * g1 * g2 / (c^2) + 4 * g3sq / (c^3 * tau^3);

    dR2dt   =  2 * g1^2            + 6 * g3sq / (c^2 * tau^4);



    J       = [dR1dc dR1dt; ...
               dR2dc dR2dt];

    % dx = delta[c; tau]

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


end


end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
%         
% 
% %%                       
%          
%         % -----------------------------------------------------------------
%         % Hack the code here in order to make CP contour map figures:         
% %         for to_make_code_fold = 1:1        
%   
% 
%             % % -------------------------------------------------------------------------
%             % % -------------------------------------------------------------------------
%             % % FIGURE:  t0oD distribution showing ABS_Blade modification      ----------
%             % % ------------------------------------------------------------------------- 
%             % RG = linspace(Rhub_oR,1,100);
%             % 
%             % fig;
%             %     
%             % % Initial thickness distribution: t0oD == t0oc.*CoD_2
%             %     H0 = plot(RG,interp1(RC,t0oc,RG,'spline','extrap').*InterpolateChord(RC,CoD_2,RG),'r--'); 
%             %     
%             %     
%             % % Modified thickness distribution: t0oD == t0oc.*CoD,   with same t0oc
%             %     H1 = plot(RG,interp1(RC,t0oc,RG,'spline','extrap').*InterpolateChord(RC,CoD,RG),'k');
%             %     
%             %     
%             % % Root thickness    
%             %     t0oD25 = 0.0258;    
%             %   
%             %     H2 = plot(0.25,t0oD25,'.k','markersize',18);
%             % 
%             % % Formatting    
%             %     HL = legend([H0,H1,H2],'Initial thickness','Modified thickness','ABS requirement');
%             %     
%             %     set(HL,'FontName','Times','FontSize',16,'Location','SouthWest')
%             %     
%             %     set(gca,'XTick',[0.2:0.1:1])
%             %     set(gca,'YTick',[0:0.005:0.03])
%             %     
%             %     set(gca,'XTickLabel',{'0.2','','0.4','','0.6','','0.8','','1'})
%             %     set(gca,'YTickLabel',{'0','','0.01','','0.02','','0.03'})
%             %     
%             %     axis([0.2 1 0 0.03])
%             % 
%             %     set(gca,'FontName','Times','FontSize',18)
%             %     xlabel('r / R','FontName','Times','FontSize',18)
%             %     ylabel('t0 / D','FontName','Times','FontSize',18)
%             %     
%             %     % saveas(gcf,'fig_ABS_t0oD','epsc')
%             %     
%             % % -------------------------------------------------------------------------
%             % % -------------------------------------------------------------------------
% 
% 
%             %%
%             % % -------------------------------------------------------------------------
%             % % FIGURE:  Plot Case2 contours of (-CPmax) for each blade section  
%             % % ------------------------------------------------------------------------- 
% %             clr,
% %             load DDG51design1_last_iteration
%                     % -----------------------------------------------------------------
%                     % METHOD 8: Epps/Viquez  -- 3 CASES METHOD         -- optimizes CoD and t0oc
%                     % -----------------------------------------------------------------                       
%                     
%                 % -------------------------------------------------------------------------
%                 % Parameters:
%                 NACAa       = 0.8;               % 'NACA a=0.8' meanline
%                 RHOLE       = 0.640279919559199; % == 0.00639/(2*0.04995)^2 == leading edge radius / chord, for t0oc==1 (RLE = RHOLE*t0oc^2)
%                 CAVmargin   = 1.0;               % allowable SIGMAs = CAVmargin*SIGMAs
%                 CoDmax      = 0.65;              % maximum allowable c/D
%                 CoDmax      = Inf;              % maximum allowable c/D
%                 % -------------------------------------------------------------------------
% 
%                 SIGMA     = SIGMAh./VSTARh.^2;      % local cavitation number
% 
%                 % -------------------------------------------------------------------------
%                 % CASE 2: Find optimum using Newton solver
%                 g1   = 4/pi;
%                 g2   = pi*G'./((1+NACAa)*VSTAR);
%                 g3   =              (Gh'./VSTARh - G'./VSTAR);
%                 g3sq = (2/RHOLE) .* (Gh'./VSTARh - G'./VSTAR).^2;  % not exactly g3^2
% 
%                 for i = 1:Mp
%                     [CoD_2(i),t0oc(i)] = Chord_ctSolver(CAVmargin*SIGMA(i), g1, g2(i), g3sq(i) , CoD(i), t0oc(i));
%                 end
% 
% 
%                 % -----------------------------------------------------------------
%                 % CASE 1:
%                 a1   = 4/pi;
%                 a2   = pi*G'./((1+NACAa)*VSTAR);
%                 a3   = (Gh'./VSTARh - G'./VSTAR);
% 
%                 % Solve for CoD required by CASE 1 for given t0oc calculated using CASE 2 or 3.
%                 %          (1+CAVmargin*SIGMA-(1+a1*t0oc).^2)
%                 %         -(2*(1+a1*t0oc).*(a2+a3))
%                 %         -(a2.^2+2*a2.*a3)
%                 %
%                 CoD_1 = ( (2*(1+a1*t0oc).*(a2+a3)) + sqrt( (2*(1+a1*t0oc).*(a2+a3)).^2 + 4*(1+CAVmargin*SIGMA-(1+a1*t0oc).^2).*(a2.^2+2*a2.*a3) ) )./( 2*(1+CAVmargin*SIGMA-(1+a1*t0oc).^2) );
%                 % -----------------------------------------------------------------
% 
%                     
%             tt = [0.001:0.001:0.2];   % t0/c
%             cc = [0.005:0.005:1];     %  c/D
%             [t,c] = ndgrid(tt,cc);
%             
%             for i = 1:Mp
%             
%                 CP2 = (1 + g1*t + g2(i)./c).^2 + g3sq(i)./(c.^2.*t.^2) - 1;
%             
%                 fig;
%                     contourf(c,t,CP2,[-0.1:0.1:2])
%                 H1 = contour(c,t,CP2,SIGMA(i)*[1 1],'color','r','linewidth',2);   % contour of CP2 == SIGMA
%             
%                 % plot([0,1],t0oc(i)*[1 1],'r--','linewidth',2)
%             
%                 % plot([0.005:0.005:1],t0oD(i)./[0.005:0.005:1],'g--','linewidth',2)
%             
%                 % Case 1: Contour of CP1 == SIGMA, CoD versus tt
%                 CoD_case1 = ( (2*(1+a1*tt).*(a2(i)+a3(i))) + sqrt( (2*(1+a1*tt).*(a2(i)+a3(i))).^2 + 4*(1+CAVmargin*SIGMA(i)-(1+a1*tt).^2).*(a2(i).^2+2*a2(i).*a3(i))           ) )./( 2*(1+CAVmargin*SIGMA(i)-(1+a1*tt).^2) );
% 
%                 % Case 2: Contour of CP2 == SIGMA, CoD versus tt
%                 CoD_case2 = ( (2*(1+a1*tt).*a2(i))         + sqrt( (2*(1+a1*tt).*a2(i)).^2         + 4*(1+CAVmargin*SIGMA(i)-(1+a1*tt).^2).*(a2(i).^2+2*a3(i)^2./(RHOLE*tt.^2)) ) )./( 2*(1+CAVmargin*SIGMA(i)-(1+a1*tt).^2) );
% 
%                 CoD_case1( CoD_case1 < 0 ) = NaN;
%                 CoD_case2( CoD_case2 < 0 ) = NaN;
%                 
%                 plot( CoD_case1 , tt, 'r--')
%               % plot( CoD_case2 , tt, 'r--')
%                 
%                 % Case 2: Optimum chord length for a given thickness ratio
%                 plot((-g1*g2(i)+sqrt( (g1*g2(i))^2 + 4*(g1+g1^2*tt)*g3sq(i)./tt.^3 ))./(2*g1+2*g1^2*tt),tt,'y--','linewidth',2)
%             
%                 % Cross-over point between case 1 and case 2
%                 CoD_cross   = (g3(i)./(RHOLE*tt.^2) - g2(i)) ./ (1+g1*tt);
%                 plot( CoD_cross, tt , 'm--','linewidth',2)                
%             
%                 % Selected chord and thickness ratio
%                 plot(CoD_2(i),t0oc(i),'.','markersize',22,'color',[0 0.9 0])
%             
%                 % Formatting                                
%                 set(gca,'XTick',[0:0.2:1])
%                 set(gca,'YTick',[0:0.05:0.2])
%             
%                 set(gca,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'})
%                 set(gca,'YTickLabel',{'0','0.05','0.10','0.15','0.20'})
%             
%                 caxis([0 2]);  % h = colorbar; 
%             
%                 axis([0 1 0 0.2])
%                 box on
%             
%                 set(gca,'FontName','Times','FontSize',18)
%                 xlabel('c / D','FontName','Times','FontSize',18)
%                 ylabel('t0 / c','FontName','Times','FontSize',18)
%             
%             
%             
% %                 if i < 10
% %                     saveas(gcf,['fig_CPmax_case2_0',num2str(i)],'epsc')
% %                 else
% %                     saveas(gcf,['fig_CPmax_case2_',num2str(i)],'epsc')
% %                 end
%             end
%             % -------------------------------------------------------------------------
%             % -------------------------------------------------------------------------
%             %%                            
% 
% % 
% %             %%
% %             % % -------------------------------------------------------------------------
% %             % % FIGURE:  Plot Case1 contours of (-CPmax) for each blade section  --- first run Case 2 and Case 1 code blocks above
% %             % % ------------------------------------------------------------------------- 
% % %             clr,
% % %             load DDG51design1_last_iteration
% %                     % -----------------------------------------------------------------
% %                     % METHOD 8: Epps/Viquez  -- 3 CASES METHOD         -- optimizes CoD and t0oc
% %                     % -----------------------------------------------------------------    
% % 
% %                 % -------------------------------------------------------------------------
% %                 % Parameters:
% %                 NACAa       = 0.8;               % 'NACA a=0.8' meanline
% %                 RHOLE       = 0.640279919559199; % == 0.00639/(2*0.04995)^2 == leading edge radius / chord, for t0oc==1 (RLE = RHOLE*t0oc^2)
% %                 CAVmargin   = 1.0;               % allowable SIGMAs = CAVmargin*SIGMAs
% %                 CoDmax      = 0.65;              % maximum allowable c/D
% %                 CoDmax      = Inf;              % maximum allowable c/D
% %                 % -------------------------------------------------------------------------
% % 
% %                 SIGMA     = SIGMAh./VSTARh.^2;      % local cavitation number
% % 
% %                 % -------------------------------------------------------------------------
% %                 % CASE 2: Find optimum using Newton solver
% %                 g1   = 4/pi;
% %                 g2   = pi*G'./((1+NACAa)*VSTAR);
% %                 g3   =              (Gh'./VSTARh - G'./VSTAR);
% %                 g3sq = (2/RHOLE) .* (Gh'./VSTARh - G'./VSTAR).^2;  % not exactly g3^2
% % 
% %                 for i = 1:Mp
% %                     [CoD_2(i),t0oc(i)] = Chord_ctSolver(CAVmargin*SIGMA(i), g1, g2(i), g3sq(i) , CoD(i), t0oc(i));
% %                 end
% % 
% % 
% %                 % -----------------------------------------------------------------
% %                 % CASE 1:
% %                 a1   = 4/pi;
% %                 a2   = pi*G'./((1+NACAa)*VSTAR);
% %                 a3   = (Gh'./VSTARh - G'./VSTAR);
% % 
% %                 % Solve for CoD required by CASE 1 for given t0oc calculated using CASE 2 or 3.
% %                 %          (1+CAVmargin*SIGMA-(1+a1*t0oc).^2)
% %                 %         -(2*(1+a1*t0oc).*(a2+a3))
% %                 %         -(a2.^2+2*a2.*a3)
% %                 %
% %                 CoD_1 = ( (2*(1+a1*t0oc).*(a2+a3)) + sqrt( (2*(1+a1*t0oc).*(a2+a3)).^2 + 4*(1+CAVmargin*SIGMA-(1+a1*t0oc).^2).*(a2.^2+2*a2.*a3) ) )./( 2*(1+CAVmargin*SIGMA-(1+a1*t0oc).^2) );
% %                 % -----------------------------------------------------------------                    
% %                     
% %             tt = [0.001:0.001:0.2];
% %             cc = [0.005:0.005:1];
% %             [t,c] = ndgrid(tt,cc);
% %             
% %             for i = 1:Mp
% %             
% %                 CP1 = (1 + a1*t + a2(i)./c).^2  +  2*(1 + a1*t + a2(i)./c).*(a3(i)./c)   - 1;  % CASE 1
% %                    
% %                 fig;
% %                     contourf(c,t,CP1,[-0.1:0.1:2])
% %                 H1 = contour(c,t,CP1,SIGMA(i)*[1 1],'color','r','linewidth',2); 
% %             
% %                 % Case 1: Contour of CP1 == SIGMA, CoD versus tt
% %                 CoD_case1 = ( (2*(1+a1*tt).*(a2(i)+a3(i))) + sqrt( (2*(1+a1*tt).*(a2(i)+a3(i))).^2 + 4*(1+CAVmargin*SIGMA(i)-(1+a1*tt).^2).*(a2(i).^2+2*a2(i).*a3(i))           ) )./( 2*(1+CAVmargin*SIGMA(i)-(1+a1*tt).^2) );
% % 
% %                 % Case 2: Contour of CP2 == SIGMA, CoD versus tt
% %                 CoD_case2 = ( (2*(1+a1*tt).*a2(i))         + sqrt( (2*(1+a1*tt).*a2(i)).^2         + 4*(1+CAVmargin*SIGMA(i)-(1+a1*tt).^2).*(a2(i).^2+2*a3(i)^2./(RHOLE*tt.^2)) ) )./( 2*(1+CAVmargin*SIGMA(i)-(1+a1*tt).^2) );
% % 
% %                 CoD_case1( CoD_case1 < 0 ) = NaN;
% %                 CoD_case2( CoD_case2 < 0 ) = NaN;
% %                 
% % %                 plot( CoD_case1 , tt, 'b--')
% %               % plot( CoD_case2 , tt, 'r--')
% %                 
% % %                 % Case 2: Optimum chord length for a given thickness ratio
% % %                 plot((-g1*g2(i)+sqrt( (g1*g2(i))^2 + 4*(g1+g1^2*tt)*g3sq(i)./tt.^3 ))./(2*g1+2*g1^2*tt),tt,'y--','linewidth',2)
% %             
% %                 % Cross-over point between case 1 and case 2
% %                 CoD_cross   = (g3(i)./(RHOLE*tt.^2) - g2(i)) ./ (1+g1*tt);
% %                 plot( CoD_cross, tt , 'm--','linewidth',2)                
% %             
% %                 
% % %             
% % %                 % Selected chord and thickness ratio
% % %                 plot(CoD_1(i),t0oc(i),'xr','markersize',10,'linewidth',2)
% % %                 
% % %                 % Selected chord and thickness ratio
% % %                 plot(CoD_2(i),t0oc(i),'.','markersize',22,'color',[0 0.9 0])
% %             
% %                 % Formatting                                
% %                 set(gca,'XTick',[0:0.2:1])
% %                 set(gca,'YTick',[0:0.05:0.2])
% %             
% %                 set(gca,'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'})
% %                 set(gca,'YTickLabel',{'0','0.05','0.10','0.15','0.20'})
% %             
% %                 caxis([0 2]);  % h = colorbar; 
% %             
% %                 axis([0 1 0 0.2])
% %                 box on
% %             
% %                 set(gca,'FontName','Times','FontSize',18)
% %                 xlabel('c / D','FontName','Times','FontSize',18)
% %                 ylabel('t0 / c','FontName','Times','FontSize',18)
% %             
% %             
% %             
% % %                 if i < 10
% % %                     saveas(gcf,['fig_CPmax_case1_0',num2str(i)],'epsc')
% % %                 else
% % %                     saveas(gcf,['fig_CPmax_case1_',num2str(i)],'epsc')
% % %                 end
% %             end
% %             % -------------------------------------------------------------------------
% %             % -------------------------------------------------------------------------
% % 
% % 
% %                     %         % -----------------------------------------------------------------
% %                     %         % CASE 3:                      
% %                     %         g1   = 4/pi;
% %                     %         g2   = pi*G'./((1+NACAa)*VSTAR);
% %                     %         g3sq = (2/RHOLE) .* (Gh'./VSTARh - G'./VSTAR).^2;
% %                     % 
% %                     % %         % Find optimum using Newton solver
% %                     % %         for i = 1:Mp
% %                     % %             [CoD6h(i),t0oc(i)] = ctSolver(SIGMA(i), g1, g2(i), g3sq(i) , CoD6h(i), t0oc(i));
% %                     % %         end
% % 
% % 
% %                     %         % Plot contours of (-CPmax) for each blade section
% %                     %         close all
% %                     %         position1 = [746   664   560   420];
% %                     %         position2 = [1308  664   560   420];
% %                     %         position3 = [747   170   560   420];      
% %                     %         position4 = [1308  170   560   420];
% %                     %         
% %                     % 
% %                     %         tt = [0.001:0.001:0.2];
% %                     %         cc = [0.005:0.005:1];
% %                     %         [t,c] = ndgrid(tt,cc);
% %                     % 
% %                     %         for i = 1:Mp
% %                     % 
% %                     %             CP1 = (1 + k1*t + k2(i)./c).^2  +  2*(1 + k1*t + k2(i)./c).*(k3(i)./c)   - 1;  % CASE 1
% %                     %             
% %                     %             CP2 = (1 + g1*t + g2(i)./c).^2  +  g3sq(i)./(c.^2.*t.^2) - 1;  % CASE 2  (k3 > 0)
% %                     %             CP3 = (1 + g1*t - g2(i)./c).^2  +  g3sq(i)./(c.^2.*t.^2) - 1;  % CASE 3  (k3 < 0)
% %                     %             
% %                     %             if k3(i) >= 0
% %                     %                 v = max(CP1,CP2); 
% %                     %             else
% %                     %                 v = max(CP1,CP3); 
% %                     %             end
% %                     %             
% %                     %             fig(1); clf, hold on,
% %                     %                 set(gcf,'Position',position1)
% %                     %                 xlabel('c/D'), ylabel('t0/c'), colorbar
% %                     % 
% %                     %                     contourf(c,t,CP1,[0:0.1:2])
% %                     %                 H1 = contour(c,t,CP1,SIGMA(i)*[1 1],'color','k','linewidth',2); 
% %                     % 
% %                     %                 plot([0,1],t0oc0(i)*[1 1],'r--','linewidth',2)
% %                     % 
% %                     %                 plot([0.005:0.005:1],t0oD(i)./[0.005:0.005:1],'g--','linewidth',2)
% %                     % 
% %                     %     %             plot(CoD6h(i),t0oc(i),'m.','markersize',18)
% %                     % 
% %                     %     %             plot((-g1*g2(i)+sqrt( (g1*g2(i))^2 + 4*(g1+g1^2*tt)*g3sq(i)./tt.^3 ))./(2*g1+2*g1^2*tt),tt,'y--','linewidth',2)
% %                     % 
% %                     %             fig(2); clf,hold on,
% %                     %                 set(gcf,'Position',position2)
% %                     %                 xlabel('c/D'), ylabel('t0/c'), colorbar
% %                     % 
% %                     %                     contourf(c,t,CP2,[0:0.1:2])
% %                     %                 H1 = contour(c,t,CP2,SIGMA(i)*[1 1],'color','k','linewidth',2); 
% %                     % 
% %                     %                 
% %                     %             fig(3); clf,hold on,
% %                     %                 set(gcf,'Position',position3)
% %                     %                 xlabel('c/D'), ylabel('t0/c'), colorbar
% %                     % 
% %                     %                     contourf(c,t,CP3,[0:0.1:2])
% %                     %                 H1 = contour(c,t,CP3,SIGMA(i)*[1 1],'color','k','linewidth',2); 
% %                     %                 
% %                     %                 
% %                     %             fig(4); clf,hold on,
% %                     %                 set(gcf,'Position',position4)
% %                     %                 xlabel('c/D'), ylabel('t0/c'), colorbar
% %                     % 
% %                     %                     contourf(c,t,v,[0:0.1:2])
% %                     %                 H1 = contour(c,t,v,SIGMA(i)*[1 1],'color','k','linewidth',2); 
% %                     %                 
% %                     %             pause,
% %                     %         end
% % 
% %                             % saveas(gcf,'110223 CPmax contours','epsc')
% % 
% %                     % -------------------------------------------------------------------------
% %                     % -------------------------------------------------------------------------
% % %         end
% % %         -----------------------------------------------------------------
%                     
%       
        