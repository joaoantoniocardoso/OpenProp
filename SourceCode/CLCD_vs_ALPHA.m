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
% ================================================== CLCD_vs_ALPHA Function
% alpha      = [rad], current   angle of attack
% alpha0     = [rad], reference angle of attack
% ALPHA      = [rad], net angle of attack == alpha - alpha0
% ALPHAstall = [rad], net angle of attack at stall

% CL0        = [ ],     lift coefficient at alpha0, (1 x Mp)
% CD0        = [ ],     drag coefficient at alpha0, (1 x Mp)
% dCLdALPHA  = [1/rad], lift curve slope,           (1 x  1)
% -------------------------------------------------------------------------

function [CL,CD] = CLCD_vs_ALPHA(ALPHA,ALPHAstall,CL0,CD0,dCLdALPHA,Propeller_flag)

if     nargin == 1
    ALPHAstall = 8*pi/180;
    disp('Stall  ALPHA      has been set to 8 deg.')

    CL0 = 0;
    disp('CL for ALPHA == 0 has been set to zero.')    

    CD0 = 0;
    disp('CD for ALPHA == 0 has been set to 0.')

    dCLdALPHA = 2*pi;
    disp('dCLdALPHA         has been set to 2*pi.')
    
    Propeller_flag = 1;
    disp('Propeller_flag    has been set to 1.')
   
elseif nargin == 2
    CL0 = 0;
    disp('CL for ALPHA == 0 has been set to zero.')    

    CD0 = 0;
    disp('CD for ALPHA == 0 has been set to 0.')

    dCLdALPHA = 2*pi;
    disp('dCLdALPHA         has been set to 2*pi.')
    
    Propeller_flag = 1;
    disp('Propeller_flag    has been set to 1.')
       
elseif nargin == 3,
    CD0 = 0;
    disp('CD for ALPHA == 0 has been set to 0.')

    dCLdALPHA = 2*pi;
    disp('dCLdALPHA         has been set to 2*pi.')
    
    Propeller_flag = 1;
    disp('Propeller_flag    has been set to 1.')
    
elseif nargin == 4,
    dCLdALPHA = 2*pi;
    disp('dCLdALPHA         has been set to 2*pi.')
    
    Propeller_flag = 1;
    disp('Propeller_flag    has been set to 1.')
    
elseif nargin == 5,
    Propeller_flag = 1;
    disp('Propeller_flag    has been set to 1.')
       
end


%%
% % ----------------------------------------------------- TEST & DEBUG CODE
% clear, close all, clc, 
% ALPHA = [-90:0.1:90]*pi/180;
% ALPHAstall = 8*pi/180;
% CL0        = 0.4;
% CD0        = 0.01;
% dCLdALPHA  = 1.0 * 2*pi;
% CDoCL      = CD0/CL0;
% Propeller_flag = 0;
% % ------------------------------------------------- END TEST & DEBUG CODE

% drag / lift ratio
if CL0 == 0
    CDoCL = CD0 ./ (2*pi*ALPHAstall);
else
    CDoCL = CD0 ./ abs(CL0);
end

if dCLdALPHA == 0  % then section is not a lifting surface
    CL = CL0; 
    CD = CD0;
    
else
    
    B  = 20;                               % stall sharpness parameter

    CL =  CL0 ...
        + dCLdALPHA .*   ALPHA ...
        - dCLdALPHA .* ( ALPHA-ALPHAstall).*((1/pi)*atan(B*( ALPHA-ALPHAstall))+0.5) ...
        + dCLdALPHA .* (-ALPHA-ALPHAstall).*((1/pi)*atan(B*(-ALPHA-ALPHAstall))+0.5);    
    
    if Propeller_flag == 1

        % ----------------------------------------------- (Propeller) constant drag treatment
        % % ---- CD(ALPHA) ~= CD0  near ALPHA == 0
        A  = (2-CD0)./(pi/2-ALPHAstall);  % drag coefficient post-stall slope        

        CD =  CD0 ... 
            +         A.*( ALPHA-ALPHAstall).*((1/pi)*atan(B*( ALPHA-ALPHAstall))+0.5) ... 
            -         A.*(      -ALPHAstall).*((1/pi)*atan(B*(      -ALPHAstall))+0.5) ...  
            +         A.*(-ALPHA-ALPHAstall).*((1/pi)*atan(B*(-ALPHA-ALPHAstall))+0.5) ...
            -         A.*(      -ALPHAstall).*((1/pi)*atan(B*(      -ALPHAstall))+0.5);

    else

        % ------------------------------------------------ Turbine drag treatment 
        % ---- CD(ALPHA) ~= abs(CL0).*(CDoCL + dCLdALPHA*ALPHA)  near ALPHA == 0
        A  = (2-CDoCL.*(CL0+dCLdALPHA.*ALPHAstall))./(pi/2-ALPHAstall);  % drag coefficient post-stall slope        

        CD =  CDoCL.*abs(CL0) ...
            + CDoCL.*dCLdALPHA.*( ALPHA           ).*((1/pi)*atan(B*( ALPHA           ))+0.5) ... 
            - CDoCL.*dCLdALPHA.*( ALPHA-ALPHAstall).*((1/pi)*atan(B*( ALPHA-ALPHAstall))+0.5) ... 
            + CDoCL.*dCLdALPHA.*(      -ALPHAstall).*((1/pi)*atan(B*(      -ALPHAstall))+0.5) ... 
            +                A.*( ALPHA-ALPHAstall).*((1/pi)*atan(B*( ALPHA-ALPHAstall))+0.5) ... 
            -                A.*(      -ALPHAstall).*((1/pi)*atan(B*(      -ALPHAstall))+0.5) ...  
            + CDoCL.*dCLdALPHA.*(-ALPHA           ).*((1/pi)*atan(B*(-ALPHA           ))+0.5) ...
            - CDoCL.*dCLdALPHA.*(-ALPHA-ALPHAstall).*((1/pi)*atan(B*(-ALPHA-ALPHAstall))+0.5) ... 
            + CDoCL.*dCLdALPHA.*(      -ALPHAstall).*((1/pi)*atan(B*(      -ALPHAstall))+0.5) ...     
            +                A.*(-ALPHA-ALPHAstall).*((1/pi)*atan(B*(-ALPHA-ALPHAstall))+0.5) ...
            -                A.*(      -ALPHAstall).*((1/pi)*atan(B*(      -ALPHAstall))+0.5);
    end
end


% % % % ----------------------------------------------------- TEST & DEBUG CODE
% % %
% % %         % ----------------------------------------------- (Propeller) constant drag treatment
% % %         % % ---- CD(ALPHA) ~= CD0  near ALPHA == 0
% % %         A  = (2-CD0)/(pi/2-ALPHAstall);  % drag coefficient post-stall slope        
% % % 
% % %         CD =  CD0 ... 
% % %             +         A.*( ALPHA-ALPHAstall).*((1/pi)*atan(B*( ALPHA-ALPHAstall))+0.5) ... 
% % %             -         A.*(      -ALPHAstall).*((1/pi)*atan(B*(      -ALPHAstall))+0.5) ...  
% % %             +         A.*(-ALPHA-ALPHAstall).*((1/pi)*atan(B*(-ALPHA-ALPHAstall))+0.5) ...
% % %             -         A.*(      -ALPHAstall).*((1/pi)*atan(B*(      -ALPHAstall))+0.5);
% % % 
% % % CDP=CD;
% % % 
% % %         % ------------------------------------------------ Turbine drag treatment 
% % %         % ---- CD(ALPHA) ~= abs(CL0).*(CDoCL + dCLdALPHA*ALPHA)  near ALPHA == 0
% % %         A  = (2-CDoCL*(CL0+dCLdALPHA*ALPHAstall))/(pi/2-ALPHAstall);  % drag coefficient post-stall slope        
% % % 
% % %         CD =  CDoCL*abs(CL0) ...
% % %             + CDoCL*dCLdALPHA*( ALPHA           ).*((1/pi)*atan(B*( ALPHA           ))+0.5) ... 
% % %             - CDoCL*dCLdALPHA*( ALPHA-ALPHAstall).*((1/pi)*atan(B*( ALPHA-ALPHAstall))+0.5) ... 
% % %             + CDoCL*dCLdALPHA*(      -ALPHAstall).*((1/pi)*atan(B*(      -ALPHAstall))+0.5) ... 
% % %             +              A.*( ALPHA-ALPHAstall).*((1/pi)*atan(B*( ALPHA-ALPHAstall))+0.5) ... 
% % %             -              A.*(      -ALPHAstall).*((1/pi)*atan(B*(      -ALPHAstall))+0.5) ...  
% % %             + CDoCL*dCLdALPHA*(-ALPHA           ).*((1/pi)*atan(B*(-ALPHA           ))+0.5) ...
% % %             - CDoCL*dCLdALPHA*(-ALPHA-ALPHAstall).*((1/pi)*atan(B*(-ALPHA-ALPHAstall))+0.5) ... 
% % %             + CDoCL*dCLdALPHA*(      -ALPHAstall).*((1/pi)*atan(B*(      -ALPHAstall))+0.5) ...     
% % %             +              A.*(-ALPHA-ALPHAstall).*((1/pi)*atan(B*(-ALPHA-ALPHAstall))+0.5) ...
% % %             -              A.*(      -ALPHAstall).*((1/pi)*atan(B*(      -ALPHAstall))+0.5);
% % % 
% % %         
% % % CDT=CD;    
% % % 
% % %         % % % ---- CD(ALPHA) ~= abs(CL0).*CDoCL  near ALPHA == 0
% % %         % % A  = (2-CDoCL*CL0)/(pi/2-ALPHAstall);  % drag coefficient post-stall slope     
% % %         % % 
% % %         % % CD = abs(CL0).*CDoCL ...
% % %         % %     + A.*( ALPHA-ALPHAstall).*((1/pi)*atan(B*( ALPHA-ALPHAstall))+0.5) ...
% % %         % %     + A.*(-ALPHA-ALPHAstall).*((1/pi)*atan(B*(-ALPHA-ALPHAstall))+0.5) ...
% % %         % %     - A.*(      -ALPHAstall).*((1/pi)*atan(B*(      -ALPHAstall))+0.5) ...
% % %         % %     - A.*(      -ALPHAstall).*((1/pi)*atan(B*(      -ALPHAstall))+0.5);
% % %         % 
% % %         % 
% % %         % % % ---- CDoCL == constant
% % %         % % CD =  CDoCL*abs(CL);
% % % 
% % ----------------------------------------------------- TEST & DEBUG CODE
% % CL(find(ALPHA ==0))
% % CD(find(ALPHA ==0))
% % CD(1)
% % CD(find(ALPHA ==ALPHAstall))
% % CD(find(ALPHA ==-ALPHAstall))
% 
% % figure,
% %     plot(ALPHA,CL0+2*pi*ALPHA,'r--'), hold on,
% %     plot(ALPHA,CL,'k')
% %     plot(ALPHA,CD,'b')
% %     
% %     plot([ALPHAstall,ALPHAstall],[-2,2],'k:',[-ALPHAstall,-ALPHAstall],[-2,2],'k:')
% %     plot([-pi/2 pi/2],[0 0],'k:',[0 0],[-2,2],'k:')
% %     axis([-pi/2 pi/2 -2 2])
% %     box on
% %     xlabel('angle of attack','FontSize',16)
% %     ylabel('lift coefficient, CL','FontSize',16)
% %     set(gca,'FontSize',14)
% figure,
%     plot(ALPHA*180/pi,CL0+2*pi*ALPHA,'r--'), hold on,
% %     plot(ALPHA*180/pi,CD0 + CD0 * ALPHA.^2/ALPHAstall^2,'r--')
% %     plot(ALPHA*180/pi,CD0 + CD0 * abs(ALPHA)/ALPHAstall,'g--')
%     plot(ALPHA*180/pi,CL,'k')
%     plot(ALPHA*180/pi,CD,'b')
%     
% %     % Blade-normal force
% %     plot(ALPHA*180/pi,CL.*cos(ALPHA)+CD.*sin(ALPHA),'r')
%    
%     plot([ALPHAstall,ALPHAstall]*180/pi,[-2,2],'k:',[-ALPHAstall,-ALPHAstall]*180/pi,[-2,2],'k:')
%     plot([-90 90],[0 0],'k:',[0 0],[-2,2],'k:')
%     axis([-90 90 -2 2])
%     box on, grid on,
%     xlabel('ALPHA (net angle of attack) [deg]','FontSize',16)
%     ylabel('lift coefficient, CL','FontSize',16)
%     set(gca,'FontSize',14)
%     
% % figure, hold on, grid on,
% %     plot(CL,CDP,'b')
% %     plot(CL,CDT,'r')
% %     plot(CL0+[-1 1],CD0*[1 1],'k--')
% %     axis([CL0-0.7 CL0+0.7 CD0-0.004 CD0+0.014])
% % ------------------------------------------------- END TEST & DEBUG CODE

% ================================================================
% ================================================================