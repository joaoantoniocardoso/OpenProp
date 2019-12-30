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
% ===================================================== Align_wake Function
% Modified: 11/2/2011, Brenden Epps
% -------------------------------------------------------------------------
%
% This function aligns the wake to the given circulation distribution by
% iteratively computing:
% 	UAHIF,UTHIF   = the horseshoe influence functions
%   UASTAR,UTSTAR = the induced velocities
%   TANBIC        = the inflow velocity pitch angle
%
% This implementation does not iteratively update the duct circulation.
% Note:    UADUCT = UADIF*Gd == influence of duct on propeller
% -------------------------------------------------------------------------

function [UASTAR,UTSTAR,TANBIC,UAHIF,UTHIF,converged_flag] =    ...
         Align_wake(G,  VAC,VTC, TANBIC,RC,RV,L,Mp,Z, ...
                    Hub_flag,Rhub_oR,  Duct_flag,Rduct_oR,UADUCT,   Hvel,Hbeta)
    
                
if Duct_flag == 1
    disp('ERROR: Duct model is not implemented in Align_wake.m')
else
    Gd = 0;
    UADUCT = 0;
end

ITER  = 50;   % max number of iterations
relax = 0.9;  % relaxation parameter

% ------------------------------------------------ Implement RepairSpline.m
% To smooth data X(RC):  X_smooth = X*Bsmooth;
Bsmooth = RepairSplineMatrix(RC);

% Smooth the inflow angle for numerical stability:
TANBICsmooth = TANBIC * Bsmooth;    
% -------------------------------------------------------------------------   


% ----------------- Iterate to ALIGN WAKE to circulation distribution G
W_iter = 1;                     % iteration in the wake alignment loop
W_res  = 1;                     % residual BetaI between iterations
TANBIC_last = TANBIC;           % the last value of TANBIC

while W_iter <= ITER && any(W_res > 0.01)           % (WHILE LOOP WA1)
%         disp(['--- Wake alignment iteration: ',num2str(W_iter)])       

    % ---------------------------------------------------------------------
    [UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBICsmooth,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);
    
%     [UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBIC      ,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);
    
    
    % ---------------------------------------------------------------------
    % Update induced velocities on propeller:
    UASTAR = (UAHIF*G)';  
    UTSTAR = (UTHIF*G)';  % influence of propeller on propeller  

    % ---------------------------------------------------------------------
    TANBIC = (VAC + UADUCT + UASTAR)./(L*RC + VTC + UTSTAR);

    
    TANBIC = relax*TANBIC + (1-relax)*TANBIC_last;
    

    % Smooth the inflow angle for numerical stability:
    TANBICsmooth = TANBIC * Bsmooth;    
    
    % -------------------------------------- Prepare for the next iteration
    W_iter      = W_iter + 1;                  
    W_res       = abs((TANBIC - TANBIC_last)./TANBIC);
    TANBIC_last = TANBIC;                 

    if any(isnan(TANBIC))  || ~isreal(TANBIC) || any(abs(TANBIC) > 10)
         UASTAR = 0*RC;
         UTSTAR = 0*RC;
         TANBIC = 0*RC;
         UAHIF  = zeros(Mp,Mp);
         UTHIF  = zeros(Mp,Mp);

         disp('<WARNING>')
         disp('<WARNING> TANBIC == too large, NaN, or imaginary... crash avoided...')
         disp('<WARNING>')
         W_iter         = 999;
         converged_flag = 0;
    end     
    
    if any( (VAC + UASTAR) < 0 )
         UASTAR = 0*RC;
         UTSTAR = 0*RC;
         TANBIC = 0*RC;
         UAHIF  = zeros(Mp,Mp);
         UTHIF  = zeros(Mp,Mp);
        

         disp('<WARNING>')
         disp('<WARNING> Flow reversal...(VAC + UASTAR) < 0... crash avoided...')
         disp('<WARNING>')
         W_iter         = 999;
         converged_flag = 0;
    end     
    
%     if W_iter > ITER
%         disp('WARNING: While loop WA1 did NOT converge.'),
%     end
% 
%     disp(['The max W_res for iteration ',num2str(W_iter-1),' is: ',num2str(max(W_res))]),       
%         disp(' '),       

        if (Hvel ~=0)  & (Hbeta ~=0)
                % ------------------------------- Uncomment to plot UASTAR & UTSTAR
                figure(Hvel),
                ylabel('RED = UASTAR, BLUE = UTSTAR')
                    if     mod(W_iter-1,4) == 1
                        plot(RC,UASTAR,'r.-'), 
                    else
                        plot(RC,UASTAR,'m.-'), 
                    end               
                
                    if     mod(W_iter-1,4) == 1
                        plot(RC,UTSTAR,'b.-'), 
                    else
                        plot(RC,UTSTAR,'c.-'), 
                    end               

                % ---------------------------------------- Uncomment to plot TANBIC
                figure(Hbeta),
                ylabel('TANBIC')
                    if     mod(W_iter-1,4) == 1
                        plot(RC,TANBIC,'r.-'), hold on,
                    elseif mod(W_iter-1,4) == 2
                        plot(RC,TANBIC,'m.-'), hold on,
                    elseif mod(W_iter-1,4) == 3
                        plot(RC,TANBIC,'b.-'), hold on,
                    else
                        plot(RC,TANBIC,'k.-'), hold on,
                    end
        end
        % --------------------------- Uncomment to pause between iterations     
        pause(0.1),

end                                              % (END WHILE LOOP WA1)
    
    if W_iter > ITER
        disp('WARNING: Align_wake.m did NOT converge.'),
        converged_flag = 0;
    else
        converged_flag = 1;
    end

%     disp(['[Align_wake] The max W_res for iteration ',num2str(W_iter-1),' is: ',num2str(max(W_res))]),
% ================================================= END Align_wake Function
% =========================================================================