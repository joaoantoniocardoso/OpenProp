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
% using a Newton solver to determine:
% 	UAHIF,UTHIF   = the horseshoe influence functions
%   UASTAR,UTSTAR = the induced velocities
%   TANBIC        = the inflow velocity pitch angle
%
% This implementation does not iteratively update the duct circulation.
% Note:    UADUCT = UADIF*Gd == influence of duct on propeller
% -------------------------------------------------------------------------

function [UASTAR,UTSTAR,TANBIC,UAHIF,UTHIF,converged_flag] =    ...
         Align_wake_Newton(G,  VAC,VTC, TANBIC,RC,RV,L,Mp,Z, ...
                    Hub_flag,Rhub_oR,  Duct_flag,Rduct_oR,UADUCT,   Hvel,Hbeta)
    
                
if Duct_flag == 1
    disp('ERROR: Duct model is not implemented in Align_wake_Newton.m')
else
    Gd = 0;
    UADUCT = 0*G;
end

% ------------------------------------------------ Implement RepairSpline.m
% To smooth data X(RC):  X_smooth = X*Bsmooth;
Bsmooth = RepairSplineMatrix(RC);
% -------------------------------------------------------------------------   


% =====================================================================
% ==   FIND STATE OF SYSTEM USING NEWTON SOLVER:
% ==                   iterate to solve residual equations for unknowns
N_iter   = 1;                            
ERROR    = 1;
ERRORtol = 0.00001;  % 0.005
relax    = 0.9;
ITER     = 50;   % max number of iterations


   % -------------- Initialize linear system of equations matrices
    R  = zeros(3*Mp,1   );            % R  = vector of residuals
    J  =   eye(3*Mp,3*Mp);            % J  = matrix of derivatives
    DX = zeros(3*Mp,1   );            % DX = vector of change in unknowns 
 %  X  = zeros(3*Mp,1   );            % X  = vector of unknowns              

    % X(     1 :   Mp) =  TANBIC(:);   % X1
    % X(1*Mp+1 : 2*Mp) =  UASTAR(:);   % X2
    % X(2*Mp+1 : 3*Mp) =  UTSTAR(:);   % X3
    % -------------------------------------------------------------------------
    
    
    % --------------------------------------------------------------------- INITIALIZE UASTAR and UTSTAR
    % Smooth the inflow angle for numerical stability:
    TANBICsmooth = TANBIC * Bsmooth;
    
    % Update the horseshoe influence functions:
    
    [UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBICsmooth,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);
    
%     [UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBIC      ,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);


    % Evaluate velocities for R2 & R3
    UASTAR = (UAHIF*G)';  
    UTSTAR = (UTHIF*G)';     
    % ---------------------------------------------------------------------
    
    
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
while N_iter <= ITER && any(abs(ERROR) > ERRORtol)          % (WHILE LOOP N1)
    % disp(['------- Newton iteration: ',num2str(N_iter)])  % status message
    % disp(' '),
    
    % ---------------------------------- Store last state of the system
     TANBIClast =  TANBIC;
     UASTARlast =  UASTAR;
     UTSTARlast =  UTSTAR;  
    
    % ---------------------------------------------------------------------
    % Smooth the inflow angle for numerical stability:
    TANBICsmooth = TANBIC * Bsmooth;
    
    % Update the horseshoe influence functions:
    
    [UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBICsmooth,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);
    
%     [UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBIC      ,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);


    % Evaluate velocities for R2 & R3
    UASTARtemp = (UAHIF*G)';  
    UTSTARtemp = (UTHIF*G)';     
    % ---------------------------------------------------------------------
    

    % ------------- Solve Newton problem to update entire blade at once
    for m = 1:Mp  
        % ------------------------------------------------------ Evaluate residuals
        R(       m) =  TANBIC(m) -  (VAC(m)+UADUCT(m)+UASTAR(m))/(L*RC(m)+VTC(m)+UTSTAR(m)); % R1
        
        R(1*Mp + m) =  UASTAR(m) - UASTARtemp(m);      % R2

        R(2*Mp + m) =  UTSTAR(m) - UTSTARtemp(m);      % R3
        % -------------------------------------------------------------------------


        % --------------------------------------- Evaluate residual derivatives
        % disp('Evaluating residual derivatives matrix, J(m,n) = dR(m)/dX(n)...')
        % disp(' '),

        J(m ,1*Mp + m) =  -                           1/(L*RC(m)+VTC(m)+UTSTAR(m));
        J(m ,2*Mp + m) =   (VAC(m)+UADUCT(m)+UASTAR(m))/(L*RC(m)+VTC(m)+UTSTAR(m))^2;
        
%         for n = 1:Mp
%             J(m ,1*Mp + n) =  -                           1/(L*RC(m)+VTC(m)+UTSTAR(m));
%             J(m ,2*Mp + n) =   (VAC(m)+UADUCT(m)+UASTAR(m))/(L*RC(m)+VTC(m)+UTSTAR(m))^2;
%             
% %             J(1*Mp + m,2*Mp + n) = -UAHIF(m,n) / UTHIF(n,n);
% %             J(2*Mp + m,1*Mp + n) = -UTHIF(m,n) / UAHIF(n,n); 
%         end
        
%         J(1*Mp + m,2*Mp + m) = -UAHIF(m,m) / UTHIF(m,m);
%         J(2*Mp + m,1*Mp + m) = -UTHIF(m,m) / UAHIF(m,m);
        % ---------------------------------------------------------------------

    end  

    % ------------------------- Update Newton solver vector of unknowns
    DX = linsolve(J,-R);

    if any(isnan(DX))  || ~isreal(DX) || any(abs(DX) > 10^4)
         UASTAR = 0*RC;
         UTSTAR = 0*RC;
         TANBIC = 0*RC;
         UAHIF  = zeros(Mp,Mp);
         UTHIF  = zeros(Mp,Mp);

         disp('<WARNING>')
         disp('<WARNING> DX == NaN or imaginary... crash avoided...')
         disp('<WARNING>')
         N_iter      = 999;
    else
        TANBIC  = TANBIC + relax*DX( 1:Mp      )';
        UASTAR  = UASTAR + relax*DX((1:Mp)+  Mp)';
        UTSTAR  = UTSTAR + relax*DX((1:Mp)+2*Mp)';        
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
         N_iter         = 999;
         converged_flag = 0;
    end      
    
    % -----------------------------------------------------------------     
    
%     save temp
    
    % -------------------------------------- Prepare for the next iteration
    N_iter  = N_iter + 1;              % iteration in the N loop

    ERROR    = [ TANBIC - TANBIClast, ...
                 UASTAR - UASTARtemp, ...
                 UTSTAR - UTSTARtemp] ./ ...
               [max(abs(TANBIC),1e-4), ...
                max(abs(UASTAR),1e-4), ...
                max(abs(UTSTAR),1e-4)]; 


        if (Hvel ~=0)  & (Hbeta ~=0)
            % ------------------------------- Uncomment to plot UASTAR & UTSTAR
            figure(Hvel),
            ylabel('RED = UASTAR, BLUE = UTSTAR')

                    if     mod(N_iter-1,2) == 1
                        plot(RC,UASTAR,'r.-'), 
                    else
                        plot(RC,UASTAR,'m.-'), 
                    end               
                
                    if     mod(N_iter-1,2) == 1
                        plot(RC,UTSTAR,'b.-'), 
                    else
                        plot(RC,UTSTAR,'c.-'), 
                    end 
            % ---------------------------------------- Uncomment to plot TANBIC
            figure(Hbeta),
            ylabel('TANBIC')
                if     mod(N_iter-1,4) == 1
                    plot(RC,TANBIC,'r.-'), hold on,
                elseif mod(N_iter-1,4) == 2
                    plot(RC,TANBIC,'m.-'), hold on,
                elseif mod(N_iter-1,4) == 3
                    plot(RC,TANBIC,'b.-'), hold on,
                else
                    plot(RC,TANBIC,'k.-'), hold on,
                end

            % --------------------------- Uncomment to pause between iterations     
            pause(0.1),
        end
end                                              % (END WHILE LOOP WA1)
    
    if N_iter > ITER
        disp('WARNING: Align_wake_Newton.m did NOT converge.'),
        converged_flag = 0;
    else
        converged_flag = 1;
    end

%     disp(['[Align_wake_Newton] The max Residual for iteration ',num2str(N_iter-1),' is: ',num2str(max(ERROR))]),
% ================================================= END Align_wake Function
% =========================================================================