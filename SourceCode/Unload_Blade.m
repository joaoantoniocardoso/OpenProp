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
% ========================================================== Unload_Blade.m
%
% This code unloads the hub and tip as specified by HUF and TUF.
%
%
% -------------------------------------------------------------------------

function     [G,UASTAR,UTSTAR,TANBIC,UARING,URRING,Gd,UADUCT] = ...
                                    Unload_Blade(HUF,TUF,RC,Rhub_oR, G,  VAC,VTC, TANBIC,RV,DR,L,Mp,Z, ...
                                                 Hub_flag,ITER,Bsmooth,CD,CoD,Js,VMIV,Rhv,CTPdes,...
                                                 Duct_flag,UADUCT,XdRING,Rduct_oR,VARING,GdRING,UADIF,Gd,CDd,CTDdes,...
                                                 DAHIF_times_TANBIC,DRHIF_times_TANBIC,...
                                                 Plot_flag,Hgamma,HHgamma,Hvel,HHvel,Hbeta,HHbeta) 
                                             
    % ---------------------- Unload hub and tip as specified by HUF and TUF
    RU = (RC - Rhub_oR)./(1 - Rhub_oR);

    nH = 4;
    nT = 3;

    GH = HUF*G(1) *sqrt(1-RU.^2).*(1-RU.^2).^(2*nH-2);
    GT = TUF*G(Mp)*sqrt(1-RU.^2).*(  RU.^2).^(2*nT-2);

    G  = G - GH' - GT';

    % ------------------ Align the wake to the new circulation distribution
    [UASTAR,UTSTAR,TANBIC,UAHIF,UTHIF] =    ...
         Align_wake(G,  VAC,VTC, TANBIC,RC,RV,L,Mp,Z, ...
                    Hub_flag,Rhub_oR, Duct_flag,Rduct_oR,UADUCT, Hvel,Hbeta);    
                
                
    % ----------------------------------------------------------------- 
    if Duct_flag == 1 
        TANBICsmooth = TANBIC * Bsmooth;
        
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
    else
        UARING = 0;
        URRING = 0;
        Gd     = 0;
        UADUCT = 0;
    end
    % -----------------------------------------------------------------                 
    
    % ------- Iterate to scale G to get desired value of thrust coefficient
    CT_iter   = 1;                          % iteration in the CT loop
    CT_res    = 1;                          % residual CT
    CT_last2  = 0;                          % the CT prior to the last CT
    CT_last1  = 0;                          % the last value of CT
    GMF_last2 = 0;                          % the GMF prior to the last GMF
    GMF_last1 = 0;                          % the last value of GMF

    while CT_iter <= ITER & any(CT_res > 0.01)           % (WHILE LOOP CT1)

        if CT_iter == 1
            GMF = 1;

        elseif CT_iter == 2
            GMF = 1+(CTPdes-CT)/(5*CTPdes);

        elseif CT_iter > 2
            GMF = GMF_last1 + (GMF_last1-GMF_last2)*(CTPdes-CT_last1)/...
                              ( CT_last1- CT_last2);
        end

        G = GMF.*G;                   % GMF = G Multiplication Factor
        
        % ------------------------------------- Evaluate induced velocities
        UASTAR = (UAHIF*G)';  
        UTSTAR = (UTHIF*G)';    
        
                
        % ----------------------------------------------------------------- 
        if Duct_flag == 1            
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
        

        % ------------- Compute thrust & torque coefficients and efficiency
        [CT,CQ,CP,KT,KQ, CTH,TAU, Ja,Jw,VMWV, EFFYo, EFFY,ADEFFY,QF] = ...
          Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,UADUCT,CD,CoD,G,Z,Js,VMIV,Hub_flag,Rhub_oR,Rhv,CTDdes);
      
      
        % ---------------------------------- Prepare for the next iteration
        CT_iter   = CT_iter + 1;            % iteration in the CT loop
        CT_res    = abs( ( (CT-CTDdes) - CTPdes )./CTPdes ); % residual CT
        CT_last2  = CT_last1;               % the CT prior to the last CT
        CT_last1  = CT;                     % the last value of CT
        GMF_last2 = GMF_last1;              % the GMF prior to the last GMF
        GMF_last1 = GMF;                    % the last value of GMF   

    end                                              % (END WHILE LOOP CT1)


%     disp(' '),
%     disp('Forces after hub and tip unloading:')
%     disp(['CT   = ',num2str(CT)]),   
%     disp(['CQ   = ',num2str(CQ)]),   
%     disp(['CP   = ',num2str(CP)]),   
%     disp(['EFFY = ',num2str(EFFY)]),
%     disp(' ')  
    
    % ---------------------------------------------------------------------
    if Plot_flag == 1
        % ---------------------------------------------------------------------    
        figure( Hgamma),
           set(HHgamma,'linestyle','--'),
               HHgamma = plot(RC,G,'k.-');
        % ---------------------------------------------------------------------     

        % --------------------------------------------------------------------- 
        figure( Hvel),
           set(HHvel,'linestyle','--'),
        if Duct_flag == 1
               HHvel = plot(RC,UASTAR,'k.-',RC,UTSTAR,'k.-',RC,UADUCT,'k.-');
        else
               HHvel = plot(RC,UASTAR,'k.-',RC,UTSTAR,'k.-');
        end
        % --------------------------------------------------------------------- 
        
        % --------------------------------------------------------------------- 
        figure( Hbeta),
           set(HHbeta,'linestyle','--'),
               HHbeta = plot(RC,TANBIC,'k.-');
        % ---------------------------------------------------------------------         

        pause(0.0001),   drawnow,
    end  % if Plot_flag == 1     
    % ---------------------------------------------------------------------
    
   
end                                                
% ------------------------------------------------------------------------- 