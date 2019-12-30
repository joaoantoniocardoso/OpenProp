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


function [] = Make_3D_Blade_Image(X3D,Y3D,Z3D,Duct_flag,Rduct_oR,Cduct_oR,Xduct_oR, Hub_flag,Rhub_oR,  Js,BetaI_c,theta,TANBIV,RV,Rm, GUI_flag,Plots,PlotPanels)

if nargin <= 12
    Rm = 1;
    GUI_flag = 0;
elseif nargin == 13
    GUI_flag = 0; 
end

% -------------------------------------------------------------------------
Mp       = size(X3D,1)-1;
Np       = size(X3D,2)/2;
Z        = size(Y3D,3);
theta_Z  = 0:360/Z:360;                   % angle between blades [deg]
R        =   Rm;
D        = 2*R;
Rhub     =   Rhub_oR*Rm;
Dhub     = 2*Rhub;
Rduct    =   Rduct_oR*Rm;
Cduct    =   Cduct_oR*Rm;
Xduct    =   Xduct_oR*Rm;
Dduct    = 2*Rduct;
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
    if GUI_flag
        
        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotPanels(14));
        axes(h);
        
    else  
        Fig3_S = figure('position',[450 400 800 600],'name','Propeller Image','numbertitle','off');
        
    end
    
% -------------------------------------------------------------------------
    set(gca,'Position',[-0.025 -0.05 1.1 1.1]),  % zoom in
    hold on;
    axis equal;
    xlabel('X (3D) [m]','FontSize',12);
    ylabel('Y (3D) [m]','FontSize',12);
    zlabel('Z (3D) [m]','FontSize',12);
    title('3D Propeller Image','FontSize',16);
% -------------------------------------------------------------------------
    


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
    % ------------------------------------------ Plot the propeller surface
    for k = 1:Z
        surf(X3D(:,:,1),Y3D(:,:,k),Z3D(:,:,k));
    end

    colormap gray;
    shading interp;
    grid on;
    
    if Duct_flag == 0
        axis([-0.6*R R -1.05*R 1.05*R -1.05*R 1.05*R]);
    else
        axis([    -R R -1.25*R 1.25*R -1.25*R 1.25*R]);       %modified for duct
    end


    % -------------------------------------------------------- Plot the hub
    
%         % -------------------------------------------------------- Plot the hub
%         xxx = [0:0.05:1]; xxx = xxx*2;
%         [yh0,zh0,xh0] = cylinder(Rhub*sqrt(1-xxx.^2/2^2),50);  
%         xh0 = 2*Rhub*xh0 + 0.25*R;
%         surf(xh0,yh0,zh0);
% 
%         [yh1,zh1,xh1] = cylinder(Rhub,50);
%         xh1 = -6*Rhub*xh1 + 0.25*R;
%         surf(xh1,yh1,zh1);      

    
    if Hub_flag == 1
        Lhub = 1.25*Dhub;

        % tick = 90:-15:0;
        % [yh0,zh0,xh0] = cylinder(Rhub*sind(tick),50);
        
        xxx = [0:0.05:1]; xxx = xxx*2;
        [yh0,zh0,xh0] = cylinder(Rhub*sqrt(1-xxx.^2/2^2),50);
        xh0a = -2*Rhub*xh0 - Rhub;
        surf(xh0a,yh0,zh0);        
        
        xh0b = +2*Rhub*xh0 + Rhub;
        surf(xh0b,yh0,zh0);

        [yh1,zh1,xh1] = cylinder(Rhub,50); % xh1 = [0,1]
        xh1 = Lhub*xh1 - Rhub;           % xh1 = [-Rhub,c(1)-Rhub]
        surf(xh1,yh1,zh1);
    end
    
    
    % ----------------- Plot the suction side (green) & pressure side (red)
    for i = 1:Mp+1          % for each section along the span
        for k = 1:Z       % for each blade
            plot3(X3D(i,1:Np,1),     Y3D(i,1:Np,k),     Z3D(i,1:Np,k),     'g','Linewidth',1); % suction surface
            plot3(X3D(i,Np+1:2*Np,1),Y3D(i,Np+1:2*Np,k),Z3D(i,Np+1:2*Np,k),'r','Linewidth',1); % pressure surface
        end
    end

    for j = 1:Np          % for each point along the chord
        for k = 1:Z       % for each blade
            plot3(X3D(:,j,1),   Y3D(:,j,k),   Z3D(:,j,k),   'g','Linewidth',1); % suction surface
            plot3(X3D(:,j+Np,1),Y3D(:,j+Np,k),Z3D(:,j+Np,k),'r','Linewidth',1); % pressure surface
        end
    end
    
    % -------------------------------------------------- Plot the tip black
    i = Mp+1; % tip section
    
    for k = 1:Z
        for j = 1:Np-2
                plot3([X3D(i,1+j,k), X3D(i,2*Np-j,k)],...
                      [Y3D(i,1+j,k), Y3D(i,2*Np-j,k)],...
                      [Z3D(i,1+j,k), Z3D(i,2*Np-j,k)],'k','Linewidth',1); % tip surface        
        end
    end

    % --------------------------------- Plot the leading and trailing edges
    for k = 1:Z           % for each blade
        plot3(X3D(:,1,1), Y3D(:,1,k), Z3D(:,1,k), 'k','Linewidth',2); % trailing edge
        plot3(X3D(:,Np,1),Y3D(:,Np,k),Z3D(:,Np,k),'k','Linewidth',2); % leading edge
    end

    
%     % ------------------------------------------------ Plot the chord lines
%     for i = 1:Mp+1          % for each section along the span
%         for k = 1:Z         % for each blade
%             plot3(X3D(i,[1,Np],1),     Y3D(i,[1,Np],k),     Z3D(i,[1,Np],k),     'k--','Linewidth',1); % suction surface
%         end
%     end
    
%     % ------------------------------------------ Plot the coordinate system
%     % Axes
%     plot3([0 R],[0 0],[0 0],'y','LineWidth',2),
%     plot3([0 0],[0 R],[0 0],'r','LineWidth',2),
%     plot3([0 0],[0 0],[0 R],'b','LineWidth',2),
% 
%     % Circle at the X = 0 location on the hub
%     phi = 0:0.01:2*pi;
%     Xhc =   zeros(size(phi));
%     Yhc = - Rhub * sin(phi);
%     Zhc =   Rhub * cos(phi);
%     plot3(Xhc,Yhc,Zhc,'y','LineWidth',2),
% 
%     % Propeller reference line (i.e. the directrix)
%     for k = 1:Z
%         PRL(:,k) = [1,                0,                 0; ...
%                     0, cosd(theta_Z(k)), -sind(theta_Z(k)); ...
%                     0, sind(theta_Z(k)),  cosd(theta_Z(k))]*[0; 0; R];
% 
%         plot3([0, PRL(1,k)],[0, PRL(2,k)],[0, PRL(3,k)],'y--','LineWidth',1)
%     end

%     % ---------------------------------------------- Plot propeller helices
%     % Advance coefficient helix 1, black
%     phi = 0:0.01:pi/4;
%     thetaH = atan((Js/pi)*(R/Rhub));
%     Xh0 =   Rhub * phi * tan(thetaH);
%     Yh0 = - Rhub * sin(phi);
%     Zh0 =   Rhub * cos(phi);
% 
%     % BetaI angle helix 2, green
%     phi = 0:0.01:pi/4;
%     thetaH = BetaI_c(1)*pi/180;
%     Xh2 =   Rhub * phi * tan(thetaH);
%     Yh2 = - Rhub * sin(phi);
%     Zh2 =   Rhub * cos(phi);
% 
%     % Pitch angle helix 3, blue
%     phi = 0:0.01:pi/4;
%     thetaH = theta(1)*pi/180;
%     Xh3 =   Rhub * phi * tan(thetaH);
%     Yh3 = - Rhub * sin(phi);
%     Zh3 =   Rhub * cos(phi);
% 
%     plot3(Xh0,Yh0,Zh0,'k','LineWidth',2),
%     plot3(Xh1,Yh1,Zh1,'r','LineWidth',2),
%     plot3(Xh2,Yh2,Zh2,'g','LineWidth',2),
%     plot3(Xh3,Yh3,Zh3,'b','LineWidth',2),
%    
% % ---------------------------------------------- Plot the trailing vortices    
%     % Beta angle helix at each votex point (each trailing vortex)
%     BetaI_v = atand(TANBIV);
%     
%     for m = 1:Mp+1
%         phi = 0:0.01:2*pi;
%         thetaH = BetaI_v(m)*pi/180;
%         Xh4 = - RV(m)*R * phi * tan(thetaH);
%         Yh4 =   RV(m)*R * sin(phi);
%         Zh4 =   RV(m)*R * cos(phi);
%         
%         plot3(Xh4,Yh4,Zh4,'g','LineWidth',2),
%     end     
% %     
% %     % Beta angle image helix for each spanwise section
% %     for m = 1:Mp+1
% %         RVW   = Rhub_oR^2/RV(m);
% %         TANBW = TANBIV(1)*RV(1)/RVW;
% %         phi = 0:0.01:2*pi;
% %         thetaH = atand(TANBW)*pi/180;
% %         Xh4 = - RVW*R * phi * tan(thetaH);
% %         Yh4 =   RVW*R * sin(phi);
% %         Zh4 =   RVW*R * cos(phi);
% %         
% %         plot3(Xh4,Yh4,Zh4,'--r','LineWidth',2),
% %     end 
% %     
% % --------------------------------------------- Plot the horseshoe vortices    
%     % Beta angle helix at each votex point (each trailing vortex)
%     BetaI_v = atand(TANBIV);
%     dR = 0.005*R;
%     
%     for m = 1:Mp
%         phi = 0:0.01:2*pi;
%         thetaH = BetaI_v(m)*pi/180;
%         Xh4 = - (RV(m)+dR)*R * phi * tan(thetaH);
%         Yh4 =   (RV(m)+dR)*R * sin(phi);
%         Zh4 =   (RV(m)+dR)*R * cos(phi);
%         
%         plot3(Xh4,Yh4,Zh4,'g','LineWidth',2),
%         
%         thetaH = BetaI_v(m+1)*pi/180;
%         Xh4 = - (RV(m+1)-dR)*R * phi * tan(thetaH);
%         Yh4 =   (RV(m+1)-dR)*R * sin(phi);
%         Zh4 =   (RV(m+1)-dR)*R * cos(phi);
%         
%         plot3(Xh4,Yh4,Zh4,'g','LineWidth',2),
%     end     
%     
%     % Beta angle image helix for each spanwise section
%     for m = 1:Mp+1
%         RVW   = Rhub_oR^2/RV(m);
%         TANBW = TANBIV(1)*RV(1)/RVW;
%         phi = 0:0.01:2*pi;
%         thetaH = atand(TANBW)*pi/180;
%         Xh4 = - RVW*R * phi * tan(thetaH);
%         Yh4 =   RVW*R * sin(phi);
%         Zh4 =   RVW*R * cos(phi);
%         
%         plot3(Xh4,Yh4,Zh4,'--r','LineWidth',2),
%     end     

% ---------------------------------------------- Plot duct
%     if Duct_flag == 1
%         Duct_Ang=0;
%         fo = 0;
%         to = 0;
%         Duct_Ang = 0 * pi/180;
%         
%         colormap(jet)
%         Duct_Plot(Rduct,Cduct,fo,to,Duct_Ang,0.5-Xduct/Cduct)
% %       Duct_Plot(Rduct,Rduct,0,0,Duct_Ang*pi/180,0.5)
%         %axis equal
%     end
    

    
    view(-58,8)
    set(gca,'XTickLabel',{''},'YTickLabel',{''},'ZTickLabel',{''})
    set(gca,'TickLength',[0 0])
    xlabel(''), ylabel(''), zlabel(''), title('')
    grid off, axis off
    
    