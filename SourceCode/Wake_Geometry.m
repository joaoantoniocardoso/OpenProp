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
% =========================================== Update Wake Geometry Function
%
% This function updates the geometry of the wake.  The wake trailing
% vortices are assumed to lie along streamlines.  The output are the
% positions of the endpoints of each segment of the wake.
%
%
% This function evaluates the influence of a unit strength vortex line.
%
% [XH,YH ZH] = equal length vectors of vortex segment endpoint positions
% Xpf        = [x,y,z] position of field point at which to evaluate velocities
% epsilon    = vortex core size
% Z          = number of propeller blades
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% PLACEHOLDER: INSERT CODE TO CONSTRUCT WAKE FOLLOWING STREAMLINES
%
% function [WX,WY,WZ] = Wake_Geometry(UASTAR,UTSTAR,URSTAR,Vs,Js,RV,Z,WX,WY,WZ,epsilon,G)
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function [WX,WY,WZ] = Wake_Geometry(Mp,RC,RV,TANBIV,VAC,UASTAR,URSTAR,CTP)
% function [WX,WY,WZ] = Wake_Geometry(Mp,RV,TANBIV,C,Ns)
% function [WX,WY,WZ] = Wake_Geometry(Mp,RV,TANBIV)


%%
% RV = 0.2:0.1:1;
% TANBIV = tand(55)*0.2./RV;
% Mp = length(RV)-1;
% Z  = 3;
% Ns = 6;


% C  = 2;  % number of coils    in the helix
% Ns = 30; % number of segments in the helix

C  = 4;  % number of coils    in the helix
Ns = 80; % number of segments in the helix

WX = zeros(Mp+1,Ns+1);
WY = zeros(Mp+1,Ns+1);
WZ = zeros(Mp+1,Ns+1);
  
% ------------------------------------------------ Construct key wake helix
n   = 1:Ns+1;                 % index of segment endpoint
phi = C*2*pi*(n-1).^2/Ns^2;   % angular parameter, 2*pi = 1 revolution


% % -------------------------------------------------------------------------
% % -------------------------------------------------------------------------
% % Constant radius envelope:
% RE = RV(end)+0*phi; 
% 
% % Constat-radius helix
% for m = 1:Mp+1
%     WX(m,:) = - RV(m) * phi * TANBIV(m);   % x/R position
%     WY(m,:) =   RV(m) * sin(phi);          % y/R position
%     WZ(m,:) =   RV(m) * cos(phi);          % z/R position
% end 



% % -------------------------------------------------------------------------
% % -------------------------------------------------------------------------
% % Radius envelope: cubic function of x/R
% % RE = a0 + a1*RX + a2*RX^2 + a3*RX^3
% XU = 2;        % Say envelope reaches ultimate radius at XU == 1*D / R == 2
% RU = sqrt((1+sqrt(1+CTP))/(2*sqrt(1+CTP)));
% 
% a0 = RV(end);
% a1 = URSTAR(end)/(VAC(end)+UASTAR(end));
% a2 = (3*RU - 3*a0 - 2*a1*XU          )/XU^2;
% a3 = (  RU -   a0 -   a1*XU - a2*XU^2)/XU^3;
% 
% RX = [0:0.01:XU];                            % x/R of envelope
% RE = [a0 + a1*RX + a2*RX.^2 + a3*RX.^3,RU];  % r/R of envelope
% RX = [RX,100];                               % x/R of envelope
% 
% % Cubic radius helix
% for m = 1:Mp+1
%     WX(m,:) = - RV(m) * phi * TANBIV(m);   % x/R position
%     WY(m,:) =   interp1(RX,RE,abs(WX(m,:))) .* sin(phi);          % y/R position
%     WZ(m,:) =   interp1(RX,RE,abs(WX(m,:))) .* cos(phi);          % z/R position
% end 



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Actuator disk theory ultimate radius:
XU = 2;        % Say envelope reaches ultimate radius at XU == 1*D / R == 2
RU = sqrt((1+sqrt(1+CTP))/(2*sqrt(1+CTP)));

% Find wake radii by conservation of mass:
Qin  = (VAC+UASTAR).*(RV(2:Mp+1).^2 - RV(1:Mp).^2); % Q/pi for each annulus
Qtot = sum(Qin);                 % total flow rate / pi
Vout = Qtot / (RU^2 - RV(1)^2);  % axial outflow velocity

RO       = 0*RV;
RO(Mp+1) = RU;
for m = Mp:-1:1
    RO(m) = sqrt(RO(m+1)^2 - Qin(m)/Vout); 
end


% Radius envelope: cubic function of x/R
% RE = a0 + a1*RX + a2*RX^2 + a3*RX^3
a0 = RV;
a1 = interp1(RC,URSTAR./(VAC+UASTAR),RV,'linear','extrap');
a2 = (3*RO - 3*a0 - 2*a1*XU          )/XU^2;
a3 = (  RO -   a0 -   a1*XU - a2*XU^2)/XU^3;

RX = [0:0.01:XU];                                      % x/R of envelope
RE = [a0'*ones(size(RX)) + a1'*RX + a2'*RX.^2 + a3'*RX.^3];  % r/R of envelope 
RX = [RX,100];                                         % x/R of envelope
RE = [RE,RU*ones(size(a0'))];  % RE(m,:) corresponds to RV(m)

% Cubic radius helix
for m = 1:Mp+1
    WX(m,:) = - RV(m) * phi * TANBIV(m);   % x/R position
    WY(m,:) =   interp1(RX,RE(m,:),abs(WX(m,:))) .* sin(phi);          % y/R position
    WZ(m,:) =   interp1(RX,RE(m,:),abs(WX(m,:))) .* cos(phi);          % z/R position
end 





% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% CONSTRUCT WAKE WITH CONSTANT RADIUS HELICES
%
% DR = 0*TANBIV;
%DR = (C*2*pi/Ns)*TANBIV.*pchip(RC,(URSTAR./(VAC+UASTAR)),RV);


% % Linearly-increasing-radius helix
% for m = 1:Mp+1  % for each helix
%     WX(m,:) = - (RV(m)+DR(m)*(C/4)*sqrt(phi/(C*2*pi))) .* phi * TANBIV(m);   % x/R position scaled by prop radius
%     WY(m,:) =   (RV(m)+DR(m)*(C/4)*sqrt(phi/(C*2*pi))) .* sin(phi);          % y/R position scaled by prop radius
%     WZ(m,:) =   (RV(m)+DR(m)*(C/4)*sqrt(phi/(C*2*pi))) .* cos(phi);          % z/R position scaled by prop radius
% end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


   
% %%
% % ----------------------------------------------------------- Plot the wake
% figure,
% 
% % ---------------------------------------------- Plot the coordinate system
% R = 1;
% Rhub = min(RV);
% Z = 3;
% 
% % Axes
% plot3([0 R],[0 0],[0 0],'m','LineWidth',2), hold on,
% plot3([0 0],[0 R],[0 0],'r','LineWidth',2),
% plot3([0 0],[0 0],[0 R],'b','LineWidth',2),
% 
% 
% 
% 
% text(R+0.1,0,0,'X')
% text(0,R+0.1,0,'Y')
% text(0,0,R+0.1,'Z')
% 
% 
% % Circle at the X = 0 location on the hub
% phic = 0:0.01:2*pi;
% Xhc  =   zeros(size(phic));
% Yhc  = - Rhub * sin(phic);
% Zhc  =   Rhub * cos(phic);
% plot3(Xhc,Yhc,Zhc,'g','LineWidth',2),
% 
% % Propeller reference line (PRL, or directrix)
% theta_Z = 0:360/Z:360;            % angle between blades [deg]
% 
% for k = 1:Z
%     PRL(:,k) = [1,                0,                 0; ...
%                 0, cosd(theta_Z(k)), -sind(theta_Z(k)); ...
%                 0, sind(theta_Z(k)),  cosd(theta_Z(k))]*[0; 0; R];
% 
%     plot3([0, PRL(1,k)],[0, PRL(2,k)],[0, PRL(3,k)],'g--','LineWidth',1)
% end
% % -------------------------------------------------------------------------
% 
% % ------------------------------------------ Plot the radius envelope, R(x)
% for m = 1:size(RE,1)
%     plot3(-RX,0*RX,RE(m,:),'r--')
% end
% 
% % ------------------------------------------ Set the axes
% % Rmax = sqrt(max([1;WY(:,1).^2+WZ(:,1).^2;WY(:,end).^2+WZ(:,end).^2]));
% Rmax = max(RU,sqrt(max([1;WY(:,1).^2+WZ(:,1).^2;WY(:,end).^2+WZ(:,end).^2])));
% 
% axis([min(min(WX)) 0 -Rmax Rmax -Rmax Rmax])
% 
% % ------------------------------------------------ Plot only the last helix
% plot3(WX(end,:),WY(end,:),WZ(end,:),'k','LineWidth',2),
% 
% % % --------------------------------------------------- Plot all wake helices
% % for m = 1:Mp+1
% %     plot3(WX(m,:),WY(m,:),WZ(m,:),'k','LineWidth',2), hold on,
% % end
% 
% xlabel('x'), ylabel('y'),zlabel('z')
% box on
% axis equal
% 
% 
% plot3(WX(end,:),WY(end,:),WZ(end,:),'b.','MarkerSize',20),
% plot3(WX(end,1),WY(end,1),WZ(end,1),'r.','MarkerSize',20),
% 
% for i = 1:floor(Ns/10)
%     plot3(WX(end,i*10+1),WY(end,i*10+1),WZ(end,i*10+1),'r.','MarkerSize',20),
% end
% 
% % ------------------------------------------- Plot all Z outer wake helices
% HPData = [WX(end,:);WY(end,:);WZ(end,:)]; % Helix position data, 
% % HPData is formatted such that:  XH(s) == HPData(1,s), and a helix can be
% % plotted using:
% %              plot3(HPData(1,:),HPData(2,:),HPData(3,:),'b','linewidth',2) 
% 
% 
% theta_Z = 0:360/Z:360*(Z-1)/Z;            % angle between blades [deg]
%              
% for k = 1:Z % for each blade
%     
%     Rotation = [1,                0,                 0; ...
%                 0, cosd(theta_Z(k)), -sind(theta_Z(k)); ...
%                 0, sind(theta_Z(k)),  cosd(theta_Z(k))];
%     
%             
%     HPDr = Rotation*HPData;  % rotated helix position data
%     
%     plot3(HPDr(1,:),HPDr(2,:),HPDr(3,:),'r','linewidth',2)
% end
%     
% % % ---------------------------- Plot a highly resolved constant radius helix
% % Ns2 = 100*C;
% % WX2 = zeros(Mp+1,Ns2+1);
% % WY2 = zeros(Mp+1,Ns2+1);
% % WZ2 = zeros(Mp+1,Ns2+1);
% %   
% % n2   = 1:Ns2+1;              % index of segment endpoint
% % phi2 = C*2*pi*(n2-1)/Ns2;   % angular parameter, 2*pi = 1 revolution
% % 
% % % Constat-radius helix
% % for m = 1:Mp+1
% %     WX2(m,:) = - RV(m) * phi2 * TANBIV(m);   % x/R position scaled by prop radius
% %     WY2(m,:) =   RV(m) * sin(phi2);          % y/R position scaled by prop radius
% %     WZ2(m,:) =   RV(m) * cos(phi2);          % z/R position scaled by prop radius
% % end 
% % 
% % plot3(WX2(end,:),WY2(end,:),WZ2(end,:),'k--','MarkerSize',20),
% 
% % % ---------------------------- Plot a highly resolved cubic radius helix
% % WX3 = zeros(Mp+1,Ns2+1);
% % WY3 = zeros(Mp+1,Ns2+1);
% % WZ3 = zeros(Mp+1,Ns2+1);
% % 
% % n2   = 1:Ns2+1;              % index of segment endpoint
% % phi2 = C*2*pi*(n2-1)/Ns2;   % angular parameter, 2*pi = 1 revolution
% %
% % % Cubic radius helix
% % for m = 1:Mp+1
% %     WX3(m,:) = - RV(m) * phi2 * TANBIV(m);                   % x/R position
% %     WY3(m,:) =   interp1(RX,RE(m,:),abs(WX3(m,:))) .* sin(phi2);   % y/R position
% %     WZ3(m,:) =   interp1(RX,RE(m,:),abs(WX3(m,:))) .* cos(phi2);   % z/R position
% % 
% %     plot3(WX3(m,:),WY3(m,:),WZ3(m,:),'k--','MarkerSize',20),
% % end
% % 
% % plot3(WX3(end,:),WY3(end,:),WZ3(end,:),'k--','MarkerSize',20),
% 
% 
% 
% 
% 
% % ------------------------------------------ Set the axes
% % Rmax = sqrt(max([1;WY(:,1).^2+WZ(:,1).^2;WY(:,end).^2+WZ(:,end).^2]));
% Rmax = max(RU,sqrt(max([1;WY(:,1).^2+WZ(:,1).^2;WY(:,end).^2+WZ(:,end).^2])));
% 
% axis([min(min(WX)) 0 -Rmax Rmax -Rmax Rmax])

