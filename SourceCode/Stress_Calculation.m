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


function [MAX_STRESS, MAX_n_dot, blade_stress] = Stress_Calculation(rho, rho_p, allowable_stress, Vs, Vs_dot, n, n_dot, R, Mp, Np, RC, DR, CD, CL, VSTAR, CoD, BetaIC, theta, ...
                        Ixc, Iyc, Ixyc, Area, Xbar, Ybar, xl, yl, xu, yu)

omega     = 2*pi*n;      % [rad/s]   angular rotation rate
omega_dot = 2*pi*n_dot;  % [rad/s^2] angular acceleration                    
                    
sin_BetaIC = sin(BetaIC);
cos_BetaIC = cos(BetaIC);

D = 2*R; % [m] diameter

MQ   = zeros(1,Mp);
MT   = zeros(1,Mp);
Mxo  = zeros(1,Mp);
Myo  = zeros(1,Mp);

blade_stress  = zeros(Mp,2*Np);

xs = cat(2,xu,fliplr(xl));
ys = cat(2,yu,fliplr(yl));




% -------------------------------------------------------------------------
% Difference in radii between blade sections (for outboard sections only):
Rdif = zeros(Mp,Mp);
for i = 1:Mp;
    Rdif(i,i:Mp) = RC(i:Mp) - RC(i);
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Centrifugal force, [N]
Fc   = zeros(1,Mp);

for m = 1:Mp
   Fc(m)  = sum( rho_p*Area(m:Mp).*(DR(m:Mp)*R) .* omega^2.*(RC(m:Mp)*R)  ); % [N]
end
% -------------------------------------------------------------------------



% figure
% axes('fontweight','bold')
% hold on
% plot(xs(1,:),ys(1,:),'-ks','LineWidth',2)
% axis equal
% title('Root Section Ordinates','Fontweight','bold', 'Fontsize', 14)
% xlabel('X (m)','Fontweight','bold','Fontsize',12)
% ylabel('Y (m)','Fontweight','bold','Fontsize',12)


%   for m = 1:1   %  ASSUME THAT MAX STRESS IS AT ROOT  
for m = 1:Mp  %  COMPUTE STRESS FOR ENTIRE BLADE
    
    % moments in the "torque" and "thrust" directions due to lift and drag forces
    MQ(m) = sum(  rho*Vs^2*R^3 * (VSTAR.^2 .* (CL .* sin_BetaIC + CD .* cos_BetaIC).* CoD .* Rdif(m,:).*DR)  );
    MT(m) = sum(  rho*Vs^2*R^3 * (VSTAR.^2 .* (CL .* cos_BetaIC - CD .* sin_BetaIC).* CoD .* Rdif(m,:).*DR)  );

    % Moments about local blade (x,y) axes
    Mxo(m) = MT(m)*cos(theta(m)) + MQ(m)*sin(theta(m));
    Myo(m) = MT(m)*sin(theta(m)) - MQ(m)*cos(theta(m));
    
    
    
    
            % additional moments due to added mass forces
            % MQam = sum( (pi/4)*rho*(CoD(m)*D)^2 * (omega_dot * RC*R .* sin(theta) - Vs_dot * cos(theta) ) .* DR*R  .* sin(theta) .* Rdif(m,:)*R );
            % MTam = sum( (pi/4)*rho*(CoD(m)*D)^2 * (omega_dot * RC*R .* sin(theta) - Vs_dot * cos(theta) ) .* DR*R  .* cos(theta) .* Rdif(m,:)*R );
    
    % additional moments due to added mass forces, assuming Vs_dot is small
    MQam_per_omega_dot(m) = sum( (pi/4)*rho*(CoD(m)*D)^2 * RC*R .* sin(theta) .* DR*R  .* sin(theta) .* Rdif(m,:)*R );  % NOTE:  MQam == MQam_per_omega_dot * omega_dot;
    MTam_per_omega_dot(m) = sum( (pi/4)*rho*(CoD(m)*D)^2 * RC*R .* sin(theta) .* DR*R  .* cos(theta) .* Rdif(m,:)*R );  % NOTE:  MTam == MTam_per_omega_dot * omega_dot;
    
    % additional moment due to inertial acceleration
    MQin_per_omega_dot(m) = sum( rho_p * Area .* RC*R .* DR*R .* Rdif(m,:)*R );    % NOTE:  MQin == MQin_per_omega_dot * omega_dot;
    
    % Moments about local blade (x,y) axes
    Mx_wdot_per_omega_dot(m) = MTam_per_omega_dot(m)*cos(theta(m)) + (MQam_per_omega_dot(m)+MQin_per_omega_dot(m))*sin(theta(m));  % NOTE:  Mx_wdot == Mx_wdot_per_omega_dot * omega_dot;
    My_wdot_per_omega_dot(m) = MTam_per_omega_dot(m)*sin(theta(m)) - (MQam_per_omega_dot(m)+MQin_per_omega_dot(m))*cos(theta(m));  % NOTE:  My_wdot == My_wdot_per_omega_dot * omega_dot;
    

    
    xsdiff = xs(m,:) - Xbar(m);     
    ysdiff = ys(m,:) - Ybar(m);
    
    
    % Stress (positive indicates tension)
    blade_stress(m,:) = -Myo(m)*xsdiff/Iyc(m) - Mxo(m)*ysdiff/Ixc(m) + Fc(m)/Area(m)  - omega_dot*My_wdot_per_omega_dot(m)*xsdiff/Iyc(m) - omega_dot*Mx_wdot_per_omega_dot(m)*ysdiff/Ixc(m) ;
    
  
  
%     figure(m),
%         patch(xsdiff(m,:),ysdiff(m,:),blade_stress(m,:))
%         colormap(jet)
%         colorbar
%         grid on
%         axis equal
end


% figure, plot(RC,MT,'.-'), title('MT')
% figure, plot(RC,MQ,'.-'), title('MQ')
% figure, plot(RC,Mxo,'.-'), title('Mx')
% figure, plot(RC,Myo,'.-'), title('My')
% 
% 
% figure, hold on,
%     plot3(xs(1,:),ys(1,:),  blade_stress(1,:),'-k','LineWidth',2)
%     plot3(xs(1,:),ys(1,:),0*blade_stress(1,:),'-k','LineWidth',2)
%     fill(xs(1,:),ys(1,:),'y')
%     axis square

% -------------------------------------------------------------------------
% Find the maximum blade stress:
MAX_STRESS = max(blade_stress(:));



% -------------------------------------------------------------------------
% Find maximum omega_dot such that the max stress is less than the allowable stress (i.e. some fraction of the yield stress)
%
%   blade_stress(m,:) = -Myo(m)*xsdiff/Iyc(m) - Mxo(m)*ysdiff/Ixc(m) + Fc(m)/Area(m) + ( - My_wdot_per_omega_dot(m)*xsdiff/Iyc(m) - Mx_wdot_per_omega_dot(m)*ysdiff/Ixc(m) ) * omega_dot;

[m,J] = find( blade_stress == max(blade_stress(:)));  % m == blade section


xsdiff = xs(m,J) - Xbar(m);     
ysdiff = ys(m,J) - Ybar(m);
 

              %    blade_stress(m,:) = -Myo(m)*xsdiff/Iyc(m) - Mxo(m)*ysdiff/Ixc(m) + Fc(m)/Area(m)   +    ( - My_wdot_per_omega_dot(m)*xsdiff/Iyc(m) - Mx_wdot_per_omega_dot(m)*ysdiff/Ixc(m) ) * omega_dot ;
MAX_OMEGA_DOT =  ( allowable_stress - (-Myo(m)*xsdiff/Iyc(m) - Mxo(m)*ysdiff/Ixc(m) + Fc(m)/Area(m)) )  /  ( - My_wdot_per_omega_dot(m)*xsdiff/Iyc(m) - Mx_wdot_per_omega_dot(m)*ysdiff/Ixc(m) )  ;


MAX_n_dot =  MAX_OMEGA_DOT / (2*pi);   

end