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


function [xl,xu,yl,yu, X3D,Y3D,Z3D] = Stress_BladeGeometry(D,RV,XR,XCoD,Xt0oD,Xf0oc,skew0,rake0,XPoD,Mp,Np,Meanline,Thickness)

% For blade stress calculations, use RG = RV
RG = RV; 


% -------- Thickness, rake, chord, and radius scaled for the model diameter
CoD  = InterpolateChord(XR,XCoD ,RG);
t0   =            pchip(XR,Xt0oD,RG)*D;
f0oc =            pchip(XR,Xf0oc,RG);
skew =            pchip(XR,skew0,RG);        % [deg], angular translation along mid-chord helix
rake =            pchip(XR,rake0,RG)*D;     % [m],   translation along propeller axis (3D X-axis)
PoD  =            pchip(XR,XPoD ,RG);

theta    = atand(PoD./(pi*RG));   % Nose-tail pitch angle, [deg]

% -------------------------------------------------------------------------
R  = D/2;  % [m] propeller radius

c  = CoD.*D;                          % section chord at the RG sections [m]
r  = RG .*R;                          % radius of        the RG sections [m]
% -------------------------------------------------------------------------

% ---------------------------------------- Lay out the 2D coordinate system
% x0   [ ], x/c distance along mid-chord line to interpolate geometry data.
% x1   [m], x   distance along mid-chord line to interpolate geometry data.
%               By definition, x1 == c/2 - c*x0.

%               At the Leading  Edge: x1 =  c/2, x0 = 0
%               At the Trailing Edge: x1 = -c/2, x0 = 1
%
x0 = zeros(   1,Np);
x1 = zeros(Mp+1,Np);

for j = 1:Np                               % for each point along the chord
    x0(j) = 0.5*(1-cos(pi*(j-1)/(Np-1)));  % [0 : 1], Cosine spacing along the chord
end

for i = 1:Mp+1                     % for each radial section along the span         
    x1(i,:) = c(i)/2 - c(i)*x0;    % lifting line at mid-chord
end
% -------------------------------------------------------------------------
    
% -------------------------------------------------------------------------
% ---------------------- Find normalized 2D foil geometry (at x0 positions)
%   f0octilde = zeros(Mp+1,1);  % can either be scalars or given versus RG(1:Mp+1)
%    CLItilde = zeros(Mp+1,1);  % can either be scalars or given versus RG(1:Mp+1)
% alphaItilde = zeros(Mp+1,1);  % can either be scalars or given versus RG(1:Mp+1)
if iscell(Meanline)  % Assume Meanline is given versus XR
    
      f0octilde = zeros(length(XR),1);  % can either be scalars or given versus RG(1:Mp+1)
       CLItilde = zeros(length(XR),1);  % can either be scalars or given versus RG(1:Mp+1)
    alphaItilde = zeros(length(XR),1);  % can either be scalars or given versus RG(1:Mp+1)
    fof0_temp        = zeros(length(XR),Np); % eventually given versus [RG(1:Mp+1),x/c(1:Np)]
    dfof0dxoc_temp   = zeros(length(XR),Np); % eventually given versus [RG(1:Mp+1),x/c(1:Np)]
    tot0_temp        = zeros(length(XR),Np); % eventually given versus [RG(1:Mp+1),x/c(1:Np)]    
    
    fof0        = zeros(Mp+1,Np); % eventually given versus [RG(1:Mp+1),x/c(1:Np)]
    dfof0dxoc   = zeros(Mp+1,Np); % eventually given versus [RG(1:Mp+1),x/c(1:Np)]
    tot0        = zeros(Mp+1,Np); % eventually given versus [RG(1:Mp+1),x/c(1:Np)]    
      
    
    for i = 1:length(XR)
        [f0octilde(i), CLItilde(i), alphaItilde(i), fof0_temp(i,:), dfof0dxoc_temp(i,:), tot0_temp(i,:)] = GeometryFoil2D(Meanline{i},Thickness{i},x0);
    end

      f0octilde = pchip(XR,   f0octilde, RG);
       CLItilde = pchip(XR,    CLItilde, RG);
    alphaItilde = pchip(XR, alphaItilde, RG);
    
      
    for j = 1:Np
             fof0(:,j) = pchip(XR,      fof0_temp(:,j), RG');
        dfof0dxoc(:,j) = pchip(XR, dfof0dxoc_temp(:,j), RG');
             tot0(:,j) = pchip(XR,      tot0_temp(:,j), RG');
    end

else
    [f0octilde, CLItilde, alphaItilde, fof0_temp, dfof0dxoc_temp, tot0_temp] = GeometryFoil2D(Meanline,Thickness,x0);
    
    fof0        = zeros(Mp+1,Np); % eventually given versus [RG(1:Mp+1),x/c(1:Np)]
    dfof0dxoc   = zeros(Mp+1,Np); % eventually given versus [RG(1:Mp+1),x/c(1:Np)]
    tot0        = zeros(Mp+1,Np); % eventually given versus [RG(1:Mp+1),x/c(1:Np)]
    
    for i = 1:Mp+1
             fof0(i,:) =      fof0_temp; 
        dfof0dxoc(i,:) = dfof0dxoc_temp;
             tot0(i,:) =      tot0_temp;
    end

end
% -------------------------------------------------------------------------



% ------------------ Find meanline and thickness profiles (at x1 positions)
% f      = camber               at x1 positions
% dfdx   = slope of camber line at x1 positions
% t      = thickness            at x1 positions
t    = zeros(Mp+1,Np);
f    = zeros(Mp+1,Np);
dfdx = zeros(Mp+1,Np);

for i = 1:Mp+1                 % for each radial section along the span
       f(i,:) =  fof0(i,:)    *f0oc(i)*c(i);
    dfdx(i,:) = dfof0dxoc(i,:)*f0oc(i);

       t(i,:) =  tot0(i,:) * t0(i);
end      
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% ------------------------------------- Find 2D unroatated section profiles
% x2D  [m], x   position in 2D space on upper (x2D_u) and lower (x2D_l) foil surfaces
% y2D  [m], y   position in 2D space on upper (x2D_u) and lower (x2D_l) foil surfaces
x2D_u = zeros(Mp+1,Np);     x2D_l = zeros(Mp+1,Np);
y2D_u = zeros(Mp+1,Np);     y2D_l = zeros(Mp+1,Np);

for i = 1:Mp+1                           % for each section along the span
    for j = 1:Np                         % for each point   along the chord
        x2D_u(i,j) = x1(i,j) + (t(i,j)/2)*sin(atan(dfdx(i,j))); % 2D upper surface x
        x2D_l(i,j) = x1(i,j) - (t(i,j)/2)*sin(atan(dfdx(i,j))); % 2D lower surface x
        y2D_u(i,j) =  f(i,j) + (t(i,j)/2)*cos(atan(dfdx(i,j))); % 2D upper surface y
        y2D_l(i,j) =  f(i,j) - (t(i,j)/2)*cos(atan(dfdx(i,j))); % 2D lower surface y
    end
end



% ----------------------------------------- Put all the numbers in one list
% % Nose -> suctioin side -> tail -> pressure side -> nose
% x2D(:,   1:Np   ) = x2D_u(:,1:Np);     % The first Np values are the upper surface (suction side),
% x2D(:,1+Np:Np+Np) = x2D_l(:,Np:-1:1);  % and the second Np values are the lower surface (pressure side).
% y2D(:,   1:Np   ) = y2D_u(:,1:Np);
% y2D(:,1+Np:Np+Np) = y2D_l(:,Np:-1:1);
  
% % j = 1          == tail
% % j = 1:Np       == suction side
% % j = Np         == nose
% % j = Np + 1     == nose
% % j = Np+ 1:2*Np == pressure side
% % j = 2*Np       == tail
% % Tail -> suctioin side -> nose, nose -> pressure side -> tail
y2D = zeros(Mp+1,2*Np);     
y2D = zeros(Mp+1,2*Np);
x2D(:,   1:Np   ) = x2D_u(:,Np:-1:1);   % The first Np values are the upper surface (suction side),
x2D(:,Np+1:Np+Np) = x2D_l(:,1:Np);      % and the second Np values are the lower surface (pressure side).
y2D(:,   1:Np   ) = y2D_u(:,Np:-1:1);
y2D(:,Np+1:Np+Np) = y2D_l(:,1:Np);



% --------------------------------------- Find 2D roatated section profiles
% x2Dr [m], x position in 2D space after rotation for pitch angle
% y2Dr [m], y position in 2D space after rotation for pitch angle
x2Dr = zeros(Mp+1,2*Np);
y2Dr = zeros(Mp+1,2*Np);
% for i = 1:Mp        % for each section along the span
for i = 1:Mp+1        % for each section along the span
    x2Dr(i,:) = x2D(i,:)*cosd(theta(i)) - y2D(i,:)*sind(theta(i)); % rotated 2D upper and lower surface x
    y2Dr(i,:) = x2D(i,:)*sind(theta(i)) + y2D(i,:)*cosd(theta(i)); % rotated 2D upper and lower surface y
end

% --------------------------- Invoke skew and rake, and find 3D coordinates
% X3D [m], X position in 3D space (corresponds to y position in 2D space)
% Y2D [m], Y position in 3D space
% Z3D [m], Z position in 3D space
X3D = zeros(Mp+1,2*Np);
Y3D = zeros(Mp+1,2*Np);
Z3D = zeros(Mp+1,2*Np);

for i = 1:Mp+1        % for each section along the span
    for j = 1:2*Np    % for each point   along the upper and lower surfaces
            X3D(i,j) = - rake(i) - r(i)*(pi*skew(i)/180)*tand(theta(i)) + y2Dr(i,j);
            
            Y3D(i,j) = r(i)*sind( skew(i) - (180/pi)*x2Dr(i,j)/r(i) );
            Z3D(i,j) = r(i)*cosd( skew(i) - (180/pi)*x2Dr(i,j)/r(i) );
    end
end

% % --------------------- If left-hand screw, then mirror the Y3D coordinates
% LeftHand_flag = 0;  % 1 == left-handed propeller, 0 == right-handed propeller
% 
% if LeftHand_flag == 1
%     Y3D = -Y3D;
% end
% % -------------------------------------------------------------------------



% Output variables:
xl = x2D_l;
xu = x2D_u;
yl = y2D_l;
yu = y2D_u;
