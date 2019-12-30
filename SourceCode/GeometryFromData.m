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
% =================================== Determine Propeller Geometry Function    
%
% This function determines the geometry of the propeller.  It outputs
% the geometry as a 2D image, 3D image, and Rhino CAD file.
%
% Reference: 
%   J.S. Carlton, "Marine Propellers & Propulsion", ch. 3, 1994.
%
%   Abbott, I. H., and Von Doenhoff, A. E.; Theory of Wing Sections. 
%   Dover, 1959. 
%
% -------------------------------------------------------------------------
% Input Variables:
%
%   filename            file name prefix for all output files
%   Date_string         time and date to print on reports
%   Make2Dplot_flag     flag for whether to make 2D geometry plot
%   Make3Dplot_flag     flag for whether to make 3D geometry plot
%   Make_Rhino_flag     flag for whetehr to make a Rhino output file
%   Meanline            flag for choice of meanline  form
%   Thickness           flag for choice of thickness form
%
%   XR          [ ],    input radii / propeller radius
%   f0oc0       [ ],    input camber    / chord at each radius
%   t0oc0       [ ],    input thickness / chord at each radius
%   skew0       [deg],  input skew              at each radius
%   rake0       [ ],    input rake / diameter   at each radius
%
%   RC          [ ],    control point radii / propeller radius
%   CL          [ ],    section lift coefficients
%   BetaC      [deg],  Beta  at the control points
%   BetaIC     [deg],  BetaI at the control points
%   alphaI      [deg],  ideal angle of attack
%
%   D           [m],    propeller diameter
%   Z           [ ],    number of blades
%   Dhub        [m],    hub diameter
%   Rhub        [m],    hub radius
%
%   CoD         [ ],    chord / diameter at each control point radius
%   R           [m],    propeller radius
%   Mp          [ ],    number of radial 2D cross-sections
%   Np          [ ],    number of points in each 2D section
%
% Output Variables:
%
% The function has graphical and file outputs, in addition to the geometry 
% data structure.
%
% -------------------------------------------------------------------------

function [geometry] = GeometryFromData(pt)
%%
if     isfield(pt,'ginput'), i = pt.ginput; 
elseif isfield(pt,'gi'),     i = pt.gi; 
else disp('ERROR: need pt.ginput.'), return, end
    
% -------------------------------------------------------- Unpack variables
if isfield(pt,'date'), Date_string = pt.date; else  Date_string = date; end

if     isfield(pt      ,'filename'),  filename =       pt.filename;
elseif isfield(pt      ,'name'),      filename =       pt.name;
elseif isfield(pt.input,'filename'),  filename = pt.input.filename; 
else                                  filename = 'OpenProp';
end


if isfield(i, 'Hub_flag'),   Hub_flag =  i.Hub_flag;  else    Hub_flag = 1;  end  % 0 == no hub, 1 == hub
if isfield(i,'Duct_flag'),  Duct_flag = i.Duct_flag;  else   Duct_flag = 0;  end  % 0 == no duct, 1 == duct

% 0 == lifting line at mid-chord, 1 == lifting line at quarter-chord (i.e. skew offset == 0.25*c/r)
if isfield(i,'QuarterChord_flag'), QuarterChord_flag = i.QuarterChord_flag; else   QuarterChord_flag = 0; end
                          
if isfield(i, 'Meanline'), Meanline  = i.Meanline;   else  Meanline  = 'NACA a=0.8 (modified)';   end
if isfield(i,'Thickness'), Thickness = i.Thickness;  else  Thickness = 'NACA 65A010';             end


Z = i.Z;

if isfield(i,    'R'), R     = i.R;     else  R     = 1;    end % m
if isfield(i, 'Rhub'), Rhub  = i.Rhub;  else  Rhub  = 0.2;  end % m
if isfield(i,'Rcirc'), Rcirc = i.Rcirc; else  Rcirc = Rhub; end % m
if isfield(i,'Rroot'), Rroot = i.Rroot; else  Rroot = Rhub; end % m
if isfield(i,'Rduct'), Rduct = i.Rduct; else  Rduct = R;    end % m

if isfield(i,  'Mp'), Mp   = i.Mp;    else  Mp   = 20;  end
if isfield(i,  'Np'), Np   = i.Np;    else  Np   = 20;  end



XR  = i.XR; 
X1  =  ones(size(XR));
X0  = zeros(size(XR));

XCoD  = i.XCoD;
    
if     isfield(i,'XPoD'),    XPoD = i.XPoD;  
elseif isfield(i,'Xtheta'),  XPoD = tand(i.Xtheta).*pi.*XR;  % Xtheta [deg] design pitch angle
elseif isfield(i,'theta0'),  XPoD = tand(i.theta0).*pi.*XR;
else   disp('ERROR: need pitch or pitch angle.')
end

    
if     isfield(i,'Xt0oD'),    Xt0oD = i.Xt0oD;
elseif isfield(i,'Xt0oc'),    Xt0oD = i.Xt0oc .* XCoD;
elseif isfield(i,'t0oc0'),    Xt0oD = i.t0oc0 .* XCoD;
else   disp('ERROR: need thickness.')
end

    
if     isfield(i,'Xf0oc'),  Xf0oc = i.Xf0oc;          
elseif isfield(i,'f0oc0'),  Xf0oc = i.f0oc0; 
else   disp('ERROR: need camber.')
end 


if isfield(i,'skew0'), skew0  = i.skew0; else skew0 = X0; end % [deg]
if isfield(i,'rake0'), rake0  = i.rake0; else rake0 = X0; end

  

% ---------
if Rcirc < Rhub,  disp('ERROR: Rcirc must be >= Rhub.' ), end
if Rroot < Rcirc, disp('ERROR: Rroot must be >= Rcirc.'), end

if isfield(i,'t0circ'), t0circ = i.t0circ; else  t0circ = XCoD(1); end % m
% ---------
    
D        = 2*R;    % [m]
Dhub     = 2*Rhub; % [m]

Rhub_oR  = Rhub/R;    % [ ], hub                  radius / propeller radius
Rroot_oR = Rroot/R;   % [ ], blade root           radius / propeller radius
Rcirc_oR = Rcirc/R;   % [ ], circular section max radius / propeller radius
Rduct_oR = Rduct/R;   % [ ], duct                 radius / propeller radius
RoR      = 1;         % [ ], propeller            radius / propeller radius

% Input flags      
if isfield(i,'Make2Dplot_flag'), 
        Make2Dplot_flag = i.Make2Dplot_flag; 
else    Make2Dplot_flag = 1;  
end   
if isfield(i,'Make3Dplot_flag'), 
        Make3Dplot_flag = i.Make3Dplot_flag; 
else    Make3Dplot_flag = 1;  
end   
if isfield(i,'Make_Rhino_flag'), 
        Make_Rhino_flag = i.Make_Rhino_flag; 
else    Make_Rhino_flag = 0;  
end   
if isfield(i,'Make_SWrks_flag'), 
        Make_SWrks_flag = i.Make_SWrks_flag; 
else    Make_SWrks_flag = 0;  
end

% ---------------------------------------------------------- Model Diameter
% This is the diameter of the output SolidWorks/Rhino geometry
if isfield(i,'Dm'), Dm = i.Dm;  else Dm = D;    end % m
Rm = Dm/2;


% -------------------------------------------------------------------------
% % ---------------- Interpolate input geometry at selected radial sections
% % % RG = [0.9*Rhub_oR,RV(2:end-1),1];
% % 
% % % Interpolate input geometry at sections with cosine spacing along the span
% % if Rhub == Rroot
% %     if Hub_flag == 1 
% %         RG = 0.9*Rhub_oR + (1-0.9*Rhub_oR)*(sin((0:Mp)*pi/(2*Mp)));  % [0.9*Rhub_oR : 1]
% %     else
% %         RG = Rhub_oR + (1-Rhub_oR)*(sin((0:Mp)*pi/(2*Mp)));  % [0.9*Rhub_oR : 1]
% %         
% %         % % This would be useful for a sraight wing...
% %         % RG  = 0.5*( Rhub_oR+1) + 0.5*(1- Rhub_oR)*(sin((0:Mp)*pi/Mp - pi/2));  % [Rhub_oR : 1]
% %     end
% % else % Rhub < Rroot
% %         RG  = 0.5*(Rroot_oR+1) + 0.5*(1-Rroot_oR)*(sin((0:Mp)*pi/Mp - pi/2));  % [Rroot_oR : 1]
% %         
% %         DRG1 = RG(2)-RG(1);
% %         DRG2 = RG(3)-RG(2);
% %         
% %         RGE = linspace(Rcirc_oR,Rroot_oR,floor((Rroot_oR-Rcirc_oR)/DRG1));  % radii of elliptical sections
% % 
% %         RGC = linspace(Rhub_oR,Rcirc_oR, floor((Rcirc_oR -Rhub_oR)/DRG2)); % [Rhub_oR : Rcirc_oR] radii of circular sections
% % %         RGC = RGC(1:end-1);
% % end
% % 
% % % % % % disp('12/22/10 NOTE HACK: RG == [Rroot_oR : 1]')
% % % % % RG  = 0.5*(Rroot_oR+1) + 0.5*(1-Rroot_oR)*(sin((0:Mp)*pi/Mp - pi/2));  % [Rroot_oR : 1]

% disp('2/10/11 NOTE CORRECTION TO ABOVE HACK: RG == [Rroot_oR : 1]')
RG = Rhub_oR + (1-Rhub_oR)*(sin((0:Mp)*pi/(2*Mp))); 


% disp('success')


% -------- Thickness, rake, chord, and radius scaled for the model diameter
CoD  =  InterpolateChord(XR,XCoD,RG);
t0   = pchip(XR,Xt0oD,RG)*Dm;
f0oc = pchip(XR,Xf0oc,RG);
skew = pchip(XR,skew0,RG);        % [deg], angular translation along mid-chord helix
rake = pchip(XR,rake0,RG)*Dm;     % [m],   translation along propeller axis (3D X-axis)
PoD  = pchip(XR,XPoD ,RG);

theta    = atand(PoD./(pi*RG));   % Nose-tail pitch angle, [deg]

t0circ = t0(1);


c = CoD.*Dm;                          % section chord at the RG sections [m]
r = RG .*Rm;                          % radius of        the RG sections [m]
% -------------------------------------------------------------------------
% disp('success')
%
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
  % % Even spacing along the chord   
  % x0(j) =                 (j-1)/(Np-1);  % [0 : 1]
    
    % Cosine spacing along the chord
    x0(j) = 0.5*(1-cos(pi*(j-1)/(Np-1)));  % [0 : 1]
end

% QuarterChord_flag = 1;

for i = 1:Mp+1                     % for each radial section along the span    
    if QuarterChord_flag == 1
        x1(i,:) = c(i)/4 - c(i)*x0;    % lifting line at quarter-chord
    else        
        x1(i,:) = c(i)/2 - c(i)*x0;    % lifting line at mid-chord
    end
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ---------------------- Find normalized 2D foil geometry (at x0 positions)
%   f0octilde = zeros(Mp+1,1);  % can either be scalars or given versus RG(1:Mp+1)
%    CLItilde = zeros(Mp+1,1);  % can either be scalars or given versus RG(1:Mp+1)
% alphaItilde = zeros(Mp+1,1);  % can either be scalars or given versus RG(1:Mp+1)
if iscell(Meanline)  % Assume Meanline is given versus XR
    
    fof0_temp        = zeros(length(XR),Np); % eventually given versus [RG(1:Mp+1),x/c(1:Np)]
    dfof0dxoc_temp   = zeros(length(XR),Np); % eventually given versus [RG(1:Mp+1),x/c(1:Np)]
    tot0_temp        = zeros(length(XR),Np); % eventually given versus [RG(1:Mp+1),x/c(1:Np)]    
    
    fof0        = zeros(Mp+1,Np); % eventually given versus [RG(1:Mp+1),x/c(1:Np)]
    dfof0dxoc   = zeros(Mp+1,Np); % eventually given versus [RG(1:Mp+1),x/c(1:Np)]
    tot0        = zeros(Mp+1,Np); % eventually given versus [RG(1:Mp+1),x/c(1:Np)]    
      
    
    for i = 1:length(XR)
        [junk1, junk2, junk3, fof0_temp(i,:), dfof0dxoc_temp(i,:), tot0_temp(i,:)] = GeometryFoil2D(Meanline{i},Thickness{i},x0);
    end    
      
    for j = 1:Np
             fof0(:,j) = pchip(XR,      fof0_temp(:,j), RG');
        dfof0dxoc(:,j) = pchip(XR, dfof0dxoc_temp(:,j), RG');
             tot0(:,j) = pchip(XR,      tot0_temp(:,j), RG');
    end

else
    [junk1, junk2, junk3, fof0_temp, dfof0dxoc_temp, tot0_temp] = GeometryFoil2D(Meanline,Thickness,x0);
    
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



% % -------------------------------------- Compute leading edge radius points
% phiLEC = atan(dfdxLE);
% NLE    = 3; % must be odd to capture leading edge point
% phiLEs = 3*pi/8;
% phiLE  = phiLEs:(pi-2*phiLEs)/(NLE-1):pi-phiLEs; 
% xLEC  = x1(:,1)' - rLE.*cos(phiLEC);
% yLEC  =            rLE.*sin(phiLEC);
% 
% xLE = zeros(Mp+1,NLE);
% yLE = zeros(Mp+1,NLE);
% 
% for i = 1:Mp+1                           % for each section along the span
%     xLE(i,:) = xLEC(i) + rLE(i)*sin(phiLE+phiLEC(i));
%     yLE(i,:) = yLEC(i) - rLE(i)*cos(phiLE+phiLEC(i));
% end


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


% % % Arrange points as follows:
% % %     Tail -> suctioin side -> leading edge (with radius and nose) -> pressure side -> tail
% % % j = 1                               == [1         point ] tail
% % % j = 1              : Np-1           == [Np-1      points] suction side (tail to point aft of leading edge radius)
% % % j = Np-1+1         : Np-1+(NLE-1)/2 == [(NLE-1)/2 points] suction side along leading edge radius
% % % j = Np+(NLE-1)/2                    == [1         point ] nose
% % % j = Np+(NLE-1)/2+1 :    Np-1+NLE    == [(NLE-1)/2 points] pressure side along leading edge radius
% % % j = Np-1+NLE+1     : 2*(Np-1)+NLE   == [Np-1      points] pressure side (point aft of leading edge radius to tail)
% % % j = 2*(Np-1)+NLE                    == [1         point ] tail
% % %
% % % j =            1   : Np+(NLE-1)/2   == suction  side (tail to nose)
% % % j = Np+(NLE-1)/2   : 2*(Np-1)+NLE   == pressure side (nose to tail)
% % 
% % x2D(:,1              :    Np-1     ) = x2D_u(:,Np:-1:2); % [Np-1 points] suction  side (tail to point aft of leading edge radius)
% % x2D(:,Np-1+1         :    Np-1 +NLE) =   xLE(:,NLE:-1:1);     % [NLE  points] leading edge radius  
% % x2D(:,Np-1+NLE+1     : 2*(Np-1)+NLE) = x2D_l(:,2:Np);    % [Np-1 points] pressure side (point aft of leading edge radius to tail)
% % 
% % y2D(:,1              :    Np-1     ) = y2D_u(:,Np:-1:2); % [Np-1 points] suction  side (tail to point aft of leading edge radius)
% % y2D(:,Np-1+1         :    Np-1 +NLE) =   yLE(:,NLE:-1:1);     % [NLE  points] leading edge radius  
% % y2D(:,Np-1+NLE+1     : 2*(Np-1)+NLE) = y2D_l(:,2:Np);    % [Np-1 points] pressure side (point aft of leading edge radius to tail)



% %--------------------------------------- plot unrotated blade
%  Fig2_S = figure('units','normalized','position',[0.31 .06 .4 .3],'name',...
%         'Blade Image','numbertitle','off');
%     style=['r' 'g' 'b' 'm' 'k'];
%     str_prefix = {'r/R = '};
%     flag=1;
%     for i = 1:ceil(Mp/5):Mp     % for five radial sections from root to tip
%         plot(x2D(i,:),y2D(i,:),style(flag));
%         
%         for j = 1:Np
%             plot([x2D_l(i,j),x2D_u(i,j)],[y2D_l(i,j),y2D_u(i,j)],style(flag));
%         end
%         
%         str_legend(flag)=strcat(str_prefix,num2str(RC(i)));
%         hold on;
%         flag = flag+1;
%     end
%     legend(str_legend,'location','northwest');
%     axis equal;     grid on;
%     title('2D Blade Image');  xlabel('X (2D) [m]');  ylabel('Y (2D) [m]');
% %---------------------------------------






% --------------------------------------- Find 2D roatated section profiles
% x2Dr [m], x position in 2D space after rotation for pitch angle
% y2Dr [m], y position in 2D space after rotation for pitch angle
% x2Dr = zeros(Mp+1,2*(Np-1)+NLE);
% y2Dr = zeros(Mp+1,2*(Np-1)+NLE);
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
% X3D = zeros(Mp+1,2*(Np-1)+NLE,Z);
% Y3D = zeros(Mp+1,2*(Np-1)+NLE,Z);
% Z3D = zeros(Mp+1,2*(Np-1)+NLE,Z);
X3D = zeros(Mp+1,2*Np,Z);
Y3D = zeros(Mp+1,2*Np,Z);
Z3D = zeros(Mp+1,2*Np,Z);
theta_Z  = 0:360/Z:360;                   % angle between blades [deg]

% for i = 1:Mp        % for each section along the span
for i = 1:Mp+1        % for each section along the span
%     for j = 1:2*(Np-1)+NLE    % for each point   along the upper and lower surfaces
    for j = 1:2*Np    % for each point   along the upper and lower surfaces
        for k = 1:Z   % for each blade
            X3D(i,j,k) = - rake(i) - r(i)*(pi*skew(i)/180)*tand(theta(i)) + y2Dr(i,j);
            
            Y3D(i,j,k) = r(i)*sind(skew(i) - (180/pi)*x2Dr(i,j)/r(i) - theta_Z(k));
            Z3D(i,j,k) = r(i)*cosd(skew(i) - (180/pi)*x2Dr(i,j)/r(i) - theta_Z(k));
        end
    end
end

% % ---- Find axial length of blade
% max(max(max(X3D))) - min(min(min(X3D)))



% --------------------- If left-hand screw, then mirror the Y3D coordinates
LeftHand_flag = 0;  % 1 == left-handed propeller, 0 == right-handed propeller

if LeftHand_flag == 1
    Y3D = -Y3D;
end
% -------------------------------------------------------------------------

% save geometry X3D Y3D Z3D

% --------------------------------------------- Compute Expanded Area Ratio
EAR  = (2*Z/pi)*trapz(linspace(Rhub_oR,R/R,100),interp1(RG,CoD,linspace(Rhub_oR,R/R,100),'spline','extrap'));           

% ---------------------------------------- Compute Blade Thickness Fraction
BTF = interp1(XR,Xt0oD,0,'linear','extrap');  % BTF = t0oD at the prop centerline


% =========================================================================
% ============================= Pack up geometry data at the control points
t0oc    = t0 ./ c;                   % [ ],   t0/c

geometry.Meanline  = Meanline;
geometry.Thickness = Thickness;

geometry.Z         = Z;
geometry.D         = D;                          % [m]
geometry.Dhub      = Dhub;                       % [m]

geometry.EAR       = EAR;
geometry.BTF       = BTF;

geometry.RG        = RG;                         % r/R
geometry.CoD       = CoD;  
geometry.t0oc      = t0oc;           
geometry.f0oc      = f0oc;      
geometry.theta     = theta;
geometry.PoD       = PoD; 
geometry.skew      = skew;                       % [deg]
geometry.rake      = rake;                       % [m] 

% =========================================================================
% =========================================================================
% =========================================== Create plots and text outputs
%
% ----------------------------------------- Create 2D Propeller Blade Image
Make2Dplot_flag = 0;

if Make2Dplot_flag

    Make_2D_Blade_Image(RG,x2Dr,y2Dr)

%     filename_2D = strcat(filename,'_2D_Blade_Image');
%     saveas(gcf,filename_2D,'jpg')
end



% ----------------------------------------------- Create 3D Propeller Image
if Make3Dplot_flag
    Js = 1; % temporary hack...
    BetaIG = theta;
    RV = RG;
    TANBIV = tand(theta);
    Duct_flag = 0;
    Cduct_oR = 1;
    Xduct_oR = 0;
    GUI_flag = 0;
    Plots      = [];
    PlotPanels = [];

    
    Make_3D_Blade_Image(X3D,Y3D,Z3D,Duct_flag,Rduct_oR,Cduct_oR,Xduct_oR, Hub_flag,Rhub_oR,  Js,BetaIG,theta,TANBIV,RV,Rm, GUI_flag,Plots,PlotPanels)

    if Duct_flag == 1

        [xd,rd,Xd,Yd,Zd] = Duct_Plot_120329(RC,TANBIC,G,VARING,RV,Z,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR,Cduct_oR,Xduct_oR,Gd,Meanline_d,Thickness_d,x0,t0oc_duct,R,Mp,Np);

        surf(Xd,Yd,Zd)
    end    

%%
%     % ------------------------------------------------------ Save the image
%     view(2)
%     saveas(gcf,[filename,'_3D_Propeller_Image','1'],'jpg')
% 
%     view(3)
%     saveas(gcf,[filename,'_3D_Propeller_Image','2'],'jpg')

end                                              % (END IF Make3Dplot_flag)



%%
% --------------------------------------------------- Make SolidWorks files
Make_SWrks_flag = 1;

if Make_SWrks_flag
    % Make SolidWorks.txt files, with coordinates for a single blade
    
    filename_SolidWorks = strcat(filename,'_SolidWorks_v18.txt');    
    Export_SolidWorks_v18(filename_SolidWorks,Np,Mp,Z,X3D,Y3D,Z3D);
    
%     filename_SolidWorks = strcat(filename,'_SolidWorks_v14.txt');
%     Export_SolidWorks_v14(filename_SolidWorks,Np,Mp,Z,X3D,Y3D,Z3D);
end                                              % (END IF Make_SWrks_flag)
%%

% -------------------------------------------------------- Make Rhino files 
if Make_Rhino_flag
    % Make _RhinoBlade.txt, with coordinates for a single blade and
    % commands to make Z blades
    
    filename_Rhino = strcat(filename,'_RhinoProp.txt');
        
    Export_Rhino_v1(filename_Rhino,rake,R,XR,skew0,Mp,Np,X3D,Y3D,Z3D,Z,Rhub);        

end                                              % (END IF Make_Rhino_flag)
%%
% ---------------------------------------------- Make OpenProp_Geometry.txt
filename_geometry = strcat(filename,'_Geometry.txt');
fid = fopen(filename_geometry,'w');

fprintf(fid,'\t\t\t\t\t %s \n\n',filename_geometry);
fprintf(fid,'\t\t\t\t\t Propeller Geometry Table\n\n');
fprintf(fid,'Date and time: %s\n\n',Date_string);

fprintf(fid,'Number of Blades \t = %.0f\n',        Z);
fprintf(fid,'Propeller Diameter \t = %.4f m\n',    D);
fprintf(fid,'Propeller Hub Diameter \t = %.4f m\n',Dhub);

    if iscell(Meanline)
        
            fprintf(fid,'Meanline  Type: ')
            
        for j = 1:length(Meanline)
            
            fprintf(fid,[Meanline{j},', ']);
        end
            fprintf(fid,'\n');
            
            
            fprintf(fid,'Thickness  Type: ')
            
        for j = 1:length(Thickness)
            
            fprintf(fid,[Thickness{j},', ']);
        end
            fprintf(fid,'\n');        
    else
        fprintf(fid,['Meanline  Type: ',Meanline,'\n']);
        fprintf(fid,['Thickness Type: ',Thickness,'\n']);
    end

fprintf(fid,' \n');
fprintf(fid,' \n');

fprintf(fid,' r/R\t P/D\t Skew\t Xs/D\t  c/D\t  f0/c\t  t0/c\n');
for i = 1:Mp+1
    fprintf(fid, '%.4f\t %.4f\t %.4f\t %.4f\t %.4f\t %.4f\t %.4f\n'...
    ,RG(i),PoD(i),skew(i),rake(i)/D,CoD(i),f0oc(i),t0oc(i));
end

fprintf(fid,' \n');
fprintf(fid,'\nr/R \t [ ], radial position of control points / propeller radius.\n');
fprintf(fid,'P/D \t [ ], section pitch / diameter.\n');
fprintf(fid,'c/D \t [ ], section chord-length / diameter.\n');
fprintf(fid,'fo/C \t [ ], section camber / section chord-length.\n');
fprintf(fid,'to/C \t [ ], section thickness / section chord-length.\n');

fclose(fid);

% =============================== END Determine Propeller Geometry Function
% =========================================================================