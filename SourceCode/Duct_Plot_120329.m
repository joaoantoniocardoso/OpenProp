
% ------------------------------------------------------------------------- 
function [x,r,X,Y,Z] = Duct_Plot_120329(RC,TANBIC,G,VARING,RV,Z,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR,Cduct_oR,Xduct_oR,Gd,Meanline_d,Thickness_d,x0,t0oc_duct,R,Mp,Np)
% -------------------------------------------------------------------------


Plot_flag = 0;

Rduct = Rduct_oR*R;
Cduct = Cduct_oR*R;
Xduct = Xduct_oR*R;

t0oc_duct = 0.1;

% Gd = 0.1095/ 2;


% ------------------ Find propeller influence at the duct quarter chord
disp(' '), disp('Computing rotor-duct interaction...be patient...'), disp(' '),
[DAHIFq_times_TANBIC,junk,DRHIFq_times_TANBIC] = Horseshoe_intr_110830(Xduct_oR-Cduct_oR/4,Rduct_oR ,RC,ones(size(RC)),RV,Z,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR); 


% ------------------------------------------ Flow at duct quarter chord
for m = 1:Mp;                               
    DAHIFq(:,m) = DAHIFq_times_TANBIC(:,m) / TANBIC(m);
    DRHIFq(:,m) = DRHIFq_times_TANBIC(:,m) / TANBIC(m);
end     

UARINGq = (DAHIFq*G)';  
URRINGq = (DRHIFq*G)'; 


VSRINGq = sqrt((VARING+UARINGq)^2+(URRINGq)^2);

CLd     = 4*pi*Gd / (VSRINGq*Cduct_oR);

BetaIDq  = atand( -URRINGq/(VARING+UARINGq) );      % [deg] inflow angle

% ----------------------------------------------- Duct pitch and camber

Meanline_d = 'parabolic';

% ---------------------- Find normalized 2D foil geometry (at x0 positions)  
[f0octilde, CLItilde, alphaItilde, fof0, dfof0dxoc, tot0] = GeometryFoil2D(Meanline_d,Thickness_d,x0);


% ----- Scale camber ratio and ideal angle of attack by 2D lift coefficient
alphaI = alphaItilde * CLd / CLItilde;        % [deg], ideal angle of attack
  f0oc =   f0octilde * CLd / CLItilde;        % f0/c, scaled for CL at RG    

 
     

% ------------------ Find meanline and thickness profiles (at x1 positions)
% f      = camber               at x1 positions
% dfdx   = slope of camber line at x1 positions
% t      = thickness            at x1 positions
%
% Note: x0 == 0 at leading edge, and x0 == 1 at trailing edge
% Note: X3D is positive upstream of propeller plane
% Note: Duct midchord (x0 == 0.5) is at X3D == -Xduct

x1 = -Xduct + Cduct * (0.5-x0) ;  


       f =  fof0     * f0oc     *Cduct;
    dfdx = dfof0dxoc * f0oc;

       t =  tot0     * t0oc_duct*Cduct;
      
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% ------------------------------------- Find 2D unroatated section profiles
% x2D  [m], x   position in 2D space on upper (x2D_u) and lower (x2D_l) foil surfaces
% y2D  [m], y   position in 2D space on upper (x2D_u) and lower (x2D_l) foil surfaces
%
x2D_u = x1 + (t/2).*sin(atan(dfdx)); % 2D upper surface x
x2D_l = x1 - (t/2).*sin(atan(dfdx)); % 2D lower surface x
y2D_u =  f + (t/2).*cos(atan(dfdx)); % 2D upper surface y
y2D_l =  f - (t/2).*cos(atan(dfdx)); % 2D lower surface y

% ----------------------------------------- Put all the numbers in one list 
% % j = 1          == tail
% % j = 1:Np       == suction side
% % j = Np         == nose
% % j = Np + 1     == nose
% % j = Np+ 1:2*Np == pressure side
% % j = 2*Np       == tail
% % Tail -> suctioin side -> nose, nose -> pressure side -> tail
clear x2D y2D

x2D(   1:Np   ) = x2D_u(Np:-1:1);   % The first Np values are the upper surface (suction side),
x2D(Np+1:Np+Np) = x2D_l(1:Np);      % and the second Np values are the lower surface (pressure side).
y2D(   1:Np   ) = y2D_u(Np:-1:1);
y2D(Np+1:Np+Np) = y2D_l(1:Np);
      
      

if Plot_flag == 1
%--------------------------------------- plot unrotated blade
 fig;
        plot(x2D,y2D,'g');
        
        for j = 1:Np
            plot([x2D_l(j),x2D_u(j)],[y2D_l(j),y2D_u(j)],'k');
        end
        
       axis equal;     grid on;
    xlabel('X (2D) [m]');  ylabel('Y (2D) [m]');
%---------------------------------------
end


% ---------------------------------------------- Find pitch angle and pitch
theta    = -(BetaIDq + alphaI);               % Nose-tail pitch angle, [deg]


% --------------------------------------- Find 2D roatated section profiles
x2Dr = x2D*cosd(theta) - y2D*sind(theta); % [m] rotated 2D upper and lower surface x
y2Dr = x2D*sind(theta) + y2D*cosd(theta); % [m] rotated 2D upper and lower surface y


if Plot_flag == 1

       plot(x2Dr,y2Dr,'b');
       plot(x2Dr(1:Np),y2Dr(1:Np),'r');  % inside surface
        
        for j = 2:Np-1
            plot([x2Dr(j),x2Dr(2*Np+1-j)],[y2Dr(j),y2Dr(2*Np+1-j)],'k');
        end
end    
     
% find y at x=0 for inside surface
y0 = pchip(x2Dr(1:Np),y2Dr(1:Np), 0);
    
if Plot_flag == 1
    plot(0,y0,'.b')

    plot(x2Dr,y2Dr - y0,'m')
end

%%    
% rotated section, placed at correct radius    
x = x2Dr;
r = Rduct + y0 - y2Dr;

if Plot_flag == 1
plot(x,r,'k')
end


Nx = length(x);
Nt = 50;
t  = linspace(0,2*pi,Nt);

X = zeros(Nx,Nt);
Y = zeros(Nx,Nt);
Z = zeros(Nx,Nt);

for n = 1:Nt  
    
    X(:,n) = x;
    Y(:,n) = r*cos(t(n));
    Z(:,n) = r*sin(t(n));
    
end


if Plot_flag == 1

fig;
    surf(X,Y,Z)
    axis equal
end
% ------------------------------------------------------------------------- 