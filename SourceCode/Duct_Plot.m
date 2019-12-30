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
% ===================================================== ductPlot Function

%Plots duct

%Variables:
% %   c [m]:            chordlength
% %   alpha [radians]:  angle of attack
% %   vrRad [m]:        vortex ring radius (duct radius to meanline)
% %   ductRef:          chordwise reference position on duct
% %   xDuct [m]:        global propeller x-coord of ductRef
% %   fo:               max camber (% of chordlength)
% %   to:               max thickness (% of chordlength)

% % Notes: X-axis positive in streamwise direction (i.e. downstream).

function [] = Duct_Plot(vrRad,c,fo,to,alpha,ductRef)


% % %Read in meanline f(x) and thickness t(x) distribution data
% %         %Read foil data (x, f/fo, t/to) from text file
% %         %Foil_data.txt contains parabolic meanline (f/fo)
% %         %and elliptical thickness(t/to) data
% %
% % [x_over_c,f_over_fo,t_over_to]=textread('foil_data.txt','%f%f%f',...
% %                                         'headerlines',3);
% % %x_over_c range is -c/2 to c/2 (this is converted to x=0 to x=c below)

% % x=x_over_c*c + c/2;     %dimensionalizes x with a range of 0 to c
% %                         %range of 0 to c is needed for cosine spacing

x = [0 .5 .75 1.25 2.5 5 7.5 10 15 20 25 30 35 40 45 50 ...
     55 60 65 70 75 80 85 90 95 100]./100;
x=x'*c;

% ------------------------- Use NACA a=0.8 meanline
    f_over_fo = [0 .287 .404 .616 1.077 1.841 2.483 3.043 3.985 4.748 ...
           5.367 5.863 6.248 6.528 6.709 6.79 6.77 6.644 6.405  ...
           6.037 5.514 4.771 3.683 2.435 1.163 0]./6.79;
    f_over_fo = f_over_fo';

% ------------------ Use NACA 65A010 thickness form
    t_over_to = [0 .765 .928 1.183 1.623 2.182 2.65 3.04 3.658 4.127 ...
              4.483 4.742 4.912 4.995 4.983 4.863 4.632 4.304     ...
              3.899 3.432 2.912 2.352 1.771 1.188 .604 .021]./4.995;
    t_over_to =  t_over_to';

f=fo*c*f_over_fo;         %camber distribution
t=to*c*t_over_to;         %thickness distribution

fpp=spline(x,f);        %spline camber data
% % fPpp=fnder(fpp);    %splines slope of fpp
tpp=spline(x,t);        %spline thickness data

%alt method not using FNDER (Spline Toolbox)
xl=length(x);
theta=zeros(xl,1);
for m=1:xl
    if m==xl
        theta(m)=theta(m-1);
    else
        theta(m)=atan((ppval(fpp,x(m+1))-ppval(fpp,x(m)))/(x(m+1)-x(m)));
    end
end

%Generate 2-D flat cross-section
% % theta   = atan(ppval(fPpp,x));
x_upper = x - t/2.*sin(theta);
y_upper = f + t/2.*cos(theta);
x_lower = x + t/2.*sin(theta);
y_lower = f - t/2.*cos(theta);

% % % %Plot 2-D section with 0 degrees angle of attack
% % % plot(x_lower,y_lower)
% % % hold on
% % % plot(x_upper,y_upper)
% % % xlabel('X-axis');ylabel('Y-axis');
% % % title('Duct section with 0 degrees angle of attack');
% % % axis equal
% % % figure

%Reposition section
var1=ductRef*c;                 %offset (x-dir), ductRef at x=0
var2=ppval(fpp,ductRef*c);      %offset (y-dir) for camber, ductRef at y=0
var3=0.5*ppval(tpp,ductRef*c);  %offset (y-dir) for thick, ductRef at y=0
                                %ductRef is point where blade and duct meet
x_upper = -(x_upper - var1);
y_upper =   y_upper - var2 + var3;
x_lower = -(x_lower - var1);
y_lower =   y_lower - var2 + var3;

%Rotate for angle of attack and place section at correct radius
x_upper_rot = x_upper*cos(-alpha) - y_upper*sin(-alpha);
y_upper_rot = x_upper*sin(-alpha) + y_upper*cos(-alpha) + vrRad;
x_lower_rot = x_lower*cos(-alpha) - y_lower*sin(-alpha);
y_lower_rot = x_lower*sin(-alpha) + y_lower*cos(-alpha) + vrRad;

% % % %Plot duct section rotated
% % % plot(x_lower_rot,y_lower_rot)
% % % hold on
% % % plot(x_upper_rot,y_upper_rot)
% % % xlabel('X-axis');ylabel('Y-axis');
% % % title(['Duct section (repositioned) with ',num2str(alpha*180/pi),...
% % %                                          ' degrees angle of attack']);
% % % axis equal
% % % figure

%Build all sections (upper and lower surfaces) for complete 3-D duct
z=zeros(length(x),1);
[thetaU,phiU,RU]=cart2sph(y_upper_rot,x_upper_rot,z);
[thetaL,phiL,RL]=cart2sph(y_lower_rot,x_lower_rot,z);

nds=50;                         %# of duct sections for plotting
for n=1:nds
%     phi=0+pi:1*pi/(nds-1):2.1*pi+pi;    %180 deg coverage for duct
    phi=0:2*pi/(nds-1):2.1*pi;    %360 deg coverage for duct
    [x_u_3D(n,:),y_u_3D(n,:),z_u_3D(n,:)]=sph2cart(phi(n),thetaU,RU);
    [x_l_3D(n,:),y_l_3D(n,:),z_l_3D(n,:)]=sph2cart(phi(n),thetaL,RL);

% % %     %Plot duct sections individually
% % %     plot3(z_u_3D(n,:),x_u_3D(n,:),y_u_3D(n,:))
% % %     hold on
% % %     plot3(z_l_3D(n,:),x_l_3D(n,:),y_l_3D(n,:))
end
% % % xlabel('X-axis');ylabel('Y-axis');zlabel('Z-axis')
% % % title(['Duct with ',num2str(alpha*180/pi),' degrees angle of attack'])
% % % axis equal
% % % figure

%Plot duct as 3-D surface
surfl(z_u_3D,x_u_3D,y_u_3D)
hold on
surfl(z_l_3D,x_l_3D,y_l_3D)
% % % xlabel('X-axis');ylabel('Y-axis');zlabel('Z-axis')
% % % title(['Duct with ',num2str(alpha*180/pi),' degrees angle of attack'])
% % % axis equal

%shading interp
%colormap(copper)
% ================================================= END ductPlot Function
% =========================================================================


