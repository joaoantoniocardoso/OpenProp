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

function [Ixc, Iyc, Ixyc, Area, Xbar, Ybar, xl, yl, xu, yu] = Stress_Moment_of_Inertia(xl,xu,yl,yu)

% xu(i,:) corresponds to section at RV(i)
[Mp1,Np] = size(xu);
Mp = Mp1-1;

for m=1:Mp1
    yshift = abs(min(yl(m,:)));     %Distance to shift all y points so that all are positive
    yu(m,:) = yu(m,:) + yshift;     %Shift of upper surface y points
    yl(m,:) = yl(m,:) + yshift;     %Shift of lower surface y points
    
    xshift = abs(min(min(xu(m,:)),min(xl(m,:))));   %Distance to shift all x points so that all are positive
    xu(m,:) = xu(m,:) + xshift;                     %Shift of upper surface x points
    xl(m,:) = xl(m,:) + xshift;                     %Shift of lower surface x points
    
end
dxu = abs(diff(xu,1,2));
dxl = abs(diff(xl,1,2));
dyu = diff(yu,1,2);
dyl = diff(yl,1,2);

Ybar = zeros(1,Mp1);
Xbar = zeros(1,Mp1);
Ixc  = zeros(1,Mp1);
Iyc  = zeros(1,Mp1);
Area = zeros(1,Mp1);
Ixyc = zeros(1,Mp1);


for m=1:Mp1
    hru  = zeros(1,Np-1);
    hrl  = zeros(1,Np-1);
    htu  = zeros(1,Np-1);
    htl  = zeros(1,Np-1);
    xctu = zeros(1,Np-1);
    xctl = zeros(1,Np-1);
    
    for n=1:(Np-1)
        
        hru(n) = min(yu(m,n),yu(m,n+1));   %Height of upper surface elemental rectangle
        htu(n) = max(yu(m,n),yu(m,n+1));   %Height of upper surface elemental trapezoid
        hrl(n) = min(yl(m,n),yl(m,n+1));   %Height of lower surface elemental rectangle
        htl(n) = max(yl(m,n),yl(m,n+1));   %Height of lower surface elemental trapezoid
        
        if dyu(m,n)<0
            xctu(n) = xu(m,n) + 2*dxu(m,n)/3;      %Distance from y-axis to upper surface elemental triangle centroid
        else
            xctu(n) = xu(m,n) + dxu(m,n)/3;        %Note: Value depends on whether left or right side of triangle is higher
        end
        
        if dyl(m,n)>0
            xctl(n) = xl(m,n) + 2*dxl(m,n)/3;      %Distance from y-axis to lower surface elemental triangle centroid
        else
            xctl(n) = xl(m,n) + dxl(m,n)/3;        %Note: Value depends on whether left or right side of triangle is higher
        end
        
    end
    xcru = xu(m,1:(Np-1))+dxu(m,:)/2; %Distance from y-axis to upper surface elemental rectangle
    xcrl = xl(m,1:(Np-1))+dxl(m,:)/2; %Distance from y-axis to lower surface elemental rectangle
    
    aru = dxu(m,:).*hru;            %Elemental upper surface rectangle area
    atu = dxu(m,:).*(htu-hru)/2;    %Elemental upper surface triangle area
    arl = dxl(m,:).*hrl;            %Elemental lower surface rectangle area
    atl = dxl(m,:).*(htl-hrl)/2;    %Elemental lower surface triangle area
    
    ycru = hru/2;                   %Distance from x-axis to upper surface elemental rectangle centroid
    ycrl = hrl/2;                   %Distance from x-axis to lower surface elemental rectangle centroid
    yctu = hru+(htu-hru)/3;         %Distance from x-axis to upper surface elemental triangle centroid
    yctl = hrl+(htl-hrl)/3;         %Distance from x-axis to lower surface elemental triangle centroid
    
    
    Mxsu = sum(ycru.*aru + yctu.*atu);  %1st moment of upper surface about x axis
    Mxsl = sum(ycrl.*arl + yctl.*atl);  %1st moment of lower surface about x axis    
    Mxs = Mxsu - Mxsl;
    
    Mysu = sum(xcru.*aru + xctu.*atu);  %1st moment of upper surface about y axis
    Mysl = sum(xcrl.*arl + xctl.*atl);  %1st moment of lower surface about y axis    
    Mys = Mysu - Mysl;
    
    Au = sum(aru + atu);               %Area of upper surface (x axis to upper surface)
    Al = sum(arl + atl);               %Area of lower surface (x axis to lower surface)
    Area(m) = Au - Al;

    Ybar(m) = Mxs/Area(m);         %Distance to centroid from x-axis
    Xbar(m) = Mys/Area(m);         %Distance to centroid from y-axis
    
%     figure(m)
%     plot(xu(m,:),yu(m,:),xl(m,:),yl(m,:),'b')
%     line([min(xu(m,:)),max(xu(m,:))],[Ybar(m),Ybar(m)],'Color','r','LineWidth',2,'LineStyle','--')
%     line([Xbar(m),Xbar(m)],[min(yl(m,:)),max(yu(m,:))],'Color','r','LineWidth',2,'LineStyle','--')
%     axis equal
%     grid on
        
    ixru = dxu(m,:).*hru.^3/3;
    ixyru = aru.*ycru.*xcru;
    ixtu = dxu(m,:).*(htu-hru).^3/36 + atu.*yctu.^2;
    ixytu = atu.*yctu.*xctu;
    
    ixrl = dxl(m,:).*hrl.^3/3;
    ixyrl = arl.*ycrl.*xcrl;
    ixtl = dxl(m,:).*(htl-hrl).^3/36 + atl.*yctl.^2;
    ixytl = atl.*yctl.*xctl;
    
    iyru = hru.*dxu(m,:).^3/12 + aru.*xcru.^2;
    iytu = (htu - hru).*dxu(m,:).^3/36 + atu.*xctu.^2;
    iyrl = hrl.*dxl(m,:).^3/12 + arl.*xcrl.^2;
    iytl = (htl - hrl).*dxl(m,:).^3/36 + atl.*xctl.^2;
    
    Ix = sum(ixru + ixtu) - sum(ixrl + ixtl);
    Iy = sum(iyru + iytu) - sum(iyrl + iytl);
    Ixy = sum(ixyru + ixytu) - sum(ixyrl + ixytl);

    Ixc(m) = Ix - Area(m)*Ybar(m)^2;
    Iyc(m) = Iy - Area(m)*Xbar(m)^2;
    Ixyc(m) = Ixy - Area(m)*Xbar(m)*Ybar(m);
end

% Take mean values to get these versus RC, so xu(i,:) now corresponds to section RC(i)
Area = 0.5*(Area(1:Mp)+Area(2:Mp+1));
Ixc  = 0.5*( Ixc(1:Mp)+ Ixc(2:Mp+1));
Iyc  = 0.5*( Iyc(1:Mp)+ Iyc(2:Mp+1));
Ixyc = 0.5*(Ixyc(1:Mp)+Ixyc(2:Mp+1));
Xbar = 0.5*(Xbar(1:Mp)+Xbar(2:Mp+1));
Ybar = 0.5*(Ybar(1:Mp)+Ybar(2:Mp+1));

xu   = 0.5*(  xu(1:Mp,:)+xu(2:Mp+1,:));
xl   = 0.5*(  xl(1:Mp,:)+xl(2:Mp+1,:));
yu   = 0.5*(  yu(1:Mp,:)+yu(2:Mp+1,:));
yl   = 0.5*(  yl(1:Mp,:)+yl(2:Mp+1,:));



end