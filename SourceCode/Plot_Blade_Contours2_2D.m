% --------------------------------------------------- Plot_Blade_Contours.m 
% Created: Jerod Ketcham, 4/2/2010
% 
%
% This function plots s as a patch surface.  Data is formatted as:
%
% X3D,Y3D,Z3D are coordinates of each vertex.  
% s is the same size as X3D and contains the value at each vertex.
%
% X3D == [Mp1,Np]
% Y3D == [Mp1,Np]
% s   == [Mp1,Np]
%
% -------------------------------------------------------------------------

function [] = Plot_Blade_Contours2_2D(X3D,Y3D,s,plottitle)


[Mp1,Np1] = size(X3D);  % X3D is from the geometry.m module

Mp = Mp1-1;   
Np = Np1-1;

% -------------------------------------------------------------------------
% Concatenate matrices to create vertex matrix for patch function.
%
% patch('Faces',f,'Vertices',v,'FaceVertexCData',c,...)
%
% size(f) == [Tf,4] for quadralateral patches
% size(v) == [Tv,2] for 2D space
% size(c) == [Tv,1]
% 
% Let j == f(i,k).  Then the k^th vertex of face i is vertex v(j,:).
%
% -------------------------------------------------------------------------
% For the propeller:
%
Nv = Np+1;         % number of vertices per blade section
Nf = Np;           % number of    faces per blade section

Mv = Mp+1;         % number of vertices radially
Mf = Mp;           % number of    faces radially

Tv =  Mv   * Nv;     % number of vertices in total
Tf = (Mv-1)*(Nv-1);  % number of    faces in total

f = zeros(Tf,4); % faces    (indices of the vertices to use)
v = zeros(Tv,2); % vertices (x,y)
c = zeros(Tv,1); % colors

% -------------------------------------------------------------------------
% The vertices are ordered such that c(j) == c(m,n) for j = m + (n-1)*Mv
v = [X3D(:),Y3D(:)];
c =    s(:);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Create face matrix for patch function- this tells the patch function how
% to connect the vertices to create a face.  This code uses a
% square/rectangular face
%
i = 0;
for m = 1:Mf
    for n = 1:Nf
        % point (m,n) is the start of face i
        i = i + 1;
    
        f(i,1) = m   + (n-1)*Mv;
        f(i,2) = m   + (n  )*Mv;
        f(i,3) = m+1 + (n  )*Mv;
        f(i,4) = m+1 + (n-1)*Mv;
    end
end
% -------------------------------------------------------------------------

% save temp
% disp('hi')
% pause,

% -------------------------------------------------------------------------
% Create Figure
figure1 = figure('Color',[1 1 1]);

% Create axes
% axes('Visible','off','Parent',figure1)%,'CLim',[-1.5e+6 1.5e+6]);
view([0 90]);


%Call patch function
H = patch('Faces',f,'Vertices',v,'FaceVertexCData',c,'FaceColor','interp','FaceLighting','gouraud');
hold on
% Make patch edges transparent, and plot the blade outline
set(H,'EdgeAlpha',0.05)
plot(X3D(:,  1),Y3D(:,  1),'k','LineWidth',2)
plot(X3D(:,Np1),Y3D(:,Np1),'k','LineWidth',2)
plot(X3D(  1,:),Y3D(  1,:),'k','LineWidth',2)
plot(X3D(Mp1,:),Y3D(Mp1,:),'k','LineWidth',2)


% axis('CLim',[-2e+008 2e+008])
colorbar('FontName','Times','FontSize',18)
title(plottitle,'FontName','Times','FontSize',18)
 
grid on,
 box on,
axis equal
set(gca,'CLim',[-1 1],'FontName','Times','FontSize',18)
