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
% Z3D == [Mp1,Np]
% s   == [Mp1,Np]
%
% -------------------------------------------------------------------------

function [] = Plot_Blade_Contours2(X3D,Y3D,Z3D,s,plottitle)

[Mp1,TwoNp] = size(X3D);  % X3D is from the geometry.m module

Mp = Mp1-1;   
Np = TwoNp/2;

% Fix bug: add an additional row, so s is the stress at each vertex
s(Mp+1,:) = 0*s(Mp,:);

% -------------------------------------------------------------------------
% Fix buggy issue with current GeometryBlade code:  
%   Leading edge is given twice in (X3D,Y3D,Z3D).  
%   Therefore, they are size (Mp+1,2*Np), but they should be (Mp+1,2*Np-1)
%
X3D = X3D(:,[1:Np,Np+2:2*Np]);  % now size (Mp+1,2*Np-1)
Y3D = Y3D(:,[1:Np,Np+2:2*Np]);  % now size (Mp+1,2*Np-1)
Z3D = Z3D(:,[1:Np,Np+2:2*Np]);  % now size (Mp+1,2*Np-1)
  s =   s(:,[1:Np,Np+2:2*Np]);  % now size (Mp+1,2*Np-1)
% -------------------------------------------------------------------------

  

% -------------------------------------------------------------------------
% Concatenate matrices to create vertex matrix for patch function.
%
% patch('Faces',f,'Vertices',v,'FaceVertexCData',c,...)
%
% size(f) == [Tf,4] for quadralateral patches
% size(v) == [Tv,3] for 3D space
% size(c) == [Tv,1]
% 
% Let j == f(i,k).  Then the k^th vertex of face i is vertex v(j,:).
%
% -------------------------------------------------------------------------
% For the propeller:
%
Nv = 2*Np-1;         % number of vertices per blade section
Nf = 2*Np-2;         % number of    faces per blade section

Mv =   Mp+1;         % number of vertices radially
Mf =   Mp;           % number of    faces radially

Tv =  Mv   * Nv;     % number of vertices in total
Tf = (Mv-1)*(Nv-1);  % number of    faces in total

f = zeros(Tf,4); % faces    (indices of the vertices to use)
v = zeros(Tv,3); % vertices (x,y,z)
c = zeros(Tv,1); % colors

% -------------------------------------------------------------------------
% The vertices are ordered such that c(j) == c(m,n) for j = m + (n-1)*Mv
v = [X3D(:),Y3D(:),Z3D(:)];
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


% -------------------------------------------------------------------------
% Create Figure
figure1 = figure('Color',[1 1 1]);
% 
% % Create axes
axes('Visible','off','Parent',figure1)%,'CLim',[-1.5e+6 1.5e+6]);
view([90 2]);
% axis off

% grid('on');

%Call patch function
patch('Faces',f,'Vertices',v,'FaceVertexCData',c,'FaceColor','interp','FaceLighting','gouraud')
% axis('CLim',[-2e+008 2e+008])
colorbar('FontWeight','bold')
title(plottitle,'FontSize',14,'FontWeight','bold')
hold on, axis equal
