% -------------------------------------------------------------------------
% Created: Brenden Epps, 8/12/10
%
% Make SolidWorks.txt files, with coordinates for a single blade.
%
% Use this with SolidWorks macro v18.
%
% Blade geometry:
%   X3D(i,j,k) [m], X position in 3D space
%   Y2D(i,j,k) [m], Y position in 3D space
%   Z3D(i,j,k) [m], Z position in 3D space
%
%   i = 1:Mp+1      % for each section along the span
%   j = 1:2*Np      % for each point   along the upper and lower surfaces
%   k = 1:Z         % for each blade
%
% -------------------------------------------------------------------------

function [] = Duct_Export_SolidWorks_v18_OneLine(filename_SolidWorks,Np,X3D,Y3D,Z3D)


fid = fopen(filename_SolidWorks,'w');

[Np,Nt] = size(X3D);  
 Np     = Np/2;


fprintf(fid,'%g, ' ,Np);

   
    fprintf(fid,'DuctSection, ');

    % for each point along the suction and pressure surfaces
    % (trailing edge -> leading edge -> trailing edge, close the curve)
    for j = [1:Np,Np+2:2*Np-1,1] % (2*Np-1 points) does not double print the leading edge
        fprintf(fid,'%f,%f,%f, ',X3D(j,1),Y3D(j,1),Z3D(j,1));
    end


fclose(fid);

end % function