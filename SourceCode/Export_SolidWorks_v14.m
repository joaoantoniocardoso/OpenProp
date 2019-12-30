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

function [] = Export_SolidWorks_v14(filename_SolidWorks,Np,Mp,Z,X3D,Y3D,Z3D)


fid = fopen(filename_SolidWorks,'wt');  % 1-Aug-2013 BEPPS: 'w' changed to 'wt'  so newline character '\n' appears properly on Windows machines

% Prop Parameters at beginning of file
fprintf(fid,'%g, ' ,Np); 
fprintf(fid,'%g, ' ,Mp);
fprintf(fid,'%g,\n',Z);


    % Output curves defining each 2D section along the span
    % for each section along the span
    for i = 1:Mp+1       
        fprintf(fid,strcat('SectionCurve',num2str(i),',\n'));
        
        % for each point along the suction and pressure surfaces
        % (trailing edge -> leading edge -> trailing edge, close the curve)
%         for j = [1:Np,Np+2:2*Np,1] % (2*Np   points) does not double print the leading edge but does double print the trailing edge to close the curves        
%         for j = [1:Np,Np+2:2*Np-1,1] % (2*Np-1 points) does not double print the leading edge but does double print the trailing edge to close the curves                    
        for j = [1:Np,Np+2:2*Np] % (2*Np-1 points) does not double print the leading edge
%         for j = 1:2*(Np-1)+NLE % each curve contains (2*(Np-1)+NLE) points
            fprintf(fid,'%f,%f,%f,\n',X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1));
        end
    end
    
   
    % Make guide curves
    n = 0;   
    % for 7 points along the chord
    for j = [1 floor(Np/3) floor(2*Np/3) Np floor(4*Np/3) floor(5*Np/3) 2*Np];
%     for j = [1 floor(1*(2*(Np-1)+NLE)/6) floor(2*(2*(Np-1)+NLE)/6) ...
%                floor(3*(2*(Np-1)+NLE)/6) floor(4*(2*(Np-1)+NLE)/6) ...
%                floor(5*(2*(Np-1)+NLE)/6) floor(6*(2*(Np-1)+NLE)/4)];   
        n = n + 1;
        
        
        fprintf(fid,strcat('GuideCurve',num2str(n),',\n'));
%         for i = 1:Mp  % for each section along the span except the last one
        for i = 1:Mp+1  % for each section along the span
            
%             if i == Mp+1 && j == 2*Np
%                 fprintf(fid,'%f,%f,%f',X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1));
%                 continue
%             end
            
               
%         plot3(X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1),'.b','markersize',20)
%         pause,
            
            fprintf(fid,'%f,%f,%f,\n',X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1));
        end
        
    end
    

    
    % Output duplicate trailing edge guide curves:    
    % Guide curve 1:
    fprintf(fid,'TEGuideCurve1,\n');    
        j = 1;
    for i = 1:Mp+1  % for each section along the span        
        fprintf(fid,'%f,%f,%f,\n',X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1));
    end

    % Guide curve 7:
    fprintf(fid,'TEGuideCurve7,\n');    
        j = 2*Np;
    for i = 1:Mp+1  % for each section along the span        
        fprintf(fid,'%f,%f,%f,\n',X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1));
    end    
    
    
    % Output duplicate tip section profile:
    i = Mp+1;
    fprintf(fid,strcat('TipSectionCurve',num2str(i),',\n'));        
    % for each point along the suction and pressure surfaces
    % (trailing edge -> leading edge -> trailing edge)
    for j = [1:Np,Np+2:2*Np] % (2*Np-1 points) does not double print the leading edge
        fprintf(fid,'%f,%f,%f,\n',X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1));
    end    
    
    % Output tip curves
    for j = 1:Np-2
%     for j = 1:(Np-1+(NLE-1)/2-1)
        fprintf(fid,strcat('TipCurve',num2str(j),',\n'));
        i=Mp+1;
        fprintf(fid,'%f,%f,%f,\n',X3D(i,   1+j,1),Y3D(i,   1+j,1),Z3D(i,   1+j,1));
        fprintf(fid,'%f,%f,%f,\n',X3D(i,2*Np-j,1),Y3D(i,2*Np-j,1),Z3D(i,2*Np-j,1));
%         fprintf(fid,'%f,%f,%f,\n',X3D(i,             1+j,1),Y3D(i,             1+j,1),Z3D(i,             1+j,1));
%         fprintf(fid,'%f,%f,%f,\n',X3D(i,(2*(Np-1)+NLE)-j,1),Y3D(i,(2*(Np-1)+NLE)-j,1),Z3D(i,(2*(Np-1)+NLE)-j,1));
    end
    
    
    % Output duplicate root section profile:
    i = 1;
    fprintf(fid,strcat('RootSectionCurve',num2str(i),',\n'));        
    % for each point along the suction and pressure surfaces
    % (trailing edge -> leading edge -> trailing edge)
    for j = [1:Np,Np+2:2*Np] % (2*Np-1 points) does not double print the leading edge
        fprintf(fid,'%f,%f,%f,\n',X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1));
    end
    
    
    % Output root curves
    for j = 1:Np-2
%     for j = 1:(Np-1+(NLE-1)/2-1)
        fprintf(fid,strcat('RootCurve',num2str(j),',\n'));
        i=1;
        fprintf(fid,'%f,%f,%f,\n',X3D(i,   1+j,1),Y3D(i,   1+j,1),Z3D(i,   1+j,1));
        fprintf(fid,'%f,%f,%f,\n',X3D(i,2*Np-j,1),Y3D(i,2*Np-j,1),Z3D(i,2*Np-j,1));
%         fprintf(fid,'%f,%f,%f,\n',X3D(i,             1+j,1),Y3D(i,             1+j,1),Z3D(i,             1+j,1));
%         fprintf(fid,'%f,%f,%f,\n',X3D(i,(2*(Np-1)+NLE)-j,1),Y3D(i,(2*(Np-1)+NLE)-j,1),Z3D(i,(2*(Np-1)+NLE)-j,1));
    end
    
    
    % Output trailing edge curves for each 2D section along the span
    for i = 1:Mp+1         
        fprintf(fid,strcat('TECurve',num2str(i),',\n'));
        j=1;
        fprintf(fid,'%f,%f,%f,\n',X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1));
        j=2*Np;
%         j = 2*(Np-1)+NLE;
        fprintf(fid,'%f,%f,%f,\n',X3D(i,j,1),Y3D(i,j,1),Z3D(i,j,1));
    end    

    fclose(fid);
end % function