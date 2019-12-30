% -------------------------------------------------------- Make Rhino files
% Modified: 9/25/09 by Jordan Stanway and Brenden Epps
% 
% This code makes a script that you can run in Rhino.  Here are the steps:
%   1) At the Command: prompt, type "ReadCommandFile"
%           -- Locate the script file, e.g. OpenProp_RhinoProp.txt 
%   2) When the Document Properties window opens, set:
%           -- Model Units: meters
%           -- Absolute tolerance: 0.00001
%           -- Relative tolerance: 0.1
%           -- Angle    tolerance: 0.1
%   3) Watch as Rhino reads in all the points, makes each section, fills 
%      the tip section, lofts the remaining sections, joins the surfaces,
%      makes Z blades from the key blade, and makes the hub
%   4) If your propeller does not loft or join automatically, then try
%      increasing or decreasing the tolerance values.
% 
%
% Make _RhinoBlade.txt, with coordinates for a single blade and
% commands to make Z blades
    
function [ ] = Export_Rhino_v1(filename_Rhino,rake,R,XR,skew0,Mp,Np,X3D,Y3D,Z3D,Z,Rhub) 


        sprintf(['This generates a script that you can run in Rhino. \n\n',...
          '1) At the Command: prompt, type "ReadCommandFile" \n',...
          '  -- Locate the script file, e.g. OpenProp_RhinoProp.txt \n',...
          '2) When the Document Properties window opens, set: \n',...
          '  -- Model Units: meters \n',...
          '  -- Absolute tolerance: 0.00001 \n',...
          '  -- Relative tolerance: 0.1 \n',...
          '  -- Angle    tolerance: 0.1 \n',...
          '3) Watch as Rhino reads in all the points, makes each section,\n',... 
          '   fills the tip section, lofts the remaining sections, joins\n',...
          '   the surfaces, makes Z blades from the key blade, and makes\n',...
          '   the hub. \n',...
          '4) If your propeller does not loft or join automatically, try\n',...
          '   increasing or decreasing the tolerance values. \n'])
    
      
    
    % Initialize file
    fid = fopen(filename_Rhino,'wt');   % 1-Aug-2013 BEPPS: 'w' changed to 'wt'  so newline character '\n' appears properly on Windows machines
    
    
    fprintf(fid,'!_SetActiveViewport Perspective\n');
    
    % In order for the surface to loft correctly, you probably will
    % need to manually set the "absolute precision" of Rhino to be
    % "10^-5 units" and manually change the "model units" to meters.  
    % Note: the coordinates output from OpenProp are in meters.
    % This command should pause the script fro
    fprintf(fid,'_DocumentPropertiesPage Units \n');
          
    % Define Rhino curve type (choose one)
    % curve_cmd   = 'Curve \n';
    curve_cmd   = 'InterpCrv \n';

    
    % Compute where the blade tip should be
    tip_x = -rake(end) - R*pchip(XR,skew0,1)*(pi/180);
    tip_y =         R*sind(pchip(XR,skew0,1));
    tip_z =         R*cosd(pchip(XR,skew0,1));
    tip = [tip_x, tip_y, tip_z];    
 
    for i = 1:Mp+1                        % For each section along the span
        fprintf(fid,curve_cmd);               % Print curve command in file
            
        % For each point   along the upper and lower surfaces:
        for j = [1:Np,Np+2:2*Np] % (2*Np-1 points) does not double print the leading edge
            fprintf(fid,'%.9f,%.9f,%.9f\n',X3D(i,j),Y3D(i,j),Z3D(i,j));  % print to file with 9 decimal places
        end
        
        % If the first and last points in the section are identical, then 
        % do nothing, else close the curve by adding another point the 
        % same as the first one.
        if strcmp(sprintf('%.9f,%.9f,%.9f\n',X3D(i,1)  ,Y3D(i,1)  ,Z3D(i,1)),...
                  sprintf('%.9f,%.9f,%.9f\n',X3D(i,end),Y3D(i,end),Z3D(i,end)))
            % disp(sprintf('%i start and end are identical, not adding point to close', i));
        else
            fprintf(fid,'%.9f,%.9f,%.9f\n',X3D(i,1),Y3D(i,1),Z3D(i,1));
        end
    end
    
    % Extrude the tip section curve to the "tip" point
    fprintf(fid,'SelNone\n');
    fprintf(fid,'SelLast\n');
    fprintf(fid,'ExtrudeCrv Mode=ToPoint \n');
    fprintf(fid,'%.9f,%.9f,%.9f\n',tip(1),tip(2),tip(3));
    fprintf(fid,'enter\n');

    % Loft the other sections to the tip section
    fprintf(fid,'SelNone\n');
    fprintf(fid,'SelClosedCrv\n');
    fprintf(fid,'-Loft Type=Tight Simplify=None \n');
    fprintf(fid,'enter\n');
    fprintf(fid,'enter\n');
    fprintf(fid,'enter\n');
    fprintf(fid,'SelNone\n');
    fprintf(fid,'SelSrf\n');
    fprintf(fid,'Join\n');
    fprintf(fid,'SelNone\n');
    fprintf(fid,'Zoom All Extents\n');
    fprintf(fid,'enter\n');
 
    
    % ---------- Commands to make Z blades:
    fprintf(fid,'SelPolysrf\n');
    fprintf(fid,'Rotate3D\n');
    fprintf(fid,'0,0,0\n');
    fprintf(fid,'1,0,0\n');
    fprintf(fid,'Copy=Yes\n');
    % copy blades
    for k=2:Z
        fprintf(fid,'%f\n',(k-1)*(360/Z));
    end
    fprintf(fid,'enter\n');
    
    
    % ----------- Commands to make hub
    fprintf(fid,'Circle Vertical 0,0,0 \n');
    fprintf(fid,'%f \n',Rhub);
    fprintf(fid,'0,0,%f\n',Rhub);
    fprintf(fid,'0,%f,0\n',Rhub);
    % to choose direction of the circle
        % ** this doesn't seem to work all the time... :-(
    fprintf(fid,'SelNone\n');
    fprintf(fid,'SelLast\n');
    fprintf(fid,'ExtrudeCrv BothSides=Yes Cap=Yes DeleteInput=Yes \n');
    
    Lhub = 2*R;
    fprintf(fid,'%f \n',Lhub);
    fprintf(fid,'Zoom All Extents\n');
    fprintf(fid,'enter\n');
    
    
    
    % Close the file:
    fclose(fid);
    
end