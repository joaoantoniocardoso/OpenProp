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
% ========================================================= Repair Function
%
% This function repairs a discontinuous function.  It searches for
% discontinuous points using interp1(...) and replaces them as necessary.
%
% Input variables:
%   RC    = control point radii / R
%   X     = {G, UASTAR, UTSTAR, BETAIC} function to be repaired
%   name  = variable name string
%   H     = variable plot handle (H == 0 for no plot)
%
% Output variable:
%   XX    = repaired variable
%
% 
% This version deletes points if too much different than interp1(...) and
% replaces them using interp1(...)
%
% -------------------------------------------------------------------------

function [XX] = Repair(RC,X,name,H)

if     nargin == 2
    name = '[unknown]';
    H    = 0;
        
elseif nargin == 3,
    H    = 0;
end

% disp(['Repairing: ',name])

Mp   = length(RC); 
Xrms = sqrt(mean(X.^2));
tolerance = (mean(diff(RC))/sqrt(mean(RC.^2)))^2;

% --------------------------------------------- Identify leading bad values
not_smooth = 1;
m  = 1;

while not_smooth == 1 && m < Mp
    
    Xpredicted = interp1(RC(m+1:end),X(m+1:end),RC(m),'spline','extrap');

    if Xrms < 0.1 
        error = abs( Xpredicted - X(m) );
    else
        error = abs((Xpredicted - X(m))/Xrms);    
    end
    
    if error < tolerance
        not_smooth = 0;  % then the function IS smooth at RC(m)
        
        Mbad = m - 1;               % the previous m was the last bad value
    else
        % disp(['    ',name,'(',num2str(m),') is discontinuous.'])
        
        m = m + 1;                  % check the next point
    end
end

% if Mbad == 0
%     disp(['    Note: ',name,'(1)  is continuous.'])
% end
%
% -------------------------------------------------------------------------

% -------------------------------------------- Identify trailing bad values
not_smooth = 1;
n  = 1;

while not_smooth == 1 && n < Mp
    
    Xpredicted = interp1(RC(Mbad+1:end-n),X(Mbad+1:end-n),RC(n),'spline','extrap');

    if Xrms < 0.1 
        error = abs( Xpredicted - X(n) );
    else
        error = abs((Xpredicted - X(n))/Xrms);    
    end

    if error < tolerance
        not_smooth = 0;   % then the function IS smooth at RC(n)
        
        Nbad = n - 1;               % the previous n was the last bad value
    else
        % disp(['    ',name,'(',num2str(Mp + 1 - n),') is discontinuous.'])
        
        n = n + 1;                  % check the next point
    end
end

% if Nbad == 0
%     disp(['    Note: ',name,'(',num2str(Mp),') is continuous.'])
% end
% -------------------------------------------------------------------------

% ------------------------------------------------------ Display bad values
if H ~= 0
    figure(H), 
        hold on,
        H1 = plot(RC              ,X              ,'c.');
        H2 = plot(RC(1:Mbad)      ,X(1:Mbad)      ,'k.');
        H3 = plot(RC(Mp+1-Nbad:Mp),X(Mp+1-Nbad:Mp),'k.');
end
% -------------------------------------------------------------------------

% ------------------------------------------------------ Replace bad values
% XXpp = csaps(RC(Mbad+1:end-Nbad),X(Mbad+1:end-Nbad),0.99);
% XX   = ppval(fnxtr(XXpp),RC);
% Comment the above two lines and uncomment the following line if you do
% not have the Matlab Spline Toolbox (which has the csaps function)
XX = interp1(RC(Mbad+1:end-n),X(Mbad+1:end-n),RC,'spline','extrap');
% -------------------------------------------------------------------------
    
% ------------------------------------------------------ Display new values
if H ~= 0
    figure(H),
        H4 = plot(RC,XX,'g.');
    %     H4 = plot(RC(1:Mbad)      ,XX(1:Mbad)      ,'g.');
    %     H5 = plot(RC(Mp+1-Nbad:Mp),XX(Mp+1-Nbad:Mp),'g.');    

    pause(0.1),
%     disp('Press any button to continue repair.')
%     pause,
    set(H1,'visible','off')
    set(H2,'visible','off')
    set(H3,'visible','off')
    set(H4,'visible','off')
%     set(H5,'visible','off')
end
%--------------------------------------------------------------------------
% ===================================================== END Repair Function


