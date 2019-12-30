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
% =================================================== RepairSpline Function
% Created: Brenden Epps, 12/13/2010
%
% This function repairs a discontinuous function by fitting a low-order
% spline to the data.
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
% -------------------------------------------------------------------------

% % -------------------------------------------------------------------------
% % ---------------------------------------------------------- Implementation
% % Spline inputs:   
% %       n  = 5;  number of splines Ms == n + 1, NOTE: n must be ODD or else spline is always zero at TC == 0!
% %       k  = 4;  polynomial order, k == 4 for cubic B-splines
% %
% % ---------------------- Create smooth B-spline represntation of unknown, X
% % X == size [1,Mp]      variable
% % A == size [Ms,1]      B-spline amplitudes
% % B == size [Mp,Ms]     B-spline basis functions
% %
% % X(1,TC) = ( B(TC,spline) * A(spline) )'
% %
% % A = (Binv * X');
% %
% % Xnew = (B * A)' == ( B * (Binv*X') )' == (Binv*X')' * B' == X * Binv'*B'
% %
% % Xnew = X * Bsmooth
% % -------------------------------------------------------------------------

% % ------------------------------------------------ Implement RepairSpline.m
% % To smooth the data X(RC):  X_smooth = (SplineBasis*pinv(SplineBasis)*X(:))';
% %
% TC = (RC(:) - RC(1)) / (RC(end) - RC(1)); % Map RC into TC in the interval [0,1]
% n = 5;                                    % number of splines == n + 1
% k = 4;                                    % polynomial order, k == 4 for cubic B-splines
% SplineBasis = Bspline_basis(TC,n,k);      % B-spline basis functions [Mp,n+1]
% %
% Bsmooth = (SplineBasis*pinv(SplineBasis))';
% %
% % Now: X_smooth = X * Bsmooth;
% % -------------------------------------------------------------------------


function [XX] = RepairSpline(RC,X,name,H)

[a,b] = size(X);

if     nargin == 2
    name = '[unknown]';
    H    = 0;
        
elseif nargin == 3,
    H    = 0;
end

RC = RC(:);
X  =  X(:);

% disp(['Repairing: ',name])

Mp   = length(RC);



% ----------------------------------------------------------- Spline inputs   
n  = 3;                      % number of splines Ms == n + 1, NOTE: n must be ODD or else spline is always zero at TC == 0!
k  = 4;                      % polynomial order, k == 4 for cubic B-splines
Nk = k+n+1;                  % number of knots
% -------------------------------------------------------------------------


% ------------------------------- Find spline parameter TC, as in S = S(TC)
% i.e. Map RC into TC in the interval [0,1], where RC is assumed to be 
%      given in equal intervals.

TC = (RC - RC(1)) / (RC(end) - RC(1));                        % size [Mp,1]
% -------------------------------------------------------------------------

% % % -------------------------------------------------------------------------
% % % ------------------------ Find spline knot sequence using averaging method
% % knot         = zeros(Nk,1);
% % knot(n+2:Nk) = ones(k,1);
% % 
% % for j = 1:n+1-k
% %     knot(k+j) = sum(TC(1+(j:j+k-2)))/(k-1);
% % end
% % % -------------------------------------------------------------------------


% ----------------------------------- Evaluate the B-spline basis functions
% % [B] = Bspline_basis(TC,knot,k);  % B-spline basis functions [Mp,Ms], B(TC,spline)

[B] = Bspline_basis(TC,n,k);    % B-spline basis functions [Mp,Ms], B(TC,spline)

% figure; hold on, grid on, box on,
%     for i = 1:n+1
%         plot(TC,B(:,i))
%     end
% -------------------------------------------------------------------------


% ----------------------- Solve linear system for spline curve amplitudes
%
% X(TC) = B(TC,spline) * A(spline)
%
% ---------------------------------- Solve for spline amplitudes, A(spline)
% A = linsolve(B,X)  % this only works if the number of splines equals the number of data points (i.e. no data smoothing)
A = pinv(B)*X;

% ----------------------------------------- Interpolate and smooth the data
XX = B*A;
% -------------------------------------------------------------------------

% ------------------------------- Evaluate the spline on a finer resolution
TT = linspace(0,1,100);  % fine grid of spline parameters

[BB] = Bspline_basis(TT,n,k);  % size [length(TT),Ms]

RRR = TT * (RC(end) - RC(1)) + RC(1);
XXX = BB*A;
% -------------------------------------------------------------------------


% ------------------------------------------------------------ Display plot
if H ~= 0
    figure(H), 
        hold on,
        H1 = plot(RRR,XXX,'k');
        H2 = plot(RC,X ,'c.');
        H3 = plot(RC,XX,'g.');

    pause(0.1),
%     pause,
    set(H1,'visible','off')
    set(H2,'visible','off')
    set(H3,'visible','off')
end
%--------------------------------------------------------------------------


%----------------------------------------------------- Reshape if necessary
[aa,bb] = size(XX);

if aa == b & bb == a
    XX = XX';
end
%--------------------------------------------------------------------------

% ===================================================== END Repair Function


