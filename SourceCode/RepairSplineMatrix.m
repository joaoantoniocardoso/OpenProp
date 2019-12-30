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
% MIT OpenProp_v3.2.0
% Last modified: 11/17/2011 Brenden Epps
% =========================================================================
%
% This function creates the smoothing matrix Bsmooth.
%
%--------------------------------------------------------------------------

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



function [Bsmooth] = RepairSplineMatrix(RC)


% ------------------------------------------------ Implement RepairSpline.m
% To smooth the data X(RC):  X_smooth = (SplineBasis*pinv(SplineBasis)*X(:))';
%
TC = (RC(:) - RC(1)) / (RC(end) - RC(1)); % Map RC into TC in the interval [0,1]
n = 5;                                    % number of splines == n + 1
k = 4;                                    % polynomial order, k == 4 for cubic B-splines
SplineBasis = Bspline_basis(TC,n,k);      % B-spline basis functions [Mp,n+1]
%
Bsmooth = (SplineBasis*pinv(SplineBasis))';
%
% Now: X_smooth = X*Bsmooth;
% -------------------------------------------------------------------------