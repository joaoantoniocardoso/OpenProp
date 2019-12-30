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


% ------------------------------------------------ Evaluate basis functions
% Created: Brenden Epps, 4/10/10
%
% Evaluates the B-spline basis functions at points t.
%
% Inputs:
%   t       [Mt,1]  vector of field points
%   n       [1, 1]  spine parameter, Ms == n+1 == number of splines
%   k       [1, 1]  spline order, k == 4 for cubic splines
%
% Output:
%   B       [Mt,Ms]  B(i,j) == basis funciton j evaluated at t(i)
%
%   D1      [Mt,Ms] D1(i,j) == 1st derivative of basis funciton j evaluated at t(i)
%   D2      [Mt,Ms] D2(i,j) == 2nd derivative of basis funciton j evaluated at t(i)
%
%   tstar   [Ms,1]  nodes (that represent the "t" location of the vertices 
%                          of the spline) which are given by the averages 
%                          of successive  k-1  knots:
%
%            tstar(i) == ( t(i+1) + ... + t(i+k-1) ) / (k-1)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% To evaluate a B-spline given spline amplitudes:
%   A       [Ms,1]  spline amplitude vector
%   S(t)    [Mt,1]  spline values at t (function that A represents), where
%
%            S    == B*A    
%            S(i) == sum(B(i,1:Ms).*A(1:Ms),   S(i) === S(t(i))
%
% ----------- RUN THIS CODE -----------
% n = 6;         % number of spline segments
% t = 0:0.01:1;  % field points
% 
% [B, D1, D2, knot, tstar] = Bspline_basis(t,n,4);
% 
% A = [1 2 3 4 5 3 2]'; % n+1 spline amplitudes (vertices)
% 
% S = B*A;  
% 
% figure, hold on,
%     plot(tstar,A,'o--g')
%     plot(t    ,S,'r')
% --------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% To approximate a function F(t) with spline S(t):
%   F(t)    [Mt,1] function values at t
%
%            A    == linsolve(B,F)
%            S    == B*A  ~~ F, 
%
% Note if Ms == Mt, then S == F, otherwise S approximates F smoothly
%
% ----------- RUN THIS CODE -----------
% n  = 4;
% k  = 4;
% t  = [0 : 0.1 : 1]';  % field points
% F  = t .* (1-t);      % function
% F1 = 1 - 2*t;         % 1st derivative
% F2 =   - 2;           % 2nd derivative
%    
% [B, D1, D2, knot, tstar] = Bspline_basis(t,n,k);   
% 
% A = linsolve(B,F);  % if the number of splines is less than the number of data points, then possibly need to run: A = pinv(B)*F;
% 
% S  = B *A;
% S1 = D1*A;
% S2 = D2*A;
% 
% % Evaluate the spline on a finer resolution
% tt   = linspace(0,1,100); 
% 
% [BB, DD1, DD2, ~, ttstar] = Bspline_basis(tt,n,k);
% 
% SS  = BB *A;
% SS1 = DD1*A;
% SS2 = DD2*A;
% 
% 
% % Display plots
% figure(1), hold on, grid on, box on,
%     plot(t,F ,'r.-')
%     plot(t,S ,'g.','markersize',16);
%     plot(tt,SS,'k-');
%     plot(tstar, A,'o--k','markersize',10);  % vertices
% 
% figure(2), hold on, grid on, box on,
%     plot(t,F1 ,'r.-')
%     plot(t,S1 ,'g.','markersize',16);
%     plot(tt,SS1,'k-');    
%  
% figure(3), hold on, grid on, box on,
%     plot(t,F2 ,'r.-')
%     plot(t,S2 ,'g.','markersize',16);
%     plot(tt,SS2,'k-');
%     
%     
% figure(4); hold on, grid on, box on,
%     for i = 1:n+1
%         plot(t,B(:,i))
%     end    
% 
%     for i = 1:n+1
%         plot(tt,BB(:,i),'k')
%     end
% --------------------------------------
% -------------------------------------------------------------------------

function [B, D1, D2, knot,tstar] = Bspline_basis(t,n,k,KnotMethod)

if nargin < 4
    KnotMethod = 'Uniform';
end

if length(n) > 1  % then function called as Bspline_basis(t,knot,k)
    knot = n;
    
    Mk   = length(knot); % number of knots 
    
    n    = Mk - k - 1;   % spline parameter
else
    Mk = k+n+1;  % number of knots
    
    
    if strcmp(KnotMethod , 'Uniform')
        % ----------------------------- Form a uniform knot sequence on [0,1]
        knot            = zeros(Mk,1);           % knot sequence: [0,...,0,x,x,x,1,...,1]  
        knot(k:n+2)     = linspace(0,1,n-k+3)';  %  linspace provides:  (0,x,x,x,1)  so  n-k+1 "x" entries, k  "0" entries, and k "1" entries
        knot(n+3:k+n+1) = ones(k-1,1);           %                                   
        % -------------------------------------------------------------------      
        
    elseif strcmp(KnotMethod , 'Cosine')
        % ------------------------- Form a cosine-spaced knot sequence on [0,1]
        knot            = zeros(Mk,1);           % knot sequence: [0,...,0,x,x,x,1,...,1]
        Np              = n-k+3;
        knot(k:n+2)     = 0.5*(  1 - cos( pi*[0:(Np-1)]'/(Np-1) )  );
        knot(n+3:k+n+1) = ones(k-1,1);
        % -------------------------------------------------------------------

        
    elseif strcmp(KnotMethod , 'Bezier')
        % -------------------------------- Form a Bezier knot sequence on [0,1]
        if k ~= n+1, 
            disp('Error: Bezier curve requires  k == n+1'), return, 
        else
            knot            = zeros(Mk,1);           % knot sequence
            knot(k+1:k+n+1) = ones(k,1);
        end
        % -------------------------------------------------------------------
     
        
    elseif strcmp(KnotMethod , 'Average')
        % ------------------------ Form a knot sequence based on site averaging
        if n ~= length(t) - 1, 
            disp('ERROR: site averaging requires  n == length(t) - 1  '), return, 
        else

            knot            = zeros(Mk,1);           % knot sequence
            knot(n+2:k+n+1) = ones(k,1);

            for j = 1:n+1-k
                knot(k+j) = sum(t(j:j+k-2))/(k-1);
            end
        end
    else
        disp('ERROR: unknown KnotMethod')
        return
    end
    
    % ---------------------------------------------------------------------
end
    
maxknot = max(knot);

Mt = length(t);
Ms = n+1;        % number of basis splines

% -------------------------------------------------------------------------
N  = zeros(Mk-1,k,Mt);  % B-spline basis functions of order 1,...,k
B  = zeros(Mt,Ms);      % B-spline basis functions of order k

N1 = zeros(Mk-1,k,Mt);  % B-spline basis functions of order 1,...,k: 1st derivative
N2 = zeros(Mk-1,k,Mt);  % B-spline basis functions of order 1,...,k: 2nd derivative

D1 = zeros(Mt,Ms);      % B-spline basis functions of order k: 1st derivative
D2 = zeros(Mt,Ms);      % B-spline basis functions of order k: 2nd derivative


% Evaluate the B-spline pointwise in t
for point = 1:Mt
    
    % Evaluate order 1 basis functions:
    j = 1;
    
        for i = 1:Mk-1            % for each knot panel [knot(i),knot(i+1)]
            if (knot(i) <= t(point) && t(point) < knot(i+1)) ||  (knot(i) < t(point) && knot(i+1) == maxknot && t(point) == maxknot)
                N(i,j,point) = 1;
            end
        end
        
%         i = Mk-1;
%             if knot(i) == t(point) && t(point) == knot(i+1)
%                 N(i,j,point) = 1;
%             end
    
    % Apply recursion relation to evaluate higher order basis functions:
    for j = 2:k
        % --------------------------------- Evaluate spline basis functions
        for i = 1:Mk-j            % for each knot panel [knot(i),knot(i+1)]

                % % This leads to NaN values when  [ knot(i+j-1) == knot(i) ] or [ knot(i+j) == knot(i+1) ]
                % N(i,j,point) = (t(point)  - knot(i) )/(knot(i+j-1)-knot(i)  ) *  N(i  ,j-1,point) ...
                %               +(knot(i+j) - t(point))/(knot(i+j)  -knot(i+1)) *  N(i+1,j-1,point);
                % 
                % N1(i,j,point) = (j-1                 )/(knot(i+j-1)-knot(i)  ) *  N(i  ,j-1,point) ...
                %               -(j-1                 )/(knot(i+j)  -knot(i+1)) *  N(i+1,j-1,point);
                % 
                % N2(i,j,point) = (j-1                 )/(knot(i+j-1)-knot(i)  ) * N1(i  ,j-1,point) ...
                %               -(j-1                 )/(knot(i+j)  -knot(i+1)) * N1(i+1,j-1,point);
                 
                      
%             if     N(i,j-1,point) ~= 0 && N(i+1,j-1,point) ~= 0
%                 N(i,j,point) = (t(point)  - knot(i) )/(knot(i+j-1)-knot(i)  ) *  N(i  ,j-1,point) ...
%                               +(knot(i+j) - t(point))/(knot(i+j)  -knot(i+1)) *  N(i+1,j-1,point);   
%                       
%                N1(i,j,point) = (j-1                 )/(knot(i+j-1)-knot(i)  ) *  N(i  ,j-1,point) ...
%                               -(j-1                 )/(knot(i+j)  -knot(i+1)) *  N(i+1,j-1,point);
% 
%                N2(i,j,point) = (j-1                 )/(knot(i+j-1)-knot(i)  ) * N1(i  ,j-1,point) ...
%                               -(j-1                 )/(knot(i+j)  -knot(i+1)) * N1(i+1,j-1,point);             
%                       
%             elseif N(i,j-1,point) == 0 && N(i+1,j-1,point) ~= 0
%                 N(i,j,point) = 0 ...
%                               +(knot(i+j) - t(point))/(knot(i+j)  -knot(i+1)) *  N(i+1,j-1,point);   
%                       
%                N1(i,j,point) =-(j-1                 )/(knot(i+j)  -knot(i+1)) *  N(i+1,j-1,point);
% 
%                N2(i,j,point) =-(j-1                 )/(knot(i+j)  -knot(i+1)) * N1(i+1,j-1,point);               
%                       
%             elseif N(i,j-1,point) ~= 0 && N(i+1,j-1,point) == 0
%                 N(i,j,point) = (t(point)  - knot(i) )/(knot(i+j-1)-knot(i)  ) *  N(i  ,j-1,point) ...
%                               + 0;              
%                       
%                N1(i,j,point) = (j-1                 )/(knot(i+j-1)-knot(i)  ) *  N(i  ,j-1,point);
% 
%                N2(i,j,point) = (j-1                 )/(knot(i+j-1)-knot(i)  ) * N1(i  ,j-1,point);     

            if     (knot(i+j-1)-knot(i)) ~= 0  &&   (knot(i+j)-knot(i+1)) ~= 0
                N(i,j,point) = (t(point)  - knot(i) )/(knot(i+j-1)-knot(i)  ) *  N(i  ,j-1,point) ...
                              +(knot(i+j) - t(point))/(knot(i+j)  -knot(i+1)) *  N(i+1,j-1,point);   
                      
               N1(i,j,point) = (j-1                 )/(knot(i+j-1)-knot(i)  ) *  N(i  ,j-1,point) ...
                              -(j-1                 )/(knot(i+j)  -knot(i+1)) *  N(i+1,j-1,point);

               N2(i,j,point) = (j-1                 )/(knot(i+j-1)-knot(i)  ) * N1(i  ,j-1,point) ...
                              -(j-1                 )/(knot(i+j)  -knot(i+1)) * N1(i+1,j-1,point);             
                      
            elseif (knot(i+j-1)-knot(i)) == 0 && (knot(i+j)-knot(i+1)) ~= 0
                N(i,j,point) = 0 ...
                              +(knot(i+j) - t(point))/(knot(i+j)  -knot(i+1)) *  N(i+1,j-1,point);   
                      
               N1(i,j,point) =-(j-1                 )/(knot(i+j)  -knot(i+1)) *  N(i+1,j-1,point);

               N2(i,j,point) =-(j-1                 )/(knot(i+j)  -knot(i+1)) * N1(i+1,j-1,point);               
                      
            elseif (knot(i+j-1)-knot(i)) ~= 0 && (knot(i+j)-knot(i+1)) == 0
                N(i,j,point) = (t(point)  - knot(i) )/(knot(i+j-1)-knot(i)  ) *  N(i  ,j-1,point) ...
                              + 0;              
                      
               N1(i,j,point) = (j-1                 )/(knot(i+j-1)-knot(i)  ) *  N(i  ,j-1,point);

               N2(i,j,point) = (j-1                 )/(knot(i+j-1)-knot(i)  ) * N1(i  ,j-1,point); 
               
               
               
            %elseif N(i,j-1,point) == 0 && N(i+1,j-1,point) == 0
            %    N(i,j,point) = 0;   
            end
            
        end
        
%         % ----------------------------- Evaluate basis function derivatives
%         for i = 1:Mk-j            % for each knot panel [knot(i),knot(i+1)]
%             N1(i,j,point) = (j/(knot(i+j)-knot(i)))*N(i,j-1,point)
%             
%             
%         end
    end

    % ------- Record spline basis functions
     B(point,:) =  N(1:Ms,k,point);
    D1(point,:) = N1(1:Ms,k,point);
    D2(point,:) = N2(1:Ms,k,point);
    
    % ------- Record spline derivatives
     
    
    
end

% if t(end) == 1
%     B(end,end) = 1;
% end
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Determine the "nodes", which are the averages of successive  k-1  knots.
% In other words, the nodes are the points:
%
%       tstar(i) = ( t(i+1) + ... + t(i+k-1) ) / (k-1)
%
%   which represent the location of the vertices of the spline.
%
% Ms == Mk - k == n + 1

tstar = zeros(Ms,1);

for i = 1:Ms
    
    tstar(i) = mean( knot(i+[1:k-1]) );
end

