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

% -------------------------------------------------------------------------
% 
% Created: Brenden Epps, 2/18/2011
%
% This function fits a cubic B-spline to a given chord distribution and 
% interpolates this chord at the desired radii.
%
% INPUTS
%   XR      [1,Mx] r/R given radii
%   XCoD    [1,Mx] c/D given chord at XR
%   RG      [1,Mp] r/R locations to interpolate chord 
%
% OUTPUTS
%   CoD     [1,Mp] c/D at RG
%
% -------------------------------------------------------------------------

function CoD = InterpolateChord(XR,XCoD,RG)


if (abs(XR(end)-1) < 1e-4) && (XCoD(end) > 0.01)  % if XR == 1 and XCoD > 0.01 at tip  (i.e. finite chord at the tip)

    CoD  = pchip(XR,XCoD, RC);   % allow for finite chord at tip
    
else
    % Assume near-zero chord at tip, and interpolate using a B-spline curve fit to the expanded blade outline
    
    [M,N] = size(RG);


    XR   =   XR(:)';
    XCoD = XCoD(:)';
    RG   =   RG(:);

    % Number of input radii
    Mx = length(XR);


    if XR(Mx) == 1

        % Record tip chord length
        CoDtip   = XCoD(Mx);


        % Set tip chord lengtht to zero to form spline
        XCoD(Mx) = 0;    

    else
        % Create another data point at the tip:
          XR(Mx+1) = 1;
        XCoD(Mx+1) = 0; 
            CoDtip = 0;

            Mx     = Mx+1;
    end

    % -------------------------------------------------------------------------
    % Put all the numbers in one list: this is data to fit spline to
    xr2 = [XR  ,   XR(Mx-1:-1:1)]';
    xc2 = [XCoD,-XCoD(Mx-1:-1:1)]';


    % ------------------------------------------------------------------------- 
    % ----------------------------------------------------------- Spline inputs   
    Md = length(xr2);  % number of  spanwise data sites, Mx == m+1
    m  = Md - 1;      % Mx spline basis functions / control points
    k  = 4;           % polynomial order (k == 4 for cubic spline)
    Mk = k+m+1;       % number of  spanwise knots


    % ------------------------- Find spline parameters using centripital method
    dseg = ( diff(xr2).^2 + diff(xc2).^2 ).^(1/4);

    dtot = sum(dseg);

    ubar = [0;cumsum(dseg)]/dtot;    % size [Md,1]        


    % -------------------------------- Find spline knots using averaging method
    uknot         = zeros(Mk,1);  % spanwise knot sequence
    uknot(m+2:Mk) = ones(k,1);

    for j = 1:m+1-k, 
        uknot(k+j) = sum(ubar(1+(j:j+k-2)))/(k-1);
    end


    % ----------------------------------- Evaluate the B-spline basis functions
    [BC, DC] = Bspline_basis(ubar,uknot,k);    % size [Md,Md], BC(ubar,uspline)


    % ------------------------------- Solve linear system for spline amplitudes
    % Let: 
    %       xr2(ubar) = BC(ubar,uspline) * Axr2(uspline)
    % and
    %       xc2(ubar) = BC(ubar,uspline) * Axc2(uspline)
    %
    % where:
    %
    %   size(xr2)  = [Md,1], data matrix 
    %   size(xc2)  = [Md,1], data matrix 
    %   size(BC)   = [Md,Md], Md == length(ubar) == number of u-splines 
    %   size(Axr2) = [Md,1], xr2 spline amplitudes 
    %   size(Axc2) = [Md,1], xc2 spline amplitudes 
    %
    % --------------------------------- Solve for spline amplitudes, P(uspline)
    Axr2 = linsolve(BC,xr2);
    Axc2 = linsolve(BC,xc2);

    %
    % % -------------------------------------------------------- Check solution
    % xr2 - BC * Axr2
    % xc2 - BC * Axc2
    % -------------------------------------------------------------------------

    % -------------------------------------------------------------------------
    % ------------------------------------- Evaluate spline on finer resolution
    Mdd = 401;

    ubarFINE = linspace(0,1,Mdd)';

    [BC, DC] = Bspline_basis(ubarFINE,uknot,k);    % size [Md,Md], BC(ubar,uspline)

    xr3 = BC * Axr2;
    xc3 = BC * Axc2;

    % Fine resolution spline radius and chord length data: (xr,xc)
    xr =   xr3(1:(Mdd-1)/2+1);
    xc =   xc3(1:(Mdd-1)/2+1);

    %         % Check blade is symmetric        
    %         xr3(1:(Mdd-1)/2+1) - xr3(Mdd:-1:(Mdd-1)/2+1)
    %         xc3(1:(Mdd-1)/2+1) + xc3(Mdd:-1:(Mdd-1)/2+1)



    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    % Interpolate spline values at desired radii
    % CoDraw  = pchip1d(xr,xc ,RG);
    CoDraw  = interp1(xr,xc ,RG,'pchip','extrap');




    % Offset distribution, used to modify chord at the tip
    CoDoffset = (RG-RG(1))/(1-RG(1)) * CoDtip;


    % Final chord distribution, CoD(RG)
    CoD = sqrt(CoDraw.^2+CoDoffset.^2);


    % Reshape CoD if necessary
    if size(CoD,1) == N
        CoD = CoD';
    end
    % -------------------------------------------------------------------------


    % fig; 
    % 
    %     plot(xr,xc,'g','linewidth',2)   % output (r,c) data
    % 
    %     plot(RG,CoDraw,'.b')
    %     plot(XR,XCoD,'.k')        
    %     plot(RG,CoD,'.r')
    % 
    %     axis equal

    % -------------------------------------------------------------------------

end % if (abs(XR(end)-1) < 1e-4) && (XCoD(end) > 0.01)  % if XR == 1 and XCoD > 0.01 at tip  (i.e. finite chord at the tip)

%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function [B, D1, D2, knot] = Bspline_basis(t,n,k)

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
%   B       [Mt,Ms] == [Mt,n+1], where
%                         
%            B(i,j) == basis funciton j evaluated at t(i)
%
% To evaluate a B-spline given spline amplitudes:
%   A       [Ms,1]  spline amplitude vector
%   S(t)    [Mt,1]  spline values at t (function that A represents), where
%
%            S    == B*A    
%            S(i) == sum(B(i,1:Ms).*A(1:Ms),   S(i) === S(t(i))
%
% To approximate a function F(t):
%   F(t)    [Mt,1] function values at t
%
%            A    == linsolve(B,F)
%            S    == B*A  ~~ F, 
%
% Note if Ms == Mt, then S == F, otherwise S approximates F smoothly
%
% -------------------------------------------------------------------------

if length(n) > 1  % then function called as Bspline_basis(t,knot,k)
    knot = n;
    
    Mk   = length(knot); % number of knots 
    
    n    = Mk - k - 1;   % spline parameter
else
    Mk = k+n+1;                              % number of knots
    
    % ----------------------------- Form a uniform knot sequence on [0,1]
    knot            = zeros(Mk,1);           % knot sequence
    knot(k:n+2)     = linspace(0,1,n-k+3)';
    knot(n+3:k+n+1) = ones(k-1,1);
    % -------------------------------------------------------------------
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

    % Apply recursion relation to evaluate higher order basis functions:
    for j = 2:k
        % --------------------------------- Evaluate spline basis functions
        for i = 1:Mk-j            % for each knot panel [knot(i),knot(i+1)]

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
               
            end
            
        end
    end

    % ------- Record spline basis functions
     B(point,:) =  N(1:Ms,k,point);
    D1(point,:) = N1(1:Ms,k,point);
    D2(point,:) = N2(1:Ms,k,point);
      
end
% -------------------------------------------------------------------------




% -------------------------------------------------------------------------  
% ------------------------------------------------------------------------- 
function vv = pchip1d(x,v,xx)

% ----------------------------------------------------------------- pchip1d
%
% This function should replicate Matlab PCHIP.  It is copied/modified from
% pchiptx.
%
%
% PCHIPTX  Textbook piecewise cubic Hermite interpolation.
%
%  vv = pchiptx(x,v,xx) finds the shape-preserving piecewise cubic
%  interpolant P(x), with P(x(i)) = v(i), and returns vv(ii) = P(xx(ii)).
%
%  See PCHIP, SPLINETX.
%
%   INPUTS: 
%       x  = size [Nx,1]
%       v  = size [Nx,1]
%       xx = size [Nxx,1]
%
%   OUTPUTS:
%       vv = size [Nxx,1]
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Size of matrices for memory allocation
Nx  = numel(x);
Nxx = numel(xx);
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------  
% Estimate derivative dv/dx at (x)
vx = pchipslopes(x,v);   
% -------------------------------------------------------------------------  


% -------------------------------------------------------------------------  
% Interpolate at (xx)
vv = zeros(size(xx));   

% for each field point... 
for ii = 1:Nxx
    
    % ----------------------------------------------------------------
    % Find index (i) corresponding to neighboring baseline grid point,
    % such that x(i) <= xx(ii) <= x(i+1)
    i = find(x <= xx(ii) ,1,'last');   
    
    if i == Nx, i = Nx-1; end  % Error check, so i < Nx 
    % ----------------------------------------------------------------

    % Interval length, h
    h = x(i+1)-x(i);

    % Sub-interval position, s
    s = (xx(ii) - x(i))/h;

    % Evaluate Hermite polynomials
    H1 = (1 - 3*s^2 + 2*s^3  );
    H2 = (    3*s^2 - 2*s^3  );
    H3 = (      s   * (s-1)^2) * h;
    H4 = (      s^2 * (s-1)  ) * h;


    % Evaluate interpolant
    vv(ii) =   v(i)   * H1 ...
            +  v(i+1) * H2 ...
            ...
            + vx(i)   * H3 ...  
            + vx(i+1) * H4;   
end
% -------------------------------------------------------------------------  


% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
%  pchipslopes - First derivative slopes for shape-preserving Hermite cubic
%
%  pchipslopes(x,v) computes d(i) = dvdx(x(i)).
%
%  Slopes at interior points
%       delta = diff(v)./diff(x)
%
%       d(i)  = 0 if delta(i-1) and delta(i) have opposites signs or either is zero.
%       d(i)  = weighted harmonic mean of delta(i-1) and delta(i) if they have the same sign.
%
% ASSUMES v = v(x) is 1D data

function d = pchipslopes(x,v)
   h     = diff(x);
   delta = diff(v)./h;
   
   Nx   = length(h)+1;
   d    = zeros(size(h));
   i    = find(sign(delta(1:Nx-2)).*sign(delta(2:Nx-1))>0)+1;
   w1   = 2*h(i)+h(i-1);
   w2   = h(i)+2*h(i-1);
   d(i) = (w1+w2)./(w1./delta(i-1) + w2./delta(i));

%  Slopes at endpoints

   d(1)  = pchipend(h(1),h(2),delta(1),delta(2));
   d(Nx) = pchipend(h(Nx-1),h(Nx-2),delta(Nx-1),delta(Nx-2));

% -------------------------------------------------------------------------

function d = pchipend(h1,h2,del1,del2)
%  Noncentered, shape-preserving, three-point formula.
   d = ((2*h1+h2)*del1 - h1*del2)/(h1+h2);
   if sign(d) ~= sign(del1)
      d = 0;
   elseif (sign(del1)~=sign(del2))&(abs(d)>abs(3*del1))
      d = 3*del1;
   end
