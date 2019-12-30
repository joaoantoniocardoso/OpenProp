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
% ================================================= Horseshoe_intr_110830.m
% Created: 8/30/2011, Brenden Epps
% -------------------------------------------------------------------------
%
% This function computes the interaction horseshoe influence functions in 
% the axial, tangential, and radial directions: UAHIF21, UTHIF21, URHIF21.
%
% UAHIF21(m,n) = influence of n-th horseshoe vortex shed from propulsor component 1 (M1 panels) 
%                      on the m-th control point of                     component 2 (M2 panels). 
% 
% -------------------------------------------------------------------------
% 'rotor-rotor' interactions:
%
% UAHIF21(m,n) is computed at (x,r) locations (XC2,RC2(m)) for m=1:M2,
%               where RC2(m) is not necessarily equal to RC1(m),
%               and XC2 is defined positive downstream of component 1.
%                   XC2 is assumed the same for all control points.
%
% For RC2(m) ~= RC1(m),     UAHIF21(RC2(m),n) is linearly interpolated 
%                      from UAHIF21( RP(i),n) and UAHIF21(RP(i+1),n), where
%                      RP = [RC2(RC2<RV1(1)),RV1(1),RC1,RV1(M1+1),RC2(RC2>RV1(M1+1))];
%
%                   Note: This implementation does not truly handle the contracting wake case yet.
%
% -------------------------------------------------------------------------
% 'rotor-duct' interactions:
%
% UAHIF21(m,n) is computed at (x,r) locations (XC2(m),RC2) for m=1:M2,
%               where RC2 is assumed constant for all control points.
%
%               No radial interpolation necessary, since RC2 > RV1(end) is assumed
%               No hub or duct images used in rotor-duct calculations
%               UTHIF21 is zero since RC2 > RV1(end) for duct
%
% -------------------------------------------------------------------------
% Locally-constant-pitch wake assumption is made:
%
%             |
%             |                  * (XC2(m),RP1(i+1))
%             |
%             *---------------       
%             |                  * (XC2(m),RC2(m))   
%             |                  * (XC2(m),RP1(i))                       
%             |
%             *---------------
%             |
%             |
%             |        
%             *--------------- r1 = RV1(n+1);  TANBIV1 = TANBIC1(n)*RC1(n)/RV1(n+1) 
%             |
%             |  <--- RC1(n), TANBIC1(n), G1(n)
%             |
%             *--------------- r2 = RV1(n  );  TANBIV2 = TANBIC1(n)*RC1(n)/RV1(n)
%             |
%             |
%             |
%
%
% -------------------------------------------------------------------------
%
% Since, the axial and radial influence functions are proportional to 
% 1/TANBIC1, one can call...
%
%       [UAHIF21_times_TANBIC1, UTHIF21, URHIF21_times_TANBIC1] = Horseshoe_intr_110830(...,TANBIC1 == 1,...)
%
% ...and then complete the calculation by:
%
%     for n = 1:M1;                               
%         UAHIF21(:,n) = UAHIF21_times_TANBIC1(:,n) / TANBIC1(n);
%         URHIF21(:,n) = URHIF21_times_TANBIC1(:,n) / TANBIC1(n);
%     end   
%
% -------------------------------------------------------------------------
%
% Compute_flag = [Compute_UA, Compute_UT, Compute_UR] is given as input.
%                tells the code which influence function(s) to compute
%
% e.g. For the influence of one rotor on another, use [1 1 0], scalar XC2, and array RC2(n).
%      For the influence of a   rotor on a duct,  use [1 0 1], scalar RC2, and array XC2(n).
%
% Note: For a duct, use RC2 >= 1.05, since URHIF21 is nearly singular at (XC2=0,RC2=1). 
%                   With M1 == 10 - 40,    URHIF21 is accurate for all XC2 only for RC2 >= 1.05
%
% -------------------------------------------------------------------------
% Non-dimensionalization (given R1, Vs) 
%
%   UASTAR2 == UAHIF21 * G1    <----->  uastar2/Vs == { (2*pi*R1) * uahif21 } * {Gamma1 / (2*pi*R1*Vs)}
%
%   UASTAR2 == uastar2/Vs
%   UAHIF21 ==          (2*pi*R1) * uahif21   
%        G1 == Gamma1 / (2*pi*R1*Vs)
%   All lengths normalized by R1
% -------------------------------------------------------------------------


function [UAHIF21,UTHIF21,URHIF21] = Horseshoe_intr_110830(XC2,RC2,RC1,TANBIC1,RV1,Z1,Hub_flag,Rhub_oR1,Duct_flag,Rduct_oR1) 

% ------------- Error checking, must have control points given as (XC2,RC2)
if        length(RC2) > 1  &  length(XC2) == 1         % (rotor-rotor case)
                                      Interaction_type = 'rotor-rotor';    
     M2 = length(RC2);
     
    XC2 = XC2 * ones(size(RC2));
     
    Compute_flag = [1 1 0];
    
elseif    length(XC2) > 1  &  length(RC2) == 1         % (rotor-duct  case)
                                      Interaction_type = 'rotor-duct';                  
     M2 = length(XC2);
    
    if RC2 < 1.05, RC2 = 1.05; end  % prevent singularity at (XC2=0,RC2=1)
    
       RC2 = RC2 * ones(size(XC2));
    
    Compute_flag = [1 0 1];
    
elseif length(XC2) ~= length(RC2)
    disp('(Horseshoe_intr_110830) ERROR: Must have equal number of given radii and axial control point locations.')
    return
    
else  % length(XC2) == length(RC2) 
    Interaction_type = 'rotor-rotor'; 
    Compute_flag     = [1 1 1];
    M2 = length(XC2);
end
% ------------------------------------------------------------------------- 

M1 = length(RC1);

UAHIF21 = zeros(M2,M1);
UTHIF21 = zeros(M2,M1);
URHIF21 = zeros(M2,M1);

Compute_UA = Compute_flag(1);
Compute_UT = Compute_flag(2);
Compute_UR = Compute_flag(3);


% -------------------------------------------------------------------------
if strcmp(Interaction_type,'rotor-rotor')  
                
    % ---------------------------------------------------------------------
    % Radial interpolation is necessary in the case that RC2 ~= RC1
    % ---------------------------------------------------------------------
    % RP = radii to evaluate circumferential average velocities
    RP =          [RC2(RC2<RV1(1)),RV1(1),RC1,RV1(M1+1),RC2(RC2>RV1(M1+1))];
    Mfor = length([RC2(RC2<RV1(1)),RV1(1)]);
    Maft = length(                           [RV1(M1+1),RC2(RC2>RV1(M1+1))]);

    Mp = length(RP);  % == Mfor + M1 + Maft

    UAHIFtemp = zeros(Mp,M1);
    UTHIFtemp = zeros(Mp,M1);
    URHIFtemp = zeros(Mp,M1);

    % ----------------------------- Find axial horseshoe influence function
    if Compute_UA == 1                      
        for m = 1:Mp                     % for each interpolation radius, m
            for n = 1:M1                 % for each vortex  panel, n

               % Assume all XC2 == XC2(1) 

                      UA1 = HoughUA(XC2(1),RP(m), TANBIC1(n)*RC1(n)/RV1(n+1) ,RV1(n+1),Z1);  % Velocity induced at RP(m) by a unit vortex shed at RV1(n+1)
                      UA2 = HoughUA(XC2(1),RP(m), TANBIC1(n)*RC1(n)/RV1(n)   ,RV1(n)  ,Z1);  % Velocity induced at RP(m) by a unit vortex shed at RV1(n)  

                if Hub_flag == 1 

                    UA1_h = HoughUA(XC2(1),RP(m), TANBIC1(n)*RC1(n)/(Rhub_oR1^2/RV1(n+1)) , Rhub_oR1^2/RV1(n+1), Z1 ); 
                    UA2_h = HoughUA(XC2(1),RP(m), TANBIC1(n)*RC1(n)/(Rhub_oR1^2/RV1(n  )) , Rhub_oR1^2/RV1(n  ), Z1 ); 

                    UA1 = UA1 - UA1_h;
                    UA2 = UA2 - UA2_h;
                end

                if Duct_flag == 1 

                    UA1_d = HoughUA(XC2(1),RP(m), TANBIC1(n)*RC1(n)/(Rduct_oR1^2/RV1(n+1)) , Rduct_oR1^2/RV1(n+1) , Z1);
                    UA2_d = HoughUA(XC2(1),RP(m), TANBIC1(n)*RC1(n)/(Rduct_oR1^2/RV1(n  )) , Rduct_oR1^2/RV1(n  ) , Z1); 

                    UA1 = UA1 - UA1_d;
                    UA2 = UA2 - UA2_d;
                end


                UAHIFtemp(m,n) = UA1 - UA2;      

            end % for n = 1:M1
        end     % for m = 1:M2  
    end
    % ---------------------------------------------------------------------


    % ---------------------------- Find radial horseshoe influence function 
    if Compute_UR == 1
        for m = 1:Mp                     % for each interpolation radius, m
            for n = 1:M1                 % for each vortex  panel, n

               % Assume all XC2 == XC2(1) 

                      UR1 = HoughUR(XC2(1),RP(m), TANBIC1(n)*RC1(n)/RV1(n+1) ,RV1(n+1),Z1);  % Velocity induced at RP(m) by a unit vortex shed at RV1(n+1)
                      UR2 = HoughUR(XC2(1),RP(m), TANBIC1(n)*RC1(n)/RV1(n)   ,RV1(n)  ,Z1);  % Velocity induced at RP(m) by a unit vortex shed at RV1(n)  

                if Hub_flag == 1 

                    UR1_h = HoughUR(XC2(1),RP(m), TANBIC1(n)*RC1(n)/(Rhub_oR1^2/RV1(n+1)) , Rhub_oR1^2/RV1(n+1), Z1 ); 
                    UR2_h = HoughUR(XC2(1),RP(m), TANBIC1(n)*RC1(n)/(Rhub_oR1^2/RV1(n  )) , Rhub_oR1^2/RV1(n  ), Z1 ); 

                    UR1 = UR1 - UR1_h;
                    UR2 = UR2 - UR2_h;
                end

                if Duct_flag == 1 

                    UR1_d = HoughUR(XC2(1),RP(m), TANBIC1(n)*RC1(n)/(Rduct_oR1^2/RV1(n+1)) , Rduct_oR1^2/RV1(n+1) , Z1);
                    UR2_d = HoughUR(XC2(1),RP(m), TANBIC1(n)*RC1(n)/(Rduct_oR1^2/RV1(n  )) , Rduct_oR1^2/RV1(n  ) , Z1); 

                    UR1 = UR1 - UR1_d;
                    UR2 = UR2 - UR2_d;
                end


                URHIFtemp(m,n) = UR1 - UR2;      

            end % for n = 1:M1
        end     % for m = 1:M2
    end
    % ---------------------------------------------------------------------


    % ------------------------ Find tangential horseshoe influence function
    %  (Note the simple equation, since no effect of hub or duct images)   
    if Compute_UT == 1        
        for n = 1:M1                   % for each vortex panel mid-point, n

            % Index offset, since:  RP(Mfor+n) == RC1(n)
            i = Mfor + n;

            if      XC2(1) > 0

                UTHIFtemp(i,n) = Z1/RC1(n);

            elseif  XC2(1) == 0

                UTHIFtemp(i,n) = Z1/(2*RC1(n));

            else
                UTHIFtemp(i,n) = 0;
            end         

        end % for n = 1:M1  
    end
    % ---------------------------------------------------------------------


    % ---------------------------------------- Perform linear interpolation
    for m = 1:M2                                % for each control point, m

        i = find(RP > RC2(m), 1)-1;
            % error checking
                if i == 0,           i = 1;           
            elseif i == Mp,          i = Mp-1;  
            elseif isempty(i),       i = Mp-1;  
               end        

        % Percentage of distance across the element        
        s = (RC2(m) - RP(i))/(RP(i+1) - RP(i));

        % Linear interpolation
        for n = 1:M1
            UAHIF21(m,n) = (1-s)*UAHIFtemp(i,n) + s*UAHIFtemp(i+1,n);
            UTHIF21(m,n) = (1-s)*UTHIFtemp(i,n) + s*UTHIFtemp(i+1,n);
            URHIF21(m,n) = (1-s)*URHIFtemp(i,n) + s*URHIFtemp(i+1,n);
        end
    end
    % -----------------------------------------------------------------                                   
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
if strcmp(Interaction_type,'rotor-duct')

    % No radial interpolation necessary, since RC2 > RV1(end) is assumed
    % No hub or duct images used in rotor-duct calculations
    % UTHIF21 is zero since RC2 > RV1(end) for duct

    for m = 1:M2                                % for each control point, m
        for n = 1:M1                            % for each vortex  panel, n

            % -------------------------------------------------------------
            % Find axial horseshoe influence function

            UA1 = HoughUA(XC2(m),RC2(m), TANBIC1(n)*RC1(n)/RV1(n+1) ,RV1(n+1),Z1);  % Velocity induced at RC2(m) by a unit vortex shed at RV1(n+1)
            UA2 = HoughUA(XC2(m),RC2(m), TANBIC1(n)*RC1(n)/RV1(n)   ,RV1(n)  ,Z1);  % Velocity induced at RC2(m) by a unit vortex shed at RV1(n)  

            UAHIF21(m,n) = UA1 - UA2;  
            % -------------------------------------------------------------

            % -------------------------------------------------------------
            % Find radial horseshoe influence function

            UR1 = HoughUR(XC2(m),RC2(m), TANBIC1(n)*RC1(n)/RV1(n+1) ,RV1(n+1),Z1);  % Velocity induced at RC2(m) by a unit vortex shed at RV1(n+1)
            UR2 = HoughUR(XC2(m),RC2(m), TANBIC1(n)*RC1(n)/RV1(n)   ,RV1(n)  ,Z1);  % Velocity induced at RC2(m) by a unit vortex shed at RV1(n)  

            URHIF21(m,n) = UR1 - UR2;  
            % -------------------------------------------------------------
        end
    end
end
% -------------------------------------------------------------------------

end % function [UAHIF21,UTHIF21] = Horseshoe_int(...)
% =========================================================================



% =========================================================================
% =========================================================================
% HoughUA(...) function
%
% UA = (circumferrential mean axial induced velocity) * (2*pi*R)
%
% Velocity induced at control point (xc,rc) by a set of Z unit-strength
% traling vortices shed at (x=0,rv) with pitch tanbiv.
%
% xc = distance downstream of lifting line (negative for upstream)
%
% xc,rc,rv normalized by radius R
%
% NOTE: With epsilon=0, HoughUA evaluates to NaN at (xc=0, rc=rv), so this point must be avoided
%
function UA = HoughUA(xc,rc, tanbiv,rv,Z)

% Legendre function
    epsilon = 0.025;  % small number to prevent singularity at (xc=0, rc=rv)

    % z == argument for Legendre functions 
      z =  1+(xc^2+(rc-rv)^2 + epsilon^2)/(2*rc*rv);

    Q_Mhalf = Legendre_2nd_kind_minus_half_order(z);

  
% Heuman Lambda Function:
    % s ==    amplitude in elliptic integrals
    % t == k == modulus in elliptic integrals    
      s = asin(   xc   /sqrt(xc^2+(rc-rv)^2) );
      t = sqrt( 4*rc*rv/    (xc^2+(rc+rv)^2) );  

      Heuman_Lambda = Heuman(s,t);                                  
  
  
% K1 == constant in UA integral
    if (rc <= rv && xc < 0) || (rc < rv && xc >= 0)     
        K1 = pi + xc/(2*sqrt(rc*rv))*Q_Mhalf + pi/2*Heuman_Lambda;
    else
        K1 =      xc/(2*sqrt(rc*rv))*Q_Mhalf - pi/2*Heuman_Lambda;
    end
    
    K1(isnan(K1)) = 0;
   
    
% Axial induced velocity:
UA = Z*K1/(2*pi*rv*tanbiv);

end
% =========================================================================    


% =========================================================================
% HoughUR(...) function
%
% UR = (circumferrential mean radial induced velocity) * (2*pi*R)
%
% Velocity induced at control point (xc,rc) by a set of Z unit-strength
% traling vortices shed at (x=0,rv) with pitch tanbiv.
%
% xc = distance downstream of lifting line (negative for upstream)
%
% xc,rc,rv normalized by radius R
%
% NOTE: With epsilon=0, HoughUR evaluates to Inf at (xc=0, rc=rv), so this point must be avoided
%
function UR = HoughUR(xc,rc, tanbiv,rv,Z)

% Legendre function
    epsilon = 0.025;  % small number to prevent singularity at (xc=0, rc=rv)

    % z == argument for Legendre functions 
      z =  1+(xc^2+(rc-rv)^2 + epsilon^2)/(2*rc*rv);

    Q_half = Legendre_2nd_kind_half_order(z);


% Radial induced velocity
UR = -Z/(2*pi*sqrt(rc)) * (1/(rv*tanbiv)) * sqrt(rv) * Q_half;

end
% =========================================================================


% =========================================================================
% =========================================================================
% ========================================================== math Functions


% =========================================================================
% -------------------------------------------------------------------------
% Heuman: Heuman's Lambda fuction
%
% Reference: 
%   1) Abramowitz and Stegun, Handbook of Math Functions, 
%      section 17.4.39, p.595, 1972
%
%   2) Byrd and Friedman, Handbook of Elliptic Integrals for Engineers 
%      and Physicists, p37, 1954.
%
% -------------------------------------------------------------------------
%
% phi:       amplitude (radians),        (ductCMV sends 's' as phi)
% k:         modulus,                    (ductCMV sends 't' as k)
%
% Note: the "modular angle (radians)", alpha, is defined from 
%       the "modulus", k, by equation
%
%               alpha = asin(k)
%
% Abramowitz and Stegun define Heuman in terms of phi and alpha.
%
% -------------------------------------------------------------------------

function [H] = Heuman(phi,k)

% Complete elliptic integals:
[K,E] = elliptic12(pi/2,k^2);

alpha = asin(k);  % modular angle (radians)

% Incomplete elliptic integrals:
%   F  = 1st kind
%   EE = 2nd kind
                                    
[F,EE] = elliptic12(phi, sin(pi/2-alpha).^2);  % values of integrals

H     = 2/pi * (K*EE - (K-E)*F);

end % function [H] = Heuman(phi,k)
% -------------------------------------------------------------------------



% =========================================================================
% -------------------------------------------------------------------------
% Legendre_2nd_kind_half_order.m
%
% Legendre fuction of the second kind and positive half order
%
% Verified: Brenden Epps, 8/12/2011
% 
% Usage:
%           [Q] = Legendre_2nd_kind_half_order(z)
%   where:
%           z > 1
%
% Reference: 
%
%   1) Abramowitz and Stegun, Handbook of Math Functions, 1972, p.337
%
%           Equation 8.13.7 uses modulus k for elliptic integrals **
%           elliptic12 uses parameter m, where m = k^2.
%
%   2) **See also, NIST Handbook of Mathematical Functions, 2010, p. 360.
% -------------------------------------------------------------------------

% function [Q_half, dQ_half_dz] = Legendre_2nd_kind_half_order(z)
  function  Q_half              = Legendre_2nd_kind_half_order(z)

k      = sqrt(2./(z+1));

[K,E]  = elliptic12(pi/2,k.^2);

Q_half = z .* k .* K - sqrt(2*(z+1)) .* E;


% % ------- Derivative with respect to z
% dk_dz  = sqrt(2) * (-1/2) * (z+1).^(-3/2);
% dK_dk  = (E - (1-k.^2).*K)./(k.*(1-k.^2));
% dE_dk  = (E - K) ./ k;
% 
% dQ_half_dz = k .* K  +  z .* dk_dz .* K  + z .* k .* dK_dk .* dk_dz ...
%              + (2./k.^2) .* dk_dz .* E - (2./k) .* dE_dk .* dk_dz;
% % ------- 

end 
% =========================================================================



% =========================================================================
% -------------------------------------------------------------------------
% Legendre_2nd_kind_minus_half_order.m
%
% Legendre fuction of the second kind and positive half order
%
% Verified: Brenden Epps, 8/12/2011
% 
% Usage:
%           [Q] = Legendre_2nd_kind_half_order(z)
%   where:
%           z > 1
%
% Reference: 
%
%   1) Abramowitz and Stegun, Handbook of Math Functions, 1972, p.337
%
%           Equation 8.13.3 uses modulus k for elliptic integrals **
%           elliptic12 uses parameter m, where m = k^2.
%
%   2) **See also, NIST Handbook of Mathematical Functions, 2010, p. 360.
% -------------------------------------------------------------------------

% function [Q_Mhalf, dQ_Mhalf_dz] = Legendre_2nd_kind_minus_half_order(z)
  function  Q_Mhalf               = Legendre_2nd_kind_minus_half_order(z)

k      = sqrt(2./(z+1));

[K,E]  = elliptic12(pi/2,k.^2);

Q_Mhalf = k .* K;


% % ------- Derivative with respect to z
% dk_dz  = sqrt(2) * (-1/2) * (z+1).^(-3/2);
% dK_dk  = (E - (1-k.^2).*K)./(k.*(1-k.^2));
% dE_dk  = (E - K) ./ k;
% 
% dQ_Mhalf_dz = dk_dz .* K  + k .* dK_dk .* dk_dz;
% % ------- 

end 
% =========================================================================

 

% =========================================================================
% function [F,E,Z] = elliptic12(u,m,tol)
  function [F,E,Z] = elliptic12(u,m)
% ELLIPTIC12 evaluates the value of the Incomplete Elliptic Integrals
% of the First, Second Kind and Jacobi's Zeta Function.
%
%   [F,E,Z] = ELLIPTIC12(U,M,TOL) where U is a phase in radians, 0<M<1 is
%   the module and TOL is the tolerance (optional). Default value for
%   the tolerance is eps = 2.220e-16.
%
%   ELLIPTIC12 uses the method of the Arithmetic-Geometric Mean
%   and Descending Landen Transformation described in [1] Ch. 17.6,
%   to determine the value of the Incomplete Elliptic Integrals
%   of the First, Second Kind and Jacobi's Zeta Function [1], [2].
%
%       F(phi,m) = int(1/sqrt(1-m*sin(t)^2), t=0..phi);
%       E(phi,m) = int(sqrt(1-m*sin(t)^2), t=0..phi);
%       Z(phi,m) = E(u,m) - E(m)/K(m)*F(phi,m).
%
%   Tables generating code ([1], pp. 613-621):
%       [phi,alpha] = meshgrid(0:5:90, 0:2:90);                  % modulus and phase in degrees
%       [F,E,Z] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);  % values of integrals
%
%   See also ELLIPKE, ELLIPJ, ELLIPTIC12I, ELLIPTIC3, THETA, AGM.
%
%   References:
%   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
%       Dover Publications", 1965, Ch. 17.1 - 17.6 (by L.M. Milne-Thomson).
%   [2] D. F. Lawden, "Elliptic Functions and Applications"
%       Springer-Verlag, vol. 80, 1989

% GNU GENERAL PUBLIC LICENSE Version 2, June 1991
% http://www.gnu.org/licenses/gpl.html
% Everyone is permitted to copy and distribute verbatim copies of this
% script under terms and conditions of GNU GENERAL PUBLIC LICENSE.
%  
% Copyright (C) 2007 by Moiseev Igor. All rights reserved.
% 34106, SISSA, via Beirut n. 2-4,  Trieste, Italy
% For support, please reply to
%     moiseev[at]sissa.it, moiseev.igor[at]gmail.com
%     Moiseev Igor,
%     34106, SISSA, via Beirut n. 2-4,  Trieste, Italy
%
% The code is optimized for ordered inputs produced by the functions
% meshgrid, ndgrid. To obtain maximum performace (up to 30%) for singleton,
% 1-dimensional and random arrays remark call of the function unique(.)
% and edit further code.

% % Input error checking:
% if nargin<3, tol = eps; end
% if nargin<2, error('Not enough input arguments.'); end
% 
% if ~isreal(u) || ~isreal(m)
%     error('Input arguments must be real. Use ELLIPTIC12i for complex arguments.');
% end
% 
% if length(m)==1, m = m(ones(size(u))); end
if length(u)==1, u = u(ones(size(m))); end
% if ~isequal(size(m),size(u)), error('U and M must be the same size.'); end
tol = eps;

F = zeros(size(u));
E = F;              
Z = E;
m = m(:).';    % make a row vector
u = u(:).';

if any(m < 0) || any(m > 1), error('M must be in the range 0 <= M <= 1.'); end

I = uint32( find(m ~= 1 & m ~= 0) );
if ~isempty(I)
    [mu,J,K] = unique(m(I));   % extracts unique values from m
    K = uint32(K);
    mumax = length(mu);
    signU = sign(u(I));

    % pre-allocate space and augment if needed
        chunk = 7;
        a = zeros(chunk,mumax);
        c = a;
        b = a;
        a(1,:) = ones(1,mumax);
        c(1,:) = sqrt(mu);
        b(1,:) = sqrt(1-mu);
        n = uint32( zeros(1,mumax) );
        i = 1;
        while any(abs(c(i,:)) > tol)                                    % Arithmetic-Geometric Mean of A, B and C
        i = i + 1;
        if i > size(a,1)
          a = [a; zeros(2,mumax)];
          b = [b; zeros(2,mumax)];
          c = [c; zeros(2,mumax)];
        end
        a(i,:) = 0.5 * (a(i-1,:) + b(i-1,:));
        b(i,:) = sqrt(a(i-1,:) .* b(i-1,:));
        c(i,:) = 0.5 * (a(i-1,:) - b(i-1,:));
        in = uint32( find((abs(c(i,:)) <= tol) & (abs(c(i-1,:)) > tol)) );
        if ~isempty(in)
          [mi,ni] = size(in);
          n(in) = ones(mi,ni)*(i-1);
        end
        end
     
    mmax = length(I);
        mn = double(max(n));
        phin = zeros(1,mmax);     C  = zeros(1,mmax);    
        Cp = C;  e  = uint32(C);  phin(:) = signU.*u(I);
        i = 0;   c2 = c.^2;
        while i < mn                                                    % Descending Landen Transformation
        i = i + 1;
        in = uint32(find(n(K) > i));
        if ~isempty(in)    
            phin(in) = atan(b(i,K(in))./a(i,K(in)).*tan(phin(in))) + ...
                pi.*ceil(phin(in)/pi - 0.5) + phin(in);
            e(in) = 2.^(i-1) ;
            C(in) = C(in)  + double(e(in(1)))*c2(i,K(in));
            Cp(in)= Cp(in) + c(i+1,K(in)).*sin(phin(in));  
        end
        end
   
    Ff = phin ./ (a(mn,K).*double(e)*2);                                                      
    F(I) = Ff.*signU;                                               % Incomplete Ell. Int. of the First Kind
    Z(I) = Cp.*signU;                                               % Jacobi Zeta Function
    E(I) = (Cp + (1 - 1/2*C) .* Ff).*signU;                         % Incomplete Ell. Int. of the Second Kind
end

% Special cases: m == {0, 1}
m0 = find(m == 0);
if ~isempty(m0), F(m0) = u(m0); E(m0) = u(m0); Z(m0) = 0; end

m1 = find(m == 1);
um1 = abs(u(m1));
if ~isempty(m1),
    N = floor( (um1+pi/2)/pi );  
    M = find(um1 < pi/2);              
   
    F(m1(M)) = log(tan(pi/4 + u(m1(M))/2));  
    F(m1(um1 >= pi/2)) = Inf.*sign(u(m1(um1 >= pi/2)));
   
    E(m1) = ((-1).^N .* sin(um1) + 2*N).*sign(u(m1));
   
    Z(m1) = (-1).^N .* sin(u(m1));                      
end
end % function [F,E,Z] = elliptic12(u,m)
% =========================================================================

% ====================================================== END math Functions
% =========================================================================