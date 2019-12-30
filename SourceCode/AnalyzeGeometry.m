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
% Last modified: 11/2/2011 Brenden Epps
% =========================================================================

% =========================================================================
% ======================================================== Analyze Function
%
% This function computes the "states" of a given propeller/turbine "design"
% at given off-design tip speed ratios, "LAMBDAall".
%
% -------------------------------------------------------------------------
%
% Inputs:
%       ginput      propeller/turbine input geometry data structure
% 
%       LAMBDAall   tip-speed ratios, omega*R/Vs
%
%       LSGeoCorr   flag for lifting surface geometry correction. Options:
%                       'none'               (turbine   default)
%                       'Morgan1968'         (propeller default)
%                       'EckhardtMorgan1955'
%                       'vanManen1958'
% 
%       feather     [deg], feathering angle, i.e. pitch angle displacement
%                          relative to design geometry, such that pitch 
%                          angle with feathering == pitch angle + feather
%                          (for controllable-pitch rotors)
% Outputs:   
%       states      data structure with operating states at LAMBDAall
%
% -------------------------------------------------------------------------
% Method:
%       This is a continuation of OpenProp v3.1.0 AnalyzeGeometry110628.m
%
%       Newton solver state variable: (VSTAR,ALPHA,CL,G,UASTAR,UTSTAR), size [6,1] 
%       with (TANBIC,...) updated between iterations
%
%       dCLdALPHA is estimated from elliptic wing unless specified
%
% -------------------------------------------------------------------------


function [states,ginput] = AnalyzeGeometry(ginput,LAMBDAall,feather)

if nargin == 2,
    feather    = 0;    
end

% =========================================================================
feather = feather * pi/180;  % [rad] feathering angle
% =========================================================================

warning off


% ========================================================== Error checking
% 0 == Horseshoe(...,Wrench(...)), 1 == Wake_Horseshoe(...)      
if isfield(ginput,'Wake_flag'),  Wake_flag = ginput.Wake_flag;  
                          else   Wake_flag = 0;  end
                          
if Wake_flag == 1, disp('ERROR: Wake_Horseshoe(...) is no longer supported in AnalyzeGeometry110628.'), return, end


if isfield(ginput,'Duct_flag'),  Duct_flag = ginput.Duct_flag;   % 0 == no duct, 1 == duct
                          else   Duct_flag = 0;  end

if Duct_flag == 1, disp('ERROR: AnalyzeGeometry110628 can not model the ducted case.'), return, end
% =========================================================================
%%


              

% -------------------------------------------------------------------------
% If raw blade geometry (XR,XCoD,...) given, create vortex lattice geometry
% -------------------------------------------------------------------------
if ~isfield(ginput,'Rhub_oR')    
    if     isfield(ginput,'Rhub') && isfield(ginput,'R'),  ginput.Rhub_oR = ginput.Rhub/ginput.R;
    elseif isfield(ginput,'Dhub') && isfield(ginput,'D'),  ginput.Rhub_oR = ginput.Dhub/ginput.D;
    elseif isfield(ginput,'Dhub_oD'),                      ginput.Rhub_oR = ginput.Dhub_oD;      
    else
        disp('ERROR: Must specify hub diameter.')
        return
    end
end

% 0 == no hub, 1 == hub
if isfield(ginput,'Hub_flag'),    Hub_flag = ginput.Hub_flag;  
                           else   Hub_flag = 0;  end
                           

% 0 == use CLCD_vs_ALPHA.m , 1 == use CL and CD data given in ginput
if isfield(ginput,'CLCD_flag'),  CLCD_flag = ginput.CLCD_flag;  
                          else   CLCD_flag = 0;  end
                          
                          
% =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=                           
if isfield(ginput,'XR') && ~isfield(ginput,'RC')

    ginput.part2 = '----------- Vortex lattice geometry --------';
    
    if ~isfield(ginput,'Mp'),  ginput.Mp = 20;   end


    [ginput.RC,ginput.RV,ginput.DR] = LLPanelRadii(ginput.Mp,ginput.Rhub_oR,Hub_flag,Duct_flag);

    RC = ginput.RC;
    
    % --------------------------------------- Interpolate to control points
    if isfield(ginput,'XVA')
        ginput.VAC  = pchip(ginput.XR,ginput.XVA   ,RC);   % axial      inflow vel. / ship vel. at ctrl pts
    else
        ginput.VAC  = ones(size(RC));
    end
    
    if isfield(ginput,'XVT')
        ginput.VTC  = pchip(ginput.XR,ginput.XVT   ,RC);   % tangential inflow vel. / ship vel. at ctrl pts
    else
        ginput.VTC  = zeros(size(RC));
    end
    
    if isfield(ginput,'XPoD')
        ginput.PoD      = pchip(ginput.XR,ginput.XPoD     ,RC);
    else
        ginput.THETArad = pchip(ginput.XR,ginput.XTHETArad,RC);
    end

    ginput.CoD  = InterpolateChord(ginput.XR,ginput.XCoD, RC);   % section chord / propeller diameter at ctrl pts

    
    if CLCD_flag == 0  % assume geometry is given as Meanline and Thickness distributions
    
        ginput.CD   = pchip(ginput.XR,ginput.XCD   ,RC);   % section drag coefficient           at ctrl pts
        ginput.f0oc = pchip(ginput.XR,ginput.f0oc0 ,RC);
    
        if isfield(ginput,'t0oc0' ), 
            ginput.t0oc  = pchip(ginput.XR,ginput.t0oc0 ,RC);
        else
            ginput.t0oc  = pchip(ginput.XR,ginput.Xt0oD ,RC) ./ max([CoD,1e-8]);  
        end 

    
        XXR = linspace(ginput.Rhub_oR,1,100); % temp radii

        ginput.EAR = (2*ginput.Z/pi) * trapz( XXR, InterpolateChord(ginput.XR,ginput.XCoD, XXR)  ); % expanded area ratio 

        ginput.BTF = interp1(ginput.RC,ginput.t0oc.*ginput.CoD,0,'linear','extrap');                % blade thickness fraction

               XVAtemp  = pchip(ginput.XR,ginput.XVA,XXR);          
        ginput.VMIV     = 2*trapz(XXR,XXR.*XVAtemp)/(1-ginput.Rhub_oR^2);  % Volumetric Mean Inflow Velocity

        % ----------------------------------------------------------- dCLdALPHA
        if isfield(ginput,'dCLdALPHA' ), 

            dCLdALPHA  = ginput.dCLdALPHA;

            % make dCLdALPHA the same length as RC...
            if      length(dCLdALPHA) == 1
                    ginput.dCLdALPHA  = dCLdALPHA * ones(size(RC));        
            else
                    ginput.dCLdALPHA  = pchip(ginput.XR,dCLdALPHA,RC);
            end
        else
            % Propeller Aspect Ratio (PAR)
            PAR = (1-ginput.Rhub_oR)^2 / trapz( XXR , InterpolateChord(ginput.XR,ginput.XCoD,XXR) );

            % Lift curve slope:
            ginput.dCLdALPHA = (2*pi / (1 + 2/PAR)) * ones(size(RC));
        end
        % --------------------------------------------------------------------- 

        % --------------------------------------------------- 2D section performance
        if isfield(ginput,'f0octilde') && isfield(ginput,'CLItilde') && isfield(ginput,'alphaItilde')

            if length(ginput.f0octilde) == length(ginput.XR)        
                % Then we need to interpolate onto RC
                ginput.f0octilde   = pchip(ginput.XR,ginput.f0octilde,   RC);  % versus RC
                ginput.CLItilde    = pchip(ginput.XR,ginput.CLItilde,    RC);  % versus RC
                ginput.alphaItilde = pchip(ginput.XR,ginput.alphaItilde, RC);  % versus RC         
            end

        else % Assume Meanline and Thickness are given, possibly versus XR

            Meanline  = ginput.Meanline;
            Thickness = ginput.Thickness;
            XR        = ginput.XR;

            if iscell(Meanline)

                if length(Meanline) ~= length(XR)
                    disp('<ERROR>')
                    disp('<ERROR> Meanline given as cell array but different length as XR...crash.')
                    disp('<ERROR>')
                    return
                else
                      f0octilde = 0*XR; % memory allocation
                       CLItilde = 0*XR; % memory allocation
                    alphaItilde = 0*XR; % memory allocation
                    
                    for j = 1:length(XR)

                        [f0octilde(j), CLItilde(j), alphaItilde(j)] = GeometryFoil2D(Meanline{j},Thickness{j});
                    end
                end

                % Finally we need to interpolate onto RC
                ginput.f0octilde   = pchip(ginput.XR,  f0octilde,    RC);  % versus RC
                ginput.CLItilde    = pchip(ginput.XR,   CLItilde,    RC);  % versus RC
                ginput.alphaItilde = pchip(ginput.XR,alphaItilde,    RC);  % versus RC
            else
                [ginput.f0octilde, ginput.CLItilde, ginput.alphaItilde] = GeometryFoil2D(Meanline,Thickness);           
            end
        end
        % -------------------------------------------------------------------------
    else 
        ginput.VMIV = 1;
    end
% =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=    
else  % Assume all parameters are given versus RC
    
    RC        = ginput.RC;
    ginput.Mp = length(RC);
        
    if CLCD_flag == 0  % assume geometry is given as Meanline and Thickness distributions
    
        if ~isfield(ginput,'t0oc' ), 
            ginput.t0oc = ginput.t0oD ./ max([ginput.CoD,1e-8]);  
        end 


        XXR = linspace(ginput.Rhub_oR,1,100); % temp radii


        ginput.EAR = (2*ginput.Z/pi) * trapz( XXR, InterpolateChord(ginput.RC,ginput.CoD, XXR)  ); % expanded area ratio 

        ginput.BTF = interp1(ginput.RC,ginput.t0oc.*ginput.CoD,0,'linear','extrap');                % blade thickness fraction

               XVAtemp  = pchip(ginput.RC,ginput.VAC,XXR);          
        ginput.VMIV     = 2*trapz(XXR,XXR.*XVAtemp)/(1-ginput.Rhub_oR^2);  % Volumetric Mean Inflow Velocity
        
        
        % ----------------------------------------------------------- dCLdALPHA
        if isfield(ginput,'dCLdALPHA' ), 

            dCLdALPHA  = ginput.dCLdALPHA;

            % make dCLdALPHA the same length as RC...
            if      length(dCLdALPHA) == 1
                    ginput.dCLdALPHA  = dCLdALPHA * ones(size(ginput.RC));      
            end

        else    
            XXR = linspace(ginput.Rhub_oR,1,100); % temp radii

            % Propeller Aspect Ratio (PAR)
            PAR = (1-ginput.Rhub_oR)^2 / trapz( XXR , InterpolateChord(ginput.RC,ginput.CoD,XXR)  );     

            % Lift curve slope:
            ginput.dCLdALPHA = (2*pi/(1 + 2/PAR)) * ones(size(ginput.RC));
        end
        % ------------------------------------------------------------------------- 


        % --------------------------------------------------- 2D section performance
        if isfield(ginput,'f0octilde') && isfield(ginput,'CLItilde') && isfield(ginput,'alphaItilde')

             % Then we are all set!

        else % Assume Meanline and Thickness are given, possibly versus RC

            Meanline  = ginput.Meanline;
            Thickness = ginput.Thickness;

            if iscell(Meanline)

                if length(Meanline) ~= length(RC)
                    disp('<ERROR>')
                    disp('<ERROR> Meanline given as cell array but different length as RC...crash.')
                    disp('<ERROR>')
                    return
                else
                      f0octilde = 0*RC; % memory allocation
                       CLItilde = 0*RC; % memory allocation
                    alphaItilde = 0*RC; % memory allocation
                    
                    for j = 1:length(RC)

                        [f0octilde(j), CLItilde(j), alphaItilde(j)] = GeometryFoil2D(Meanline{j},Thickness{j});  % versus RC
                    end
                end
            else
                [f0octilde, CLItilde, alphaItilde] = GeometryFoil2D(Meanline,Thickness);  % scalar
            end

            ginput.f0octilde   =   f0octilde;  % either scalar or versus RC
            ginput.CLItilde    =    CLItilde;  % either scalar or versus RC
            ginput.alphaItilde = alphaItilde;  % either scalar or versus RC
        end
        % -------------------------------------------------------------------------
    else 
        VMIV = 1;
    end
    
end
% =+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------









%%
% =========================================================================
% -------------------------- Initialize output data structure "states.____"
% Sort LAMBDAall into asending order
LAMBDAall = sort(LAMBDAall);

NL = length(LAMBDAall);
Mp = ginput.Mp;


states.part1   = '------ Off-design states ------';
states.L       = LAMBDAall(:);    % [NL x 1] tip speed ratio 
states.Js      = pi./states.L;    % [NL x 1] advance coefficient 

states.RC      = ginput.RC;       % [ 1 x Mp] radius / propeller radius
states.CoD     = ginput.CoD;      % [ 1 x Mp] chord  / propeller diameter 

states.G       = zeros(NL,Mp);    % [NL x Mp] circulation distribution
states.UASTAR  = zeros(NL,Mp);    % [NL x Mp] axial      induced velocity distribution
states.UTSTAR  = zeros(NL,Mp);    % [NL x Mp] tangential induced velocity distribution
if Wake_flag == 1
states.URSTAR  = zeros(NL,Mp);    % [NL x Mp] radial     induced velocity distribution
end
states.VSTAR   = zeros(NL,Mp);    % [NL x Mp] total inflow       velocity distribution
states.TANBC   = zeros(NL,Mp);    % [NL x Mp] tangent of free-stream inflow angle
states.TANBIC  = zeros(NL,Mp);    % [NL x Mp] tangent of total inflow angle

states.ALPHA   = zeros(NL,Mp);    % [NL x Mp] ALPHA = alpha - alphaI
states.CL      = zeros(NL,Mp);    % [NL x Mp] lift coefficient distribution
states.CD      = zeros(NL,Mp);    % [NL x Mp] drag coefficient distribution


states.part2  = '------ Numerical metrics ------';
states.converged  = zeros(NL,1);
states.iteration  = zeros(NL,1);

states.part3   = '------ Performance metrics ------';
states.KT      = zeros(NL,1);     % [NL x 1]  thrust coefficient
states.KQ      = zeros(NL,1);     % [NL x 1]  torque coefficient
states.CT      = zeros(NL,1);     % [NL x 1]  thrust coefficient
states.CQ      = zeros(NL,1);     % [NL x 1]  torque coefficient
states.CP      = zeros(NL,1);     % [NL x 1]  power  coefficient

if abs(ginput.VMIV - 1) > 1e-8,  % i.e. if VMIV is not equal to 1 
    states.EFFYo  = zeros(NL,1); % open water efficiency
    states.Ja     = zeros(NL,1); % == VMIV*Js
end
states.EFFY    = zeros(NL,1);     % [NL x 1]  efficiency == VMIV*EFFYo

states.VMWV    = zeros(NL,1);               % [1 x 1]

% -------------------------------------------------------------------------

%%              
% ---------------------- Analyze lower LAMBDA states first, then higher values
if isfield(ginput,'L')
    LAMBDA0 =    ginput.L;
else
    LAMBDA0 = pi/ginput.Js;
end
LAMBDA_low  = LAMBDAall(find(LAMBDAall <= LAMBDA0,1,'last'):-1:1  );
LAMBDA_high = LAMBDAall(find(LAMBDAall <= LAMBDA0,1,'last'): 1:end);


% If all LAMBDA are greater than the design LAMBDA, then 
% analyze all starting from the lowest LAMBDA ...
if isempty(LAMBDA_low)
    LAMBDA_low  = [];
    LAMBDA_high = [LAMBDAall(1);LAMBDAall(:)];
end

% ------------------- Analyze lower LAMBDA states first, then higher values
    [s_low ] = Analyze_some(ginput,LAMBDA_low ,feather);
    [s_high] = Analyze_some(ginput,LAMBDA_high,feather);     

% ------------------------------- Output results in increasing LAMBDA order
states.G         = [     s_low.G(end:-1:1,:);     s_high.G(2:end,:)];   % [NL x Mp]
states.UASTAR    = [s_low.UASTAR(end:-1:1,:);s_high.UASTAR(2:end,:)];   % [NL x Mp] 
states.UTSTAR    = [s_low.UTSTAR(end:-1:1,:);s_high.UTSTAR(2:end,:)];   % [NL x Mp] 
if Wake_flag == 1
states.URSTAR    = [s_low.URSTAR(end:-1:1,:);s_high.URSTAR(2:end,:)];   % [NL x Mp] 
end
states.VSTAR     = [ s_low.VSTAR(end:-1:1,:); s_high.VSTAR(2:end,:)];   % [NL x Mp]     
states.TANBC     = [ s_low.TANBC(end:-1:1,:); s_high.TANBC(2:end,:)];   % [NL x Mp] 
states.TANBIC    = [s_low.TANBIC(end:-1:1,:);s_high.TANBIC(2:end,:)];   % [NL x Mp] 

states.ALPHA     = [ s_low.ALPHA(end:-1:1,:); s_high.ALPHA(2:end,:)];   % [NL x Mp]
states.CL        = [    s_low.CL(end:-1:1,:);    s_high.CL(2:end,:)];   % [NL x Mp] 
states.CD        = [    s_low.CD(end:-1:1,:);    s_high.CD(2:end,:)];   % [NL x Mp] 


states.converged   = [    s_low.converged(end:-1:1,:);    s_high.converged(2:end,:)];   % [NL x 1]
states.iteration   = [    s_low.iteration(end:-1:1,:);    s_high.iteration(2:end,:)];   % [NL x 1]


states.KT        = [    s_low.KT(end:-1:1,:);    s_high.KT(2:end,:)];   % [NL x 1]
states.KQ        = [    s_low.KQ(end:-1:1,:);    s_high.KQ(2:end,:)];   % [NL x 1]
states.CT        = [    s_low.CT(end:-1:1,:);    s_high.CT(2:end,:)];   % [NL x 1]
states.CQ        = [    s_low.CQ(end:-1:1,:);    s_high.CQ(2:end,:)];   % [NL x 1]
states.CP        = [    s_low.CP(end:-1:1,:);    s_high.CP(2:end,:)];   % [NL x 1]

if abs(ginput.VMIV - 1) > 1e-8,  % i.e. if VMIV is not equal to 1 
    
states.EFFYo     = [ s_low.EFFYo(end:-1:1,:); s_high.EFFYo(2:end,:)];   % [NL x 1]
states.Ja        = [    s_low.Ja(end:-1:1,:);    s_high.Ja(2:end,:)];   % [NL x 1]
end

states.EFFY      = [  s_low.EFFY(end:-1:1,:);  s_high.EFFY(2:end,:)];   % [NL x 1]
states.VMWV      = [  s_low.VMWV(end:-1:1,:);  s_high.VMWV(2:end,:)];   % [NL x 1]



end % of the function [states] = Analyze(design,LAMBDAall)
% =========================================================================
% =========================================================================

% =========================================================================
% ============================================== FUNCTION Analyze_some(...)
%%
function [some] = Analyze_some(ginput,LAMBDA_some,feather)

%%
% -------------------------------------------- Unpack ginput data structure   

% ------------------------------------------------------------------- Flags
Propeller_flag = ginput.Propeller_flag; % 0 == turbine, 1 == propeller
Viscous_flag   = ginput.Viscous_flag;   % 0 == viscous forces off (CD = 0), 
                                        % 1 == viscous forces on

% 0 == no hub, 1 == hub
if isfield(ginput,'Hub_flag'),   Hub_flag = ginput.Hub_flag;  
                          else   Hub_flag = 0;  end
                          
% 0 == no duct, 1 == duct
if isfield(ginput,'Duct_flag'),  Duct_flag = ginput.Duct_flag;  
                          else   Duct_flag = 0;  end

% 0 == do not display plots, 1 == display plots
if isfield(ginput,'Plot_flag'),  Plot_flag = ginput.Plot_flag;  
                          else   Plot_flag = 0;  end 

% 0 == Horseshoe(...,Wrench(...)), 1 == Wake_Horseshoe(...)      
if isfield(ginput,'Wake_flag'),  Wake_flag = ginput.Wake_flag;  
                          else   Wake_flag = 0;  end
                          
% 0 == use CLCD_vs_ALPHA.m , 1 == use CL and CD data given in ginput
if isfield(ginput,'CLCD_flag'),  CLCD_flag = ginput.CLCD_flag;  
                          else   CLCD_flag = 0;  end
                          
% Lifting Surface Geometry Corrections
if isfield(ginput,'LSGeoCorr'),  LSGeoCorr = ginput.LSGeoCorr;  
else
    if Propeller_flag == 1
        LSGeoCorr = 'Morgan1968';
    else
        LSGeoCorr = 'none';  
    end
end                 
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% If CL(ALPHA) and CD(ALPHA) data are given in ginput, then extract them.
% Each of thse matrices is size[M,N], where
% M == number of unique RELM  radii
% N == number of unique ALPHA angles
if CLCD_flag == 1  
    ALPHAdata = ginput.ALPHAdata;
     RELMdata = ginput.RELMdata;
       CLdata = ginput.CLdata;
       CDdata = ginput.CDdata;
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------                         
% '------ Performance inputs ------'                         
Z = ginput.Z;            % [1 x 1] scalar number of blades

if Propeller_flag == 1
    Js = ginput.Js;
    L = pi/Js;
else
    L     = ginput.L;
    Js    = pi/L;
end           

% -------------------------------------------------------------------------
                 
              
% -------------------------------------------------------------------------
% '------ Vortex lattice / blade geometry ------'
RC          = ginput.RC;          % [1 x Mp]   control point radii
RV          = ginput.RV;          % [1 x Mp+1] vortex  point radii
DR          = ginput.DR;          % [1 x Mp]   difference in vortex point radii
Rhub_oR     = ginput.Rhub_oR;     % [1 x 1]
VAC         = ginput.VAC;         % [1 x Mp] 
VTC         = ginput.VTC;         % [1 x Mp] 
CoD         = ginput.CoD;         % [1 x Mp], chord / diameter       at RC

if isfield(ginput,'PoD')
    PoD   = ginput.PoD;

    theta = atan(PoD./(pi*RC));     % [rad] design pitch angle
else
    theta = ginput.THETArad; % [rad]
end
 

if CLCD_flag == 0
    CD0         = ginput.CD;          % [Mp x 1], [ ] ginput drag coefficient       at RC
    f0oc        = ginput.f0oc;
    t0oc        = ginput.t0oc;
    EAR         = ginput.EAR;  % expanded area ratio
    BTF         = ginput.BTF;  % blade thickness fraction
    dCLdALPHA   = ginput.dCLdALPHA;
end

VMIV = ginput.VMIV;

Mp   = length(RC);
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% '------ Computational inputs ------'
if isfield(ginput,'ITER'), ITER = ginput.ITER;  else  ITER = 50;   end
if isfield(ginput,'Rhv' ), Rhv  = ginput.Rhv;   else  Rhv  = 0.5;  end

if isfield(ginput,'ALPHAstall'), ALPHAstall  = ginput.ALPHAstall; % [rad]
                           else  ALPHAstall  = 8*pi/180;  end
                           
if length(ALPHAstall) == 1, ALPHAstall = ALPHAstall * ones(size(RC)); end                           
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% '------ Duct parameters ------'
    Rduct_oR = 1;          % [1 x 1]
    UADUCT   = 0*RC;       % [1,Mp]   velocity induced by duct at propeller
    CTD      = 0;          % [1 x 1]  duct thrust coefficient on design
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% --------------------------------------- Initialize performance parameters 
TANBC  = (VAC)./(L*RC + VTC);

BetaC  = atan(TANBC); % [rad], beta at RC

BetaIC = BetaC;  % [rad] assumed initialization for: BetaIC 
TANBIC = TANBC;


% -------------------------------------------------------------------------
if CLCD_flag == 0
    % '------ 2D section performance ------'
       CLItilde = ginput.CLItilde;    % [ ]   either scalar or versus RC
    alphaItilde = ginput.alphaItilde; % [deg] either scalar or versus RC
      f0octilde = ginput.f0octilde;   % [ ]   either scalar or versus RC



    alphaIraw = (pi/180)*alphaItilde .* (f0oc./f0octilde); % [rad] ideal angle of attack,  versus RC
       CLIraw =             CLItilde .* (f0oc./f0octilde); % [ ]   ideal lift coefficient, versus RC

    % -------------------------------------------------------------------------
    % Apply lifting surface geometry correction LSGeoCorr
    if strcmp(LSGeoCorr,'none')

        % -------------------------------------------------------------------------
        % -------------- no camber correction or pitch correction
           CLI =    CLIraw;             %           CLI is for effective 2D flow  
        alphaI = alphaIraw;             % [rad], alphaI is for effective 2D flow
        % -------------------------------------------------------------------------

    elseif strcmp(LSGeoCorr,'EckhardtMorgan1955') 

        % -------------------------------------------------------------------------
        % ---- Apply flow curvature correction factor from Eckhardt and Morgan 1955
        % CL(curved flow) = K1K2 * CL(2D flow)
        [K1K2] = EckhardtMorgan1955(EAR,RC,TANBIC);

        %alphaI = alphaIraw;
        alphaI = alphaIraw ./ K1K2; % [rad], alphaI is for equivalent 2D flow
           CLI =    CLIraw ./ K1K2; %           CLI is for equivalent 2D flow  
        % -------------------------------------------------------------------------

    elseif strcmp(LSGeoCorr,'Morgan1968')     

        % -------------------------------------------------------------------------
        % ----------------- Apply flow curvature correction factor from Morgan 1968
        % Kc = 1; Ka = 1; Kt = 0;
        [Kc,Ka,Kt] = Morgan1968(RC,TANBIC,EAR,Z);  % extrapolate in Z
        % [Kc,Ka,Kt] = Morgan1968v1(RC,TANBIC,EAR,Z); % use Z==4

           CLI =    CLIraw       ./ Kc;          %           CLI is for equivalent 2D flow  
        alphaI = alphaIraw .* Ka ./ Kc + Kt*BTF; % [rad], alphaI is for actual     3D flow
        % -------------------------------------------------------------------------

    elseif strcmp(LSGeoCorr,'vanManen1958')  

        % -------------------------------------------------------------------------
        % -------------- Apply flow curvature correction factor from van Manen 1958
        % f_effective = Kvm * f_geometric
        %
        % alphaIraw = (pi/180)*alphaItilde * (f0oc/f0octilde); % [rad] ideal angle of attack
        %    CLIraw =             CLItilde * (f0oc/f0octilde); % [ ]   ideal lift coefficient

        % VMk = 1; VMca = 0; VMcf = 1;

        VMca  = 0.23333333;
        VMcf  = 1.21666667; 
        VMk   = vanManen1958(RC,TANBIC,EAR,Z);

           CLI =    CLIraw .* VMk/VMcf;                     %           CLI is for effective 2D flow  
        alphaI = alphaIraw .* VMk/VMcf +  f0oc * VMca/VMcf; % [rad], alphaI is for effective 2D flow
        % -------------------------------------------------------------------------
    else
                disp('ERROR: Unknown LSGeoCorr (lifting surface geometry correction).')
                return
    end
    % -------------------------------------------------------------------------

else
    alphaI = 0;
    CLI    = 0;
end

% -------------------------------------------------------------------------
% Set reference inflow angle and angle of atack:
BetaIC0 = theta - alphaI; % [rad] "design BetaIC" == (2D BetaIC) == (3D geometry theta) - (3D corrected alphaI)
CL0     = CLI;            % [ ]   "design CLI" where CL == CL0 when BetaIC == BetaIC0 
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Initialize the unknown parameters:
alpha   = theta - BetaC;  % [rad] assumed initialization for angle of attack

ALPHA   = alpha - alphaI; % [rad] == BetaIC0 - BetaIC

if CLCD_flag == 1
    
    CL = interp2(ALPHAdata,RELMdata,CLdata,ALPHA,RC,'linear');
else
    CL = CLCD_vs_ALPHA(ALPHA,ALPHAstall,CL0,CD0,dCLdALPHA,Propeller_flag);
end

VSTAR = sqrt(VAC.^2 + (L*RC+VTC).^2);  % assumed initialization for: V* / Vs

G     = (CL.*VSTAR.*CoD)'/(2*pi);      % assumed initialization for: circulation distribution

[UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBC,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);
        
UASTAR = (UAHIF*G)';
UTSTAR = (UTHIF*G)';

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------              
%%


%%
% ---------------------------- Initialize output data structure "some.____"
NL = length(LAMBDA_some);

some.part1  = '------ Off-design states ------';
some.L      = LAMBDA_some(:);  % [NL x 1] tip speed ratio 
some.Js     = pi./some.L;      % [NL x 1] advance coefficient 

some.RC     = RC;
some.G      = zeros(NL,Mp);    % [NL x Mp] circulation distribution
some.UASTAR = zeros(NL,Mp);    % [NL x Mp] 
some.UTSTAR = zeros(NL,Mp);    % [NL x Mp] 
if Wake_flag == 1
some.URSTAR = zeros(NL,Mp);    % [NL x Mp] 
end
some.VSTAR  = zeros(NL,Mp);    % [NL x Mp] 
some.TANBC  = zeros(NL,Mp);    % [NL x Mp] 
some.TANBIC = zeros(NL,Mp);    % [NL x Mp]

some.part2  = '------ 2D section performance ------';
some.ALPHA  = zeros(NL,Mp);    % [NL x Mp] 
some.CL     = zeros(NL,Mp);    % [NL x Mp] 
some.CD     = zeros(NL,Mp);    % [NL x Mp] 


some.part3  = '------ Numerical metrics ------';
some.converged  = zeros(NL,1);
some.iteration  = zeros(NL,1);


some.part4  = '------ Performance metrics ------';
some.CT     = zeros(NL,1);     % [NL x 1]
some.CQ     = zeros(NL,1);     % [NL x 1]
some.CP     = zeros(NL,1);     % [NL x 1]
some.KT     = zeros(NL,1);     % [NL x 1]
some.KQ     = zeros(NL,1);     % [NL x 1]
some.EFFY   = zeros(NL,1);     % [NL x 1]

some.EFFYo  = zeros(NL,1);     % open water efficiency
some.Ja     = zeros(NL,1);     % == VMIV*Js
some.VMWV   = zeros(NL,1); 

% ------------------------------------------------ Implement RepairSpline.m
% To smooth data X(RC):  X_smooth = X*Bsmooth;
Bsmooth = RepairSplineMatrix(RC);
% -------------------------------------------------------------------------            

%%
Crash_count = 0;

% -------------- For each LAMBDA, find the system state, and compute forces
for i = 1:length(LAMBDA_some)
%     % -------------------------------------------- Initialize the new state 
%     disp(['-------- Determining state:  LAMBDA = ',num2str(LAMBDA_some(i)),' , Js = ',num2str(pi/LAMBDA_some(i))]),
%     disp(' '),
%%    

    if Crash_count > 3 
        disp(['Crash limit of 3 reached...skipping state:  L = ',num2str(LAMBDA_some(i)),' , Js = ',num2str(pi/LAMBDA_some(i))])
        continue
    end
    
    converged = 1;
    
    L  =    LAMBDA_some(i);
    Js = pi/LAMBDA_some(i);

    TANBC  = (VAC)./(L*RC + VTC);


    % Allow twice as many iterations for the first state analyzed
    if i == 1
        ITER = 2*ITER;
    elseif i == 2
        ITER = ITER/2;
    end
    

    % -------------------------------------------------------------------------
    if Plot_flag == 1
        close all,
        % ------------------------------------------------- Uncomment to plot G
        Hgamma = fig;                  set(Hgamma,'Position',[130 680 560 420])
            HHgamma = plot(RC,G,'k.-');   
            set(gca,'XLim',[Rhub_oR 1])
            xlabel('RC [ ]','Fontsize',14,'Fontname','Times')
            ylabel('G [ ]' ,'Fontsize',14,'Fontname','Times')
        % ---------------------------------------------------------------------     

        % ------------------------------ Uncomment to plot UASTAR, UTSTAR, & URSTAR
        Hvel = fig;                      set(Hvel,'Position',[695 680 560 420])
            if Duct_flag == 1
                   HHvel = plot(RC,UASTAR,'b.-',RC,UTSTAR,'r.-',RC,UADUCT,'g.-');
            else
                   HHvel = plot(RC,UASTAR,'b.-',RC,UTSTAR,'r.-');
            end
            axis([Rhub_oR 1 -1 1])
            xlabel('RC [ ]'    ,'Fontsize',14,'Fontname','Times')
            ylabel('UASTAR (blue), UTSTAR (red)','Fontsize',14,'Fontname','Times')
        % --------------------------------------------------------------------- 

        % ---------------------------------------------------------------------
        Hbeta = fig;                    set(Hbeta,'Position',[695 160 560 420])
    %     % ------------------------------------ TANBC & TANBIC
    %         HHbeta = plot(RC,TANBC,'k.-',RC,TANBIC,'r.-'); 
    % 
    %         set(gca,'XLim',[Rhub_oR 1])
    %         xlabel('RC [ ]'    ,'Fontsize',14,'Fontname','Times')
    %         ylabel('TANBIC (red), TANBC (black)','Fontsize',14,'Fontname','Times')
    %     % --------------------------------------------------------------------- 
                % ------------------------------------------------------- ALPHA
                HHbeta = plot(RC,ALPHA*180/pi,'r.'); hold on,            
                plot([0 1],[0 0],'k:')
                plot([0 1],[ALPHAstall(1) ALPHAstall(1)]*180/pi,'k:',[0 1],-[ALPHAstall(1) ALPHAstall(1)]*180/pi,'k:')
                axis([0 1 -20 20])
                xlabel('RC [ ]'    ,'Fontsize',14,'Fontname','Times')
                ylabel('ALPHA, [deg]','Fontsize',14,'Fontname','Times')
                % ----------------------------------------------- percent error
    %             HHbeta = plot([],[],'.');
    %             plot( [Mp,Mp],[-0.1,0.1],'k:',2*[Mp,Mp],[-0.1,0.1],'k:',...
    %                 3*[Mp,Mp],[-0.1,0.1],'k:',4*[Mp,Mp],[-0.1,0.1],'k:',...
    %                 5*[Mp,Mp],[-0.1,0.1],'k:',6*[Mp,Mp],[-0.1,0.1],'k:'),
    %             hold on,
    %             text(10     ,-0.09,'VSTAR'),
    %             text(10+1*Mp,-0.09,'ALPHA'),
    %             text(10+2*Mp,-0.09,'CL'),
    %             text(10+3*Mp,-0.09,'G'),
    %             text(10+4*Mp,-0.09,'UASTAR'),
    %             text(10+5*Mp,-0.09,'UTSTAR'),
    %             plot([0, 6*Mp],[0,0],'k:')
    %             axis([0, 6*Mp, -0.1, 0.1]),
    %             ylabel('percent error, [ ]'    ,'Fontsize',14,'Fontname','Times')            
            % ---------------------------------------------------------------------       


        % ---------------------------------------------------------------------
        Hcod = fig;                      set(Hcod,'Position',[130 160 560 420])
                % ---------------------------------------------- residual error
                HHcod = plot([],[],'.');            
                plot( [Mp,Mp],[-0.1,0.1],'k:',2*[Mp,Mp],[-0.1,0.1],'k:',...
                    3*[Mp,Mp],[-0.1,0.1],'k:',4*[Mp,Mp],[-0.1,0.1],'k:',...
                    5*[Mp,Mp],[-0.1,0.1],'k:',6*[Mp,Mp],[-0.1,0.1],'k:'),
                hold on,
                text(10     ,-0.09,'VSTAR'),
                text(10+1*Mp,-0.09,'ALPHA'),
                text(10+2*Mp,-0.09,'CL'),
                text(10+3*Mp,-0.09,'G'),
                text(10+4*Mp,-0.09,'UASTAR'),
                text(10+5*Mp,-0.09,'UTSTAR'),
                plot([0, 6*Mp],[0,0],'k:')
                axis([0, 6*Mp, -0.1, 0.1]),
                ylabel('residual error, [ ]'    ,'Fontsize',14,'Fontname','Times')
            % ---------------------------------------------------------------------

        pause(0.001), drawnow,
    end  % if Plot_flag == 1
    % -------------------------------------------------------------------------     

    
    % =====================================================================
    % ==   FIND STATE OF SYSTEM USING NEWTON SOLVER:
    % ==                   iterate to solve residual equations for unknowns
    N_iter   = 1;                            
    ERROR    = 1;
    ERRORtol = 0.005;
    relax    = 0.9;
%%    
    while N_iter <= ITER && any(abs(ERROR) > ERRORtol)          % (WHILE LOOP N1)
        % disp(['------- Newton iteration: ',num2str(N_iter)])  % status message
        % disp(' '),

        % ---------------------------------- Store last state of the system
         VSTARlast =  VSTAR;
         alphalast =  ALPHA;
            CLlast =     CL;
             Glast =      G;
        UASTARlast = UASTAR;
        UTSTARlast = UTSTAR;

        % ------------------------------------- Evaluate velocities for R6 & R7
        UASTARtemp = (UAHIF*G)';  
        UTSTARtemp = (UTHIF*G)'; 

        if CLCD_flag == 1

            CLtemp        =               interp2(ALPHAdata,RELMdata,CLdata,ALPHA,RC,'linear');  % indexed versus RC
            [dCLda,dCDda] = Find_dCLCDdALPHA_data(ALPHAdata,RELMdata,CLdata,ALPHA,RC);           % indexed versus RC
        end        
        
        
        % --------------- Solve Newton problem to update each blade section
        for m = 1:Mp  
            
            % -------------- Initialize linear system of equations matrices
            R  = zeros(6,1);            % R  = vector of residuals
            J  =   eye(6,6);            % J  = matrix of derivatives
            DX = zeros(6,1);            % DX = vector of change in unknowns 
          % X  = zeros(6,1);            % X  = vector of unknowns 

            % X(1) =  VSTAR(m);   % X1
            % X(2) =  ALPHA(m);   % X2
            % X(3) =     CL(m);   % X3
            % X(4) =      G(m);   % X4
            % X(5) = UASTAR(m);   % X5
            % X(6) = UTSTAR(m);   % X6
            % -------------------------------------------------------------------------

            % ------------------------------------------------------ Evaluate residuals
            R(1) =  VSTAR(m) - sqrt((VAC(m)+UASTAR(m)).^2 + (L*RC(m)+VTC(m)+UTSTAR(m)).^2);     % R1

            R(2) =  ALPHA(m) - ( BetaIC0(m) - atan(TANBIC(m)) + feather );                      % R2
            
            
            if CLCD_flag == 1

                R(3) =     CL(m) - CLtemp(m);                                                   % R3             
            else                
                R(3) =     CL(m) - CLCD_vs_ALPHA(ALPHA(m),ALPHAstall(m),CL0(m),CD0(m),dCLdALPHA(m),Propeller_flag);      % R3
            end

            
            R(4) =      G(m) - (1/(2*pi))*VSTAR(m)*CL(m)*CoD(m);                                % R4

            R(5) = UASTAR(m) - UASTARtemp(m);                                                   % R5

            R(6) = UTSTAR(m) - UTSTARtemp(m);                                                   % R6            
            % -------------------------------------------------------------------------
                        
            
            % --------------------------------------- Evaluate residual derivatives
            % disp('Evaluating residual derivatives matrix, J...')
            % disp(' '),
            % J = eye(6,6) ==>  J(i,i) == 1 for all i = 1:6

            J(1,5) = -        (VAC(m)+UASTAR(m))/sqrt((VAC(m)+UASTAR(m)).^2 + (L*RC(m)+VTC(m)+UTSTAR(m)).^2);
            J(1,6) = -(L*RC(m)+VTC(m)+UTSTAR(m))/sqrt((VAC(m)+UASTAR(m)).^2 + (L*RC(m)+VTC(m)+UTSTAR(m)).^2);

            J(2,5) =   (1/(1+(TANBIC(m))^2))*(        1/(L*RC(m)+VTC(m)+UTSTAR(m)));
            J(2,6) = - (1/(1+(TANBIC(m))^2))*(TANBIC(m)/(L*RC(m)+VTC(m)+UTSTAR(m)));

            if CLCD_flag == 1
                J(3,2) = - dCLda(m);
            else
                J(3,2) = - Find_dCLCDdALPHA(ALPHA(m),ALPHAstall(m),CL0(m),CD0(m),dCLdALPHA(m),Propeller_flag);
            end
            
            J(4,1) = -(1/(2*pi))*   CL(m)*CoD(m);
            J(4,3) = -(1/(2*pi))*VSTAR(m)*CoD(m);

            J(5,4) = - UAHIF(m,m);
            J(6,4) = - UTHIF(m,m);
            % ---------------------------------------------------------------------

            
            % ----------------------------- Update Newton solver vector of unknowns
            DX = linsolve(J,-R);
            
            if any(isnan(DX))  || ~isreal(DX) || any(abs(DX) > 10^4)
                 DX = zeros(6,1);

                 disp('<WARNING>')
                 disp('<WARNING> DX == NaN or imaginary... crash avoided...')
                 disp('<WARNING>')
                 N_iter      = 999;
                 converged   = 0;
                 Crash_count = Crash_count + 1;
            end         
            

             VSTAR(m) =  VSTAR(m) + relax * DX(1);   % X1
             ALPHA(m) =  ALPHA(m) + relax * DX(2);   % X2
                CL(m) =     CL(m) + relax * DX(3);   % X3
                 G(m) =      G(m) + relax * DX(4);   % X4
            UASTAR(m) = UASTAR(m) + relax * DX(5);   % X5
            UTSTAR(m) = UTSTAR(m) + relax * DX(6);   % X6

        end   % ------------------------------------ END for each blade section

        % ----------------- Update unkowns not implemented in the Newton solver
        TANBIC = (VAC + UASTAR)./(L*RC + VTC + UTSTAR);

        % Smooth the inflow angle for numerical stability:
        TANBICsmooth = TANBIC * Bsmooth;
        % -----------------------------------------------------------------         
        
        
        if Wake_flag == 0 

            [UAHIF,UTHIF] = Horseshoe110628(Mp,Z,TANBICsmooth,RC,RV,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR);

        else
            % % disp('Beginning to update {UAHIF,UTHIF,URHIF}'), 
            % 
            % [CT,CQ,CP,KT,KQ,EFFY,TAU,CTP,CTH,CTD,KTP] = ...
            % Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,CD,CoD,G,Z,Js,VMIV,Hub_flag,Rhv,CTD);
            % 
            % URSTAR = (URHIF*G)';        
            % 
            % % [WX,WY,WZ]          =  Wake_Geometry(Mp,RC,RV,TANBIV,VAC,UASTAR,URSTAR,CTP);
            % [WX,WY,WZ] = Wake_Geometry2U(VAC,VTC,UASTAR,UTSTAR,RC,RV,Js);
            % 
            % %[UAHIF,UTHIF,URHIF] = Wake_Horseshoe(Mp,Z,RC,SCF,WX,WY,WZ,epsilon);
            % [UAHIF,UTHIF,URHIF] = Wake_Horseshoe(Mp,Z,TANBIV,RC,RV,SCF,Hub_flag,Rhub_oR,Duct_flag,Rduct_oR,WX,WY,WZ,epsilon);
            % %disp('Done updating {UAHIF,UTHIF,URHIF}'), 
        end 
        
        % ----------------------------------------------------------------- 
        if CLCD_flag == 0
            % Apply lifting surface geometry correction
            if strcmp(LSGeoCorr,'EckhardtMorgan1955') 

                % -----------------------------------------------------------------
                % ---- Apply flow curvature correction factor from Eckhardt and Morgan 1955
               [K1K2] = EckhardtMorgan1955(EAR,RC,TANBIC);                      % unconstrained wake geometry
               % [K1K2] = EckhardtMorgan1955(EAR,RC,mean(TANBIC./TANBC).*TANBC);  % Lerbs-constrained wake geometry

                alphaI = alphaIraw ./ K1K2; % [rad], alphaI is for equivalent 2D flow
                   CLI =    CLIraw ./ K1K2; %           CLI is for equivalent 2D flow  

                BetaIC0 = theta - alphaI; % [rad] "design BetaIC"
                CL0     = CLI;            % [ ]   "design CLI" where CL == CL0 when BetaIC == BetaIC0 
                % -----------------------------------------------------------------

            elseif strcmp(LSGeoCorr,'Morgan1968')     

                % -------------------------------------------------------------------------
                % ----------------- Apply flow curvature correction factor from Morgan 1968
                [Kc,Ka,Kt] = Morgan1968(RC,TANBIC,EAR,Z);  % extrapolate in Z
                %[Kc,Ka,Kt] = Morgan1968(RC,mean(TANBIC./TANBC).*TANBC,EAR,Z);  % extrapolate in Z, Lerbs-aligned

                   CLI =    CLIraw       ./ Kc;          %           CLI is for equivalent 2D flow  
                alphaI = alphaIraw .* Ka ./ Kc + Kt*BTF; % [rad], alphaI is for equivalent 2D flow

                BetaIC0 = theta - alphaI; % [rad] "design BetaIC"
                CL0     = CLI;            % [ ]   "design CLI" where CL == CL0 when BetaIC == BetaIC0 
                % -------------------------------------------------------------------------

            elseif strcmp(LSGeoCorr,'vanManen1958') 

                % -------------------------------------------------------------------------
                % -------------- Apply flow curvature correction factor from van Manen 1958
                VMk   = vanManen1958(RC,TANBIC,EAR,Z);

                   CLI =    CLIraw .* VMk/VMcf;                     %           CLI is for effective 2D flow  
                alphaI = alphaIraw .* VMk/VMcf +  f0oc * VMca/VMcf; % [rad], alphaI is for effective 2D flow

                % [Kc,Ka,Kt] = Morgan1968(RC,TANBIC,EAR,Z);  % extrapolate in Z
                % alphaI = alphaIraw .* VMk/VMcf +  f0oc * VMca/VMcf + Kt*BTF; % [rad], alphaI is for effective 2D flow

                BetaIC0 = theta - alphaI; % [rad] "design BetaIC"
                CL0     = CLI;            % [ ]   "design CLI" where CL == CL0 when BetaIC == BetaIC0 
                % -------------------------------------------------------------------------
            end
        end
        % ----------------------------------------------------------------- 
        
        
        % ------------------------------------------------ End update variables
        % ----------------------------------------------------------------- 

        % ------------------------------------------ Evaluate normalized residuals
        UASTARtemp = (UAHIF*G)';  
        UTSTARtemp = (UTHIF*G)'; 
        
               
        if CLCD_flag == 1
            CLtemp = interp2(ALPHAdata,RELMdata,CLdata,ALPHA,RC,'linear');
        else
            CLtemp = CLCD_vs_ALPHA(ALPHA,ALPHAstall,CL0,CD0,dCLdALPHA,Propeller_flag);
        end
                

        ERROR    = [ VSTAR -  sqrt((VAC+UASTAR).^2 + (L*RC+VTC+UTSTAR).^2), ...
                     ALPHA -  (BetaIC0 - atan(TANBIC) + feather), ...
                        CL -  CLtemp, ... 
                        G' -  (1/(2*pi))*VSTAR.*CL.*CoD, ...
                    UASTAR - UASTARtemp, ...
                    UTSTAR - UTSTARtemp] ./ ...
                   [max(abs(VSTAR ),1e-4), ...
                    max(abs(ALPHA ),1e-2), ...
                    max(abs(CL    ),1e-4), ...
                    max(abs(G'    ),1e-6), ...
                    max(abs(UASTAR),1e-4), ...
                    max(abs(UTSTAR),1e-4)]; 
        % ------------------------------------------ END evaluate normalized residuals
        

        % ------------------------------------------ Evaluate percent error
        PERROR    = [ VSTARlast -  VSTAR, ...
                     alphalast -  ALPHA, ...
                        CLlast -     CL, ...
                        Glast' -     G', ...
                    UASTARlast - UASTAR, ...
                    UTSTARlast - UTSTAR] ./ ...
                   [max(abs(VSTAR ),1e-4), ...
                    max(abs(ALPHA ),1e-2), ...
                    max(abs(CL    ),1e-4), ...
                    max(abs(G'    ),1e-6), ...
                    max(abs(UASTAR),1e-4), ...
                    max(abs(UTSTAR),1e-4)];
        % ---------------------------------------------- END evaluate percent error

        % --------------------------------------------------- Plot system state
        if Plot_flag == 1
            delete(HHgamma)
            delete(HHvel  )
            delete(HHbeta )
            delete(HHcod  )
            pause(0.0001)
            
            % ---------------------------------------------------- Figure 1
            figure(Hgamma),
                HHgamma = plot(RC,G,'b');
            % ------------------------------------------------------------- 

            % ---------------------------------------------------- Figure 2
            figure(Hvel),               
                  HHvel = plot(RC,UASTAR,'b.-',RC,UTSTAR,'r.-');              
            % --------------------------------------------------------------------- 

            % -------------------------------------------------------- Figure 3
            figure(Hbeta);
                % ------------------------------------------------------ TANBIC
                %  HHbeta = plot(RC,TANBIC0,'r.'); hold on,
                % ------------------------------------------------------- ALPHA
                HHbeta = plot(RC,ALPHA*180/pi,'r'); hold on,        
                % ----------------------------------------------- percent error
                %  HHbeta = plot(1:length(PERROR),PERROR,'.b');
            % -----------------------------------------------------------------


            % ---------------------------------------------------- Figure 4
            figure(Hcod),
                % ------------------------------------------ residual error
                HHcod = plot(1:length(ERROR),ERROR,'b'); 
                
            % ---------------------------------------------------------------------

            pause(0.00001), drawnow, pause(0.0001),
%                  pause,
        end
        % --------------------------------------------- END if Plot_flag == 1


        % -------------------------------------- Prepare for the next iteration
        N_iter  = N_iter + 1;              % iteration in the N loop

%         if N_iter-1 < 10
%             disp(['The max error for iteration  ',num2str(N_iter-1),' is: ',num2str(max(abs(ERROR)))]),  
%         else
%             disp(['The max error for iteration ' ,num2str(N_iter-1),' is: ',num2str(max(abs(ERROR)))]),  
%         end
%         % disp(' '),

    end                                                   % (END WHILE LOOP N1)
    
    if N_iter > ITER  &&  converged == 1
        disp('WARNING: Newton solver did NOT converge.'),
        disp(' ')
        
        converged   = 0;
        Crash_count = Crash_count + 1;
    end
    
    N_iter  = N_iter - 1;        
    % =================================================== END NEWTON SOLVER  
    % =====================================================================
%%
    
    
    % ------------------------------------------------------ Compute forces    
    if CLCD_flag == 1

        CL = interp2(ALPHAdata,RELMdata,CLdata,ALPHA,RC,'linear');
        CD = interp2(ALPHAdata,RELMdata,CDdata,ALPHA,RC,'linear');
    else
        [CL,CD] = CLCD_vs_ALPHA(ALPHA,ALPHAstall,CL0,CD0,dCLdALPHA,Propeller_flag);
    end
      
   
    if Viscous_flag == 0
        CD = 0*CD;
    end
    
    [CT,CQ,CP,KT,KQ, CTH,TAU, Ja,Jw,VMWV, EFFYo, EFFY,ADEFFY,QF] = ...
          Forces(RC,DR,VAC,VTC,UASTAR,UTSTAR,UADUCT,CD,CoD,G,Z,Js,VMIV,Hub_flag,Rhub_oR,Rhv,CTD);

      
%     % ---------------------------------------------------------------------     
%     if     (Propeller_flag == 1) & (KT < 0 | KQ < 0)
%         Crash_count = Crash_count + 1;
%         
%     elseif (Propeller_flag == 0) & (CP > 0)
%         Crash_count = Crash_count + 1;
%     end

    % ------------------------- Record output data in structure "some.____"
    some.converged(i) = converged;
    some.iteration(i) = N_iter;
           
         some.G(i,:) = G';       % [NL x Mp] circulation distribution
    some.UASTAR(i,:) = UASTAR;   % [NL x Mp] 
    some.UTSTAR(i,:) = UTSTAR;   % [NL x Mp]
    if Wake_flag == 1
    some.URSTAR(i,:) = URSTAR;   % [NL x Mp]
    end
     some.VSTAR(i,:) = VSTAR;    % [NL x Mp]     
     some.TANBC(i,:) = TANBC;    % [NL x Mp] 
    some.TANBIC(i,:) = TANBIC;   % [NL x Mp] 
    
     some.ALPHA(i,:) = ALPHA;    % [NL x Mp] 
        some.CL(i,:) = CL;       % [NL x Mp] 
        some.CD(i,:) = CD;       % [NL x Mp] 
        


        some.KT(i) = KT;       % [NL x 1]
        some.KQ(i) = KQ;       % [NL x 1]     
        some.CT(i) = CT;       % [NL x 1]
        some.CQ(i) = CQ;       % [NL x 1]
        some.CP(i) = CP;       % [NL x 1]
      some.EFFY(i) = EFFY;     % [NL x 1]
      
    some.EFFYo(i)  = EFFYo; % open water efficiency
       some.Ja(i)  = Ja; % == VMIV*Js

     some.VMWV(i)  = VMWV;               % [1 x 1]


    
% -------------------------------------------------------------------------
    disp(' '),
    disp(['Forces for state:  L = ',num2str(LAMBDA_some(i)),' , Js = ',num2str(pi/LAMBDA_some(i))]),
    
if Propeller_flag == 1
   
    disp(['    Js = ',num2str(Js)]),    
    disp(['    KT = ',num2str(KT)]),   
    disp(['    KQ = ',num2str(KQ)]),   
    disp(['    CT = ',num2str(CT)]),
    if abs(VMIV - 1) > 1e-8,  % i.e. if VMIV is not equal to 1 
    disp(['    Ja = ',num2str(Ja)]),
    end
    disp(['  EFFY = ',num2str(EFFY)]),   
    
else  
    disp(['CT = ' ,num2str(CT)]), 
    disp(['CP = ' ,num2str(CP)]),       
end
    disp(' ')     
% -------------------------------------------------------------------------       
   
    
    if N_iter < 10
        disp(['The max error for iteration  ',num2str(N_iter),' is: ',num2str(max(abs(ERROR)))]),  
    else
        disp(['The max error for iteration ' ,num2str(N_iter),' is: ',num2str(max(abs(ERROR)))]),  
    end
    disp(' '),    
    disp('------------------------------------------------'),
    disp(' ')
   
    
%     if converged == 0
%         disp('WARNING: Newton solver did NOT converge.'),
%         disp(' ')
%         disp('------------------------------------------------'),
%         disp(' ')
%     else    
%         disp('------------------------------------------------'),
%         disp(' ')
%     end    
    
%     disp(' ')
%     disp('Analyze(...) paused.  Press any key to continue... ')
%     pause,
    disp(' ')
    disp(' ')
    
end
    
end % function [...] = Analyze_some(...)                                
