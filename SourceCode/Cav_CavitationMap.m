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


% ----------------------------------------------------- Cav_CavitationMap.m
% Created: Brenden Epps, 9/14/2009
%
% This function creates a cavitation map of the blade, showing regions of
% cavitation in red and regions of no cavitation in green.
%
% -------------------------------------------------------------------------

function [MCPU,MCPL,sigma,BladeR,BladeC,maxcp] = Cav_CavitationMap(pt,stateindex,M,N)

% M radial    sections
% N chordwise stations

%%
if isfield(pt,'i') && ~isfield(pt,'input' ), pt.input  = pt.i; end
if isfield(pt,'d') && ~isfield(pt,'design'), pt.design = pt.d; end
if isfield(pt,'s') && ~isfield(pt,'states'), pt.states = pt.s; end

NS = length(pt.states.L); % number of states

if nargin == 1
    stateindex = 0;  % use the design state
    
else  % stateindex is given  
    if stateindex < 1 || stateindex > NS
        disp('Cav_CavitationMap  ERROR: state index is out of bounds.')
        return
    end        
end

if nargin < 3
    M = length(pt.design.RC); % == Mp
    N = 80;
    
elseif nargin < 4
    N = 80;
end
%%
input = pt.input;

% --------------------------------------------------- Unpack data structure
if isfield(input,  'Vs'), Vs   = input.Vs;    else  Vs   = 1;   end % m/s
if isfield(input,   'R'), R    = input.R;     else  R    = 1;   end % m
if isfield(input,'Rhub'), Rhub = input.Rhub;  else  Rhub = 0.2; end % m


% -------------------------------------------------------------------------
% 0 == do not optimize chord lengths, 1 == optimize chord lengths
if isfield(input,'Chord_flag'), Chord_flag = input.Chord_flag;  
                          else  Chord_flag = 0;  end
                          
% 0 == no duct, 1 == duct
if isfield(pt.input,'Duct_flag'),  Duct_flag = pt.input.Duct_flag;  
                            else   Duct_flag = 0;  end
                            

if Duct_flag == 1 
    Rduct_oR = pt.design.Rduct_oR;
else
    Rduct_oR = 1;
end
    

% ----------- Blade geometry inputs: XR, XCoD, t0oc0, Xt0oD, XVA, XVT, dXVA
% If propeller geometry or inflow is not given, assume propeller 4118 form 
% with uniform axial inflow, no swirl inflow, and section CD == 0.008.
if isfield(input,'XR'), 
    XR  = input.XR;
    
    if isfield(input,'XCoD'), 
            XCoD  = input.XCoD;  
                      
        if     isfield(input,'Xt0oD'),  Xt0oD = input.Xt0oD;          
        elseif isfield(input,'t0oc0'),  Xt0oD = input.t0oc0 .* XCoD; 
        else                            Xt0oD = 0.1         .* XCoD;      
        end
                                  
    else
        XXR    = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];    % r/R
        XXCoD  = [0.1600 0.1818 0.2024 0.2196 0.2305 0.2311 0.2173 0.1806 0.1387 0.0010];
        Xt0oc0 = [0.2056 0.1551 0.1181 0.0902 0.0694 0.0541 0.0419 0.0332 0.0324 0.0000];
              
        XCoD   = interp1(XXR,XXCoD ,XR,'pchip','extrap');
        t0oc0  = interp1(XXR,Xt0oc0,XR,'pchip','extrap');
             
        if  isfield(input,'Xt0oD'),  
            Xt0oD = input.Xt0oD;
        else
            Xt0oD = t0oc0 .* XCoD;   
        end
    end 
   
else 
    XR    = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95 1.0];    % r/R
    XCoD  = [0.1600 0.1818 0.2024 0.2196 0.2305 0.2311 0.2173 0.1806 0.1387 0.0010];          % c/D
    t0oc0 = [0.2056 0.1551 0.1181 0.0902 0.0694 0.0541 0.0419 0.0332 0.0324 0.0000];          % t0/c  
    Xt0oD = t0oc0 .* XCoD;              
end
% -------------------------------------------------------------------------


Rhub_oR = pt.design.Rhub_oR;
RC      = pt.design.RC;
CL      = pt.design.CL;     % ideal lift coefficient for each blade section
CoD     = pt.design.CoD;
t0oc    = pt.design.t0oc;

Mp      = length(RC);

if   stateindex == 0  % use the design state
    VSTAR   = pt.design.VSTAR;
    ALPHA   = 0*RC;
else
    if stateindex == NS
        VSTAR   = pt.states.VSTAR(stateindex,:);
        ALPHA   = pt.states.ALPHA(stateindex,:);  % ALPHA == alpha - alphaI
    else
        index1  = floor(stateindex);
        index2  = floor(stateindex) + 1;
        frac    = stateindex - index1;    % fraction of the way from index1 to index2 
        
        % Linear interpolation between index1 and index2
        Jstate  =    pt.states.Js(index1,:)*(1-frac) +    pt.states.Js(index2,:)*(frac);
        VSTAR   = pt.states.VSTAR(index1,:)*(1-frac) + pt.states.VSTAR(index2,:)*(frac);
        ALPHA   = pt.states.ALPHA(index1,:)*(1-frac) + pt.states.ALPHA(index2,:)*(frac);  % ALPHA == alpha - alphaI
    end
end  
    
% -------------------------------------------------------------------------
BR    = linspace(Rhub_oR,1,20);
BCoD  = pchip(RC,CoD,BR);


% ------------------------------------------------------- Cavitation inputs
if isfield(input,'rho'),  rho = input.rho;  else rho  = 1000;  end % kg/m^3
if isfield(input,'dVs'),  dVs = input.dVs;  else dVs  = 0.30;  end % m/s
if isfield(input,'H'  ),  H   = input.H;    else H    = 3.048; end % m
if isfield(input,'g'  ),  g   = input.g;    else g    = 9.81;  end % m/s^2
if isfield(input,'Patm'), Patm= input.Patm; else Patm = 101325;end % Pa
if isfield(input,'Pv'),   Pv  = input.Pv;   else Pv   = 2500;  end % Pa

sigma  = (Patm + rho*g*(H-RC) - Pv)./(0.5*rho*(VSTAR*Vs).^2); %  cavitation number (where "H-RC" accounts for a blade at the 12 o'clock position, which is the worst case scenario)
% -------------------------------------------------------------------------

% ----------------------------------- Determine blade pressure distribution
Np     = 80;     % number of points along x/c coordinate
BladeC = zeros(Mp,Np+1);  % (CoD)
BladeR = zeros(Mp,Np+1);  % (RC)

MCPL   = zeros(Mp,Np+1);  % -CP on the lower (pressure) side
MCPU   = zeros(Mp,Np+1);  % -CP on the upper (suction)  side

for i = 1:Mp  % for each blade section 
    
      [xp, MCPUtemp, MCPLtemp] = VLM2D(Np,CL(i),ALPHA(i),t0oc(i));
     
      MCPU(i,:) = MCPUtemp;
      MCPL(i,:) = MCPLtemp;

    BladeC(i,:) = (0.5-xp)*CoD(i);  % c/R = 2 * c/D
    BladeR(i,:) = RC(i);            % r/R
end

maxcp = zeros(Mp,1);

for i = 1:Mp
    maxcp(i) = max(MCPU(i,:));
end
    
   %sigma ./ (CL'/1.8 + 2*t0oc')
%    maxcp ./ sigma
   
%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Plot cavitation map using Plot_Blade_Contours function:
%
% Plot_Blade_Contours(BladeRV,BladeCV,data,plottitle)
%
%       BladeRV == [Mp+1,Np+1], r/R spanwise coordinate
%       BladeCV == [Mp+1,Np+1], x/R chordwise coordinate
%       data    == [Mp+1,Np+1], data at site (BladeRV,BladeRC)
%
% -------------------------------------------------------------------------

% -------------------- Interpolate blade pressure distribution onto (RV,CV)
BladeCV = zeros(Mp+1,Np+1);  % (CoD)
BladeRV = zeros(Mp+1,Np+1);  % (RV)
MCPUV   = zeros(Mp+1,Np+1);  % (MCPU at RV)
MCPLV   = zeros(Mp+1,Np+1);  % (MCPL at RV)

% Interpolate input geometry at sections with cosine spacing along the span
RV = Rhub_oR + (1-Rhub_oR)*(sin((0:Mp)*pi/(2*Mp)));  % [Rhub_oR : 1]

% ----------------------------------------- Interpolate chord and thickness
if Chord_flag == 0  % if chord optimization WAS NOT performed, 
                    % then interpolate chord and thickness the same as in EppsOptimizer.m


    if (abs(XR(end)-1) < 1e-4) && ( XCoD(end) <= 0.01 )  % if XR == 1 and XCoD == 0

        CoDV  = InterpolateChord(XR, XCoD,RV);   % section chord / propeller diameter at ctrl pts
    else
        CoDV  =            pchip(XR, XCoD,RV);   % section chord / propeller diameter at ctrl pts
    end
    
       t0oDV  =            pchip(XR,Xt0oD,RV);   % section thickness / propeller dia. at ctrl pts

else % chord optimization WAS performed, so we need to be careful to iterpolate in such a way as to  
     % preserve a   zero-chord-length tip (as in the case of a free-tip propeller) or 
     % preserve a finite-chord-length tip (as in the case of a ducted   propeller).
     
    if (Duct_flag == 0) || (Rduct_oR > 1.001)   

        CoDV  = InterpolateChord(RC,pt.design.CoD  ,RV);  % yields zero   chord at the tip       
    else                                                
        CoDV  =            pchip(RC,pt.design.CoD  ,RV);  % yields finite chord at the tip   
    end

        t0oDV =            pchip(RC,pt.design.t0oD ,RV);       
end
        t0ocV = t0oDV./CoDV;
% -------------------------------------------------------------------------



CLV    = pchip(RC,CL   ,RV);
ALPHAV = pchip(RC,ALPHA,RV);
t0ocV  = pchip(RC,t0oc ,RV);
sigmaV = pchip(RC,sigma,RV); 

for i = 1:Mp+1  % for each blade section 
    
      [xp, MCPUtemp, MCPLtemp] = VLM2D(Np,CLV(i),ALPHAV(i),t0ocV(i));
     
      MCPUV(i,:) = MCPUtemp;
      MCPLV(i,:) = MCPLtemp;

    BladeCV(i,:) = (0.5-xp)*2*CoDV(i);  % c/R = 2 * c/D
    BladeRV(i,:) = RV(i);               % r/R
end


% --------------------------- MCPUVoS == -CP/sigma
MCPUVoS = zeros(Mp+1,Np+1);
MCPLVoS = zeros(Mp+1,Np+1);

for i = 1:Mp+1
    MCPUVoS(i,:) = MCPUV(i,:) / sigmaV(i);
    MCPLVoS(i,:) = MCPLV(i,:) / sigmaV(i);
end


Plot_Blade_Contours2_2D(BladeRV,BladeCV,MCPUVoS,'MCPU(i,:)/sigma(i)')

xlabel('r / R','FontName','Times','FontSize',18)
ylabel('c / R','FontName','Times','FontSize',18)
title('Back  (suction side)','FontName','Times','FontSize',18)
box on

ax = findobj(get(gcf,'Children'),'Tag','Colorbar');
set(get(ax,'ylabel'),'string','-CP / sigma','FontSize',18,'fontweight','normal')




Plot_Blade_Contours2_2D(BladeRV,BladeCV,MCPLVoS,'MCPL(i,:)/sigma(i)')


xlabel('r / R','FontName','Times','FontSize',18)
ylabel('c / R','FontName','Times','FontSize',18)
title('Face  (pressure side)','FontName','Times','FontSize',18)
box on

ax = findobj(get(gcf,'Children'),'Tag','Colorbar');
set(get(ax,'ylabel'),'string','-CP / sigma','FontSize',16,'fontweight','normal')



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



%%
% % -------------------------------------------------------------------------
% % -------------------------------------------------------------------------
% %%
% % ----------------------------------------------------- Plot cavitation map
% Hfig = figure('Position',[450 160 500 700]);
%     % ------------------------------------- Cavitation map of upper surface
%     subplot(2,1,1),
%     
%     hold on, axis equal, box on,
%     xlabel('radius / R'                 ,'Fontsize',16,'Fontname','Times')
%     ylabel('chord / D'                  ,'Fontsize',16,'Fontname','Times')
%     title('Cavitation map: suction side','Fontsize',16,'Fontname','Times')
%     set(gca,'Fontsize',14,'Fontname','Times')
%     axis([0.15 1.05 -0.25 0.25])
%     
%     for i = 1:Mp
%         for j = 1:Np
%             if     MCPU(i,j) < 0.9*sigma(i)
%                 plot(BladeR(i,j),BladeC(i,j),'.','color',[0 0.9 0],'MarkerSize',20),        % Green
%                 
%             elseif MCPU(i,j) < 1.0*sigma(i)
%                 plot(BladeR(i,j),BladeC(i,j),'.','color',[0.95 0.95 0],'MarkerSize',20),    % Yellow
%                 
%             elseif MCPU(i,j) < 1.1*sigma(i)
%                 plot(BladeR(i,j),BladeC(i,j),'.','color',[1 0.5 0],'MarkerSize',20),        % Orange
%                 
%             else
%                 plot(BladeR(i,j),BladeC(i,j),'.','color',[1 0 0],'MarkerSize',20),          % Red
%             end
%         end
%     end
%     
%     HHcod = plot([BR,BR(end:-1:1),BR(1)],[-BCoD/2,BCoD(end:-1:1)/2,-BCoD(1)/2],'b-','LineWidth',2);
% 
%     % ------------------------------------- Cavitation map of lower surface
%     subplot(2,1,2),
%     
%     hold on, axis equal, box on,
%     xlabel('radius / R'                 ,'Fontsize',16,'Fontname','Times')
%     ylabel('chord / D'                  ,'Fontsize',16,'Fontname','Times')
%     title('Cavitation map: pressure side','Fontsize',16,'Fontname','Times')
%     set(gca,'Fontsize',14,'Fontname','Times')
%     axis([0.15 1.05 -0.25 0.25])
%     
%     for i = 1:Mp
%         for j = 1:Np
%             if     MCPL(i,j) < 0.9*sigma(i)
%                 plot(BladeR(i,j),BladeC(i,j),'.','color',[0 0.9 0],'MarkerSize',20),        % Green
%                 
%             elseif MCPL(i,j) < 1.0*sigma(i)
%                 plot(BladeR(i,j),BladeC(i,j),'.','color',[0.95 0.95 0],'MarkerSize',20),    % Yellow
%                 
%             elseif MCPL(i,j) < 1.1*sigma(i)
%                 plot(BladeR(i,j),BladeC(i,j),'.','color',[1 0.5 0],'MarkerSize',20),        % Orange
%                 
%             else
%                 plot(BladeR(i,j),BladeC(i,j),'.','color',[1 0 0],'MarkerSize',20),          % Red
%             end
%         end
%     end
% 
%     HHcod = plot([BR,BR(end:-1:1),BR(1)],[-BCoD/2,BCoD(end:-1:1)/2,-BCoD(1)/2],'b-','LineWidth',2);
% 
%  
%% 
% % ------------------------------------- Pressure vs. chord for all sections
% figure, 
%     hold on, box on,
%     axis([-0.1 1.1 -1 1.1*max(abs([MCPL(:);MCPU(:)]))])
%     xlabel('x / c'               ,'Fontsize',16,'Fontname','Times')
%     ylabel('-CP'                 ,'Fontsize',16,'Fontname','Times')
%     title('Pressure distribution','Fontsize',16,'Fontname','Times')
%     set(gca,'Fontsize',14,'Fontname','Times')    
%         
%     for i = 1:Mp
%         plot(xp,MCPU(i,:),'-',xp,MCPL(i,:),'-')
%     end
% 
%     legend('-CPUpper','-CPLower')

% %% 
% % ------------------------------------- Pressure vs. chord for section i
% i = 4;
% 
% figure, 
%     hold on, box on,
%     axis([-0.1 1.1 -1 1.1*max(abs([MCPL(:);MCPU(:)]))])
%     xlabel('x / c'               ,'Fontsize',16,'Fontname','Times')
%     ylabel('-CP'                 ,'Fontsize',16,'Fontname','Times')
%     title('Pressure distribution','Fontsize',16,'Fontname','Times')
%     set(gca,'Fontsize',14,'Fontname','Times')    
%         
%     plot(xp,MCPU(i,:),'-b',xp,MCPL(i,:),'b-')
%     
%     plot([0 1],sigma(i)*[1 1],'k')
%     
%     plot([0 1],(CL(i)/1.8 + 2*t0oc(i))*[1 1],'r')
% 
%     legend('-CPUpper','-CPLower')
% 
% 
% 
% %%
% sigma ./ (CL'/1.8 + 2*t0oc')
% %%
% maxcp ./ (CL'/1.8 + 2*t0oc')
% 
% %%
% maxcp ./ sigma
% 
% %%
% i=1;
% 
% % trapz(xp,MCPU(i,:))
% 
% (CL(i)/1.8 + 2*t0oc(i))
% 
% %sigma(i)
% 
% (CL(i)/1.8 + 2*t0oc(i)) / maxcp(i)

