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
% =================================================== Make_Reports Function
%
% Make graphical and text reports for a propeller design
%
% -------------------------------------------------------------------------

function [] = Make_Reports(pt)

if isfield(pt,'i') & ~isfield(pt,'input' ), pt.input  = pt.i; end
if isfield(pt,'d') & ~isfield(pt,'design'), pt.design = pt.d; end

global Plots PlotPanels;

% ------------------------------------------------- Check for script or GUI
if isfield(pt.input,'GUI_flag')
    GUI_flag = pt.input.GUI_flag;
else
    GUI_flag = 0;
end


% -------------------------------------------------- Unpack input variables
input = pt.input;


if isfield(pt,'date'), Date_string = pt.date; else  Date_string = date; end

if     isfield(pt      ,'filename'),  filename =       pt.filename;
elseif isfield(pt      ,'name'),      filename =       pt.name;
elseif isfield(pt.input,'filename'),  filename = pt.input.filename; 
else                                  filename = 'OpenProp';
end


% % ----------------------------------------------- Execute directory changes
% currentdir  = pwd;
% parentdir   = ['../OpenProp_saved_files/' filename '/'];
% cd(parentdir);


% --------------------------------------------------------- Required inputs
Z = input.Z;           % number of blades

% --------------------------------------------------------- Geometry inputs
if isfield(input,  'Np'), Np   = input.Np;    else  Np   = 20;  end
if isfield(input,  'Vs'), Vs   = input.Vs;    else  Vs   = 1;   end % m/s
if isfield(input,   'R'), R    = input.R;     else  R    = 1;   end % m
if isfield(input,'Rhub'), Rhub = input.Rhub;  else  Rhub = 0.2; end % m
if isfield(input,'rho'),  rho  = input.rho;   else rho   = 1000;end % kg/m^3


% If propeller geometry or inflow is not given, assume propeller 4118 form 
% with uniform axial inflow, no swirl inflow, and section CD == 0.008.
% Note, the real propeller 4118 has CoD == 0.001 and t0oc == 0 at the tip.
if isfield(input,'XR'), 
    XR  = input.XR;
    X1  =  ones(size(XR));
    X0  = zeros(size(XR));
    
    if isfield(input,'XCoD'),  XCoD  = input.XCoD;  else  XCoD = X0; end
    if isfield(input,'t0oc0'), t0oc0 = input.t0oc0; else t0oc0 = X0; end  
    if isfield(input,'XVA'  ),   XVA = input.XVA;   else   XVA = X1; end
    if isfield(input,'XVT'  ),   XVT = input.XVT;   else   XVT = X0; end
else 
    XR    = [];    % r/R
    X1    = 1;
    X0    = 0;
    XCoD  = [];          % c/D
    t0oc0 = [];          % t0/c          
    XVA   = [];                                      % Va/Vs
    XVT   = [];                                      % Vt/Vs
end

if isfield(input,'XCD'), XCD = input.XCD; 
    if length(XCD) == 1, XCD =   XCD*X1;  end
else                     XCD = 0.008*X1;    
end
    
% ------------------------------------------------------------------- Flags
Propeller_flag = input.Propeller_flag; % 0 == turbine, 1 == propeller
Viscous_flag   = input.Viscous_flag;   % 0 == viscous forces off (CD = 0), 
                                       % 1 == viscous forces on
                                       
% 0 == no hub, 1 == hub
if isfield(input,'Hub_flag'),  Hub_flag = input.Hub_flag;  
                        else   Hub_flag = 1;  end
                                       
% 0 == no duct, 1 == duct
if isfield(input,'Duct_flag'),  Duct_flag = input.Duct_flag;  
                         else   Duct_flag = 0;  end

% 0 == do not optimize chord lengths, 1 == optimize chord lengths
if isfield(input,'Chord_flag'), Chord_flag = input.Chord_flag;  
                          else  Chord_flag = 0;  end

% 0 == do not display plots, 1 == display plots
if isfield(input,'Plot_flag'),  Plot_flag = input.Plot_flag;  
                         else   Plot_flag = 0;  end 

% 0 == Horseshoe(...,Wrench(...)), 1 == Wake_Horseshoe(...)      
if isfield(input,'Wake_flag'),  Wake_flag = input.Wake_flag;  
                         else   Wake_flag = 0;  end

% ---------------------------------------------------------- Propeller_flag
if Propeller_flag == 1
    Js = input.Js;
    L = pi/Js;
    
    % CTdes == desired total thrust coefficient == CTPdes + CTDdes
    if     isfield(input,'THRUST'),
            CTdes = input.THRUST / (0.5*rho*Vs^2*pi*R^2);
             
    elseif isfield(input,'CTDES'), 
            CTdes = input.CTDES; 
            
    elseif isfield(input,'CT'), 
            CTdes = input.CT; 
            
    elseif isfield(input,'KTDES'),
            CTdes = input.KTDES * (8/pi)/Js^2;
            
    elseif isfield(input,'KT'),
            CTdes = input.KT    * (8/pi)/Js^2;
    else 
            CTdes = 0;
    end
    
else
    L     = input.L;
    Js    = pi/L;
    CTdes = 0;     
end

% ------------------------------------------------------------ Viscous_flag 
if Viscous_flag == 0
    XCD   = X0;
    CDoCL = 0;
end
% --------------------------------------------------------------- Duct_flag
if Duct_flag == 1
    if     isfield(input,'TAU'),      TAU   = input.TAU; else TAU = 1; end  % thrust ratio == propeller thrust / total thrust 
    
    if     isfield(input,'Rduct_oR'), Rduct_oR = input.Rduct_oR;  % duct radius
    elseif isfield(input,'Rduct')     Rduct_oR = input.Rduct /R; 
    else                              Rduct_oR = 1;             
    end
    if     isfield(input,'Cduct_oR'), Cduct_oR = input.Cduct_oR;  % duct chord length
    elseif isfield(input,'Cduct')     Cduct_oR = input.Cduct /R; 
    else                              Cduct_oR = 1;             
    end
    if     isfield(input,'Xduct_oR'), Xduct_oR = input.Xduct_oR;
    elseif isfield(input,'Xduct')     Xduct_oR = input.Xduct /R;  % duct axial position downstream
    else                              Xduct_oR = 0;               % i.e. duct mid-chord centered at propeller          
    end
    if     isfield(input,'CDd'),      CDd   = input.CDd; else CDd = 0.008; end  % duct drag coefficient
    
else
    TAU     = 1; 
    Rduct   = 0;
    Cduct   = 0;
    CDd     = 0;
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ---------------------------------------------------- Computational inputs
if isfield(input,'HUF' ), HUF  = input.HUF;   else  HUF  = 0;    end
if isfield(input,'TUF' ), TUF  = input.TUF;   else  TUF  = 0;    end
if isfield(input,'Rhv' ), Rhv  = input.Rhv;   else  Rhv  = 0.5;  end

% -------------------------------------------------------------------------
% ------------------------------------------------------- Cavitation inputs
if isfield(input,'dVs'),  dVs = input.dVs;  else dVs  = 0.30;  end % m/s
if isfield(input,'H'  ),  H   = input.H;    else H    = 3.048; end % m
if isfield(input,'g'  ),  g   = input.g;    else g    = 9.81;  end % m/s^2
if isfield(input,'Patm'), Patm= input.Patm; else Patm = 101325;end % Pa
if isfield(input,'Pv'),   Pv  = input.Pv;   else Pv   = 2500;  end % Pa


% ------------------------------------------------- Unpack design variables
CT              = pt.design.CT;
CQ              = pt.design.CQ;
CP              = pt.design.CP;
VMIV            = pt.design.VMIV;
if Propeller_flag == 1
KT              = pt.design.KT;
KQ              = pt.design.KQ;
EFFY            = pt.design.EFFY;
end
RC              = pt.design.RC;
Mp              = length(RC);
RV              = pt.design.RV;
G               = pt.design.G';
VAC             = pt.design.VAC;
VTC             = pt.design.VTC;
UASTAR          = pt.design.UASTAR;
UTSTAR          = pt.design.UTSTAR;
TANBC           = pt.design.TANBC;
TANBIC          = pt.design.TANBIC;
CoD             = pt.design.CoD;
CL              = pt.design.CL;
CD              = pt.design.CD;
if Duct_flag == 1    
    TAU             = pt.design.TAU;
    XdRING          = pt.design.XdRING;
    GdRING          = pt.design.GdRING;
    UADIF           = pt.design.UADIF;
    Gd              = pt.design.Gd;
    UARING          = pt.design.UARING;
    URRING          = pt.design.URRING;
else
    TAU             = 1;
    XdRING          = 1;
    GdRING          = 0;
    UADIF           = 0*RC;
    Gd              = 0;
    UARING          = 0;
    URRING          = 0;
    Rduct_oR        = 1;
end
Beta_c          = atand(pt.design.TANBC);  % [deg]
BetaI_c         = atand(pt.design.TANBIC); % [deg]



% -------------------------------------------------------------------------
D       = 2*R;    % [m]
Dhub    = 2*Rhub; % [m]
Rhub_oR = Rhub/R;
N       = 60*Vs/(Js*D); % [RPM]
n       = N/60;         % [rev/s]

% --------------------------------------------- Create Graphical Report
if GUI_flag
    
    set(0,'CurrentFigure',Plots);
    h = axes('parent',PlotPanels(7));
    axes(h);
    hold on;
    
    plot(RC,G,'b.-','LineWidth',2,'MarkerSize',16);
    
    xlabel('r/R','FontSize',16,'FontName','Times');
    ylabel('Gamma / 2\pi R Vs','FontSize',16,'FontName','Times');
    set(gca,     'FontSize',14,'FontName','Times');
    grid on; box on,
    
    ylimits = get(gca,'Ylim');
    
    if abs(ylimits(2)) > abs(ylimits(1))
        
        set(gca,'Ylim',[0 ylimits(2)]);
    else
        set(gca,'Ylim',[ylimits(1) 0]);
    end
    % ---
    
    set(0,'CurrentFigure',Plots);
    h = axes('parent',PlotPanels(8));
    axes(h);
    hold on;
    
    plot(RC,VAC,'-b','LineWidth',2)
    plot(RC,VTC,'--b','LineWidth',2)
    plot(RC,UASTAR,'-.r','LineWidth',2)
    plot(RC,UTSTAR,':r','LineWidth',2);
    
    legend('Va/Vs','Vt/Vs','Ua*/Vs','Ut*/Vs');
    
    xlabel('r/R','FontSize',16,'FontName','Times');
    ylabel('velocity','FontSize',16,'FontName','Times');
    set(gca,     'FontSize',14,'FontName','Times');
    grid on; box on,
    
    
    % ---
    
    set(0,'CurrentFigure',Plots);
    h = axes('parent',PlotPanels(9));
    axes(h);
    hold on;
    
    plot(RC,Beta_c,'--b','LineWidth',2)
    plot(RC,BetaI_c,'-r','LineWidth',2)
    
    
    legend('','BetaI');
    
    xlabel('r/R','FontSize',16,'FontName','Times');
    ylabel('Inflow angle  [deg]','FontSize',16,'FontName','Times');
    set(gca,     'FontSize',14,'FontName','Times');
    grid on; box on,
    
    ylimits = get(gca,'Ylim');
    set(gca,'Ylim',[0 ylimits(2)]);
    
    % ---
    
    set(0,'CurrentFigure',Plots);
    h = axes('parent',PlotPanels(10));
    axes(h);
    hold on;
    
    
    XXRC  = Rhub_oR + (1-Rhub_oR)*(sin((0:60)*pi/(2*60)));
    XXCoD = InterpolateChord(RC, CoD,XXRC);
    
    plot(XXRC, XXCoD,'b','LineWidth',2);
    plot(XXRC,-XXCoD,'b','LineWidth',2);
    
    plot(RC, CoD,'b.','MarkerSize',16);
    plot(RC,-CoD,'b.','MarkerSize',16);
    
    xlabel('r/R','FontSize',16,'FontName','Times');
    ylabel('c/R','FontSize',16,'FontName','Times');
    set(gca,     'FontSize',14,'FontName','Times');
    grid on; box on,
    
else
    
    Fig1_S = figure('units','normalized','position',[.01 .06 .4 .3],...
        'name','Graphical Report','numbertitle','off');
    
    subplot(2,2,1);
    plot(RC,G);
    xlabel('r/R');     ylabel('Gamma / 2piRVs');      grid on;
    
    if Propeller_flag == 1
        TitleString=strcat('J='  ,num2str(Js,        '%10.3f'),...
                           '; Ct='  ,num2str(CT,  '%10.3f'),...
                           '; Cq='  ,num2str(CQ,  '%10.3f'),...
                           '; Kt='  ,num2str(KT,  '%10.3f'),...
                           '; Kq='  ,num2str(KQ,  '%10.3f'),...
                           '; \eta=',num2str(EFFY,'%10.3f'),...
                           '; \tau=',num2str(TAU, '%10.3f'));
    else
        TitleString=strcat('J='     ,num2str(Js,  '%10.3f'),...
                           '; Ct='  ,num2str(CT,  '%10.3f'),...
                           '; Cq='  ,num2str(CQ,  '%10.3f'),...
                           '; \tau=',num2str(TAU, '%10.3f'));        
    end
    title(TitleString);
    
    subplot(2,2,2);
    plot(RC,VAC,'-b',RC,VTC,'--b',RC,UASTAR,'-.r',RC,UTSTAR,':r');
    xlabel('r/R');   grid on;    legend('Va/Vs','Vt/Vs','Ua*/Vs','Ut*/Vs');
    
    subplot(2,2,3);
    plot(RC,Beta_c,'--b',RC,BetaI_c,'-r');
    xlabel('r/R');   ylabel('Degrees');   grid on;  legend('Beta','BetaI');
    
    subplot(2,2,4);
    plot(RC,CoD);
    xlabel('r/R');   ylabel('c/D');       grid on;
    
end

% ------------------------------------------------- Create text reports

Date_string = datestr(now,31);      % Date and time to print on reports

% --------------------------------------------- Make OpenProp_Input.txt
filename_input = strcat(filename,'_Input.txt');
fid = fopen(filename_input,'wt');   % 1-Aug-2013 BEPPS: 'w' changed to 'wt'  so newline character '\n' appears properly on Windows machines

fprintf(fid,'\t\t\t\t\t %s \n\n',filename_input);
fprintf(fid,'\t\t\t\t\t OpenProp Input Table\n\n');
fprintf(fid,'Date and time: %s\n\n',Date_string);

fprintf(fid,'%.0f \tNumber of Vortex Panels over the Radius\n',         Mp);
fprintf(fid,'%.0f \tHub Image Flag: 1=YES, 0=NO\n',                     Hub_flag);
fprintf(fid,'%.0f \tDuct Flag:      1=YES, 0=NO\n',                     Duct_flag);
fprintf(fid,'%.3f \tDuct Diameter\n',                                   Rduct_oR*2*R);
fprintf(fid,'%.1f \tHub Vortex Radius/Hub Radius\n',                    Rhv);
fprintf(fid,'%.0f \tNumber of Blades\n',                                Z);
fprintf(fid,'%.3f \tAdvance Coefficient Based on Ship Speed, Js\n',     Js);
fprintf(fid,'%.3f \tDesired Thrust Coefficient, Ct\n',                  CTdes);
fprintf(fid,'%.3f \tDesired Thrust Ratio, tau\n',                       TAU);
fprintf(fid,'%.3f \tDuct Section Drag Coefficient, CDd\n',              CDd);
fprintf(fid,'%.0f \tHub Unloading Factor: 0 = optimum\n',               HUF);
fprintf(fid,'%.0f \tTip Unloading Factor: 1 = Reduced Loading\n',       TUF);
fprintf(fid,'r/R  \t  C/D  \t   XCD\t    Va/Vs  Vt/Vs\n');

N_R0 = length(XR);

for i = 1:N_R0
    fprintf(fid,'%6.5f  %6.5f  %6.5f  %6.2f  %6.4f\n',XR(i),XCoD(i),XCD(i),XVA(i),XVT(i));
end

fprintf(fid,'\nr/R \t [ ], input radial position / propeller radius.\n');
fprintf(fid,'c/D \t [ ], input section chord-length / propeller diameter.\n');
fprintf(fid,'Cd \t [ ], input section drag coefficient.\n');
fprintf(fid,'Va \t [ ], input axial inflow velocity / ship velocity.\n');
fprintf(fid,'Vt \t [ ], input tangential inflow velocity / ship velocity.\n');

fclose(fid);

% -------------------------------------------- Make OpenProp_Output.txt
filename_output = strcat(filename,'_Output.txt');
fid = fopen(filename_output,'wt');   % 1-Aug-2013 BEPPS: 'w' changed to 'wt'  so newline character '\n' appears properly on Windows machines

fprintf(fid,'\t\t\t\t\t %s \n\n',filename_output);
fprintf(fid,'\t\t\t\t\t OpenProp Output Table\n\n');
fprintf(fid,'Date and time: %s\n\n',Date_string);

fprintf(fid,'Js \t= %5.4f\n' , Js);
fprintf(fid,'Ct \t= %5.4f\n' , CT);
fprintf(fid,'Cq \t= %5.4f\n' , CQ);
fprintf(fid,'Cp \t= %5.4f\n' , CP);
fprintf(fid,'VMIV \t= %5.4f\n',VMIV);
if Propeller_flag == 1
fprintf(fid,'Kt \t= %5.4f\n' , KT);
fprintf(fid,'Kq \t= %5.4f\n' , KQ);
fprintf(fid,'Eff \t= %5.4f\n',EFFY);
end
fprintf(fid,'Tau \t= %5.4f\n',TAU);
fprintf(fid,'Duct Circulation \t= %5.4f\n',Gd);


%     fprintf(fid,'%5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s\n',...
%                 'r/R','G','Va','Vt','Ua','Ut','Beta','BetaI','c/D','Cd');
fprintf(fid,'\nOutput at the control points for the propeller \n\n');
fprintf(fid,'r/R\t\t G\t\t\t Va\t\t Vt\t\t Ua\t\t Ua(ring)\t Ut\t\t Beta\t BetaI\t c/D\t Cd\n');

for i = 1:length(RC)
    fprintf(fid,'%5.5f  %5.6f  %5.5f  %5.4f  %5.5f  %5.5f  %5.5f  %5.3f  %5.3f  %5.5f  %5.5f\n',...
            RC(i),G(i),VAC(i),VTC(i),UASTAR(i),Gd*UADIF(i),UTSTAR(i),Beta_c(i),BetaI_c(i),CoD(i),CD(i));
end

if Duct_flag == 1
    fprintf(fid,'\nOutput on the duct ring vortices \n\n');
    fprintf(fid,'X/R\t\t G\t\t\t UA/VS\t UR/VS\n');
    for i = 1:length(XdRING)
        fprintf(fid,'%5.5f  %5.6f  %5.5f  %5.4f\n',XdRING(i),Gd*GdRING(i),UARING(i),URRING(i));
    end
else fprintf(fid,'\nThe propeller does not have a duct.\n\n');
end

fprintf(fid,'\nJs \t [ ], advance coefficient.\n');
fprintf(fid,'Ct \t [ ], required thrust coefficient.\n');
fprintf(fid,'Cp \t [ ], power coefficient. Cp = Cq*pi/J.\n');
fprintf(fid,'Kt \t [ ], thrust coefficient. Kt = Ct*Js^2*pi/8.\n');
fprintf(fid,'Kq \t [ ], torque coefficient. Kq = Cq*Js^2*pi/16.\n');
fprintf(fid,'VMIV \t [ ], volumetric mean inflow velocity / ship velocity.\n');
fprintf(fid,'Eff \t [ ], efficiency = Ct*VMIV/Cp.\n');
fprintf(fid,'Tau \t [ ], thrust ratio = propeller thrust / total thrust.\n');

fprintf(fid,'\nr/R \t [ ], radial position of control points / propeller radius.\n');
fprintf(fid,'G  \t [ ], section circulation / 2*pi*R.\n');
fprintf(fid,'Va \t [ ], axial inflow velocity / ship velocity.\n');
fprintf(fid,'Vt \t [ ], tangential inflow velocity / ship velocity.\n');
fprintf(fid,'Ua \t [ ], induced axial velocity / ship velocity.\n');
fprintf(fid,'Ut \t [ ], induced tangential velocity / ship velocity.\n');
fprintf(fid,'beta \t [deg], flow angle.\n');
fprintf(fid,'betaI \t [deg], hydrodynamic Pitch angle.\n');
fprintf(fid,'c/D \t [ ], section chord-length / propeller diameter.\n');
fprintf(fid,'Cd \t [ ], section drag coefficient.\n');

fprintf(fid,'\nX/R \t [ ], axial location of duct vortex rings / propeller radius.\n');
fprintf(fid,'G  \t [ ], duct vortex ring circulation / 2*pi*R.\n');
fprintf(fid,'UA/VS \t [ ], axial inflow induced by propeller / ship velocity.\n');
fprintf(fid,'UR/VS \t [ ], radial inflow induced by propeller / ship velocity.\n');


fclose(fid);

% ----------------------------- Calculate Propeller Section Performance
w = 2*pi*n;                                              % angular velocity [rad/s]
dV = dVs*Vs;

for k = 1:Mp
    Vstar(k)  = sqrt((VAC(k)+UASTAR(k))^2 + (w*R*RC(k)/Vs+VTC(k)+UTSTAR(k))^2)*Vs; % inflow velocity [m/s]

    Gamma(k)  = G(k)*2*pi*R*Vs;                          % circulation [m^2/s]

    CL(k)     = 2*Gamma(k) / (Vstar(k)*CoD(k)*D);        % lift coefficient

    Sigma(k)  = (101000+rho*9.81*(H-RC(k)*R)-2500)/...   % cavitation number
                (rho*Vstar(k)^2/2);

    dBetai(k) = atand((TANBIC(k)*w*RC(k)*R+dV)/(w*RC(k)*R))... % inflow angle variation [deg]
               -atand((TANBIC(k)*w*RC(k)*R-dV)/(w*RC(k)*R));
end


% --------------------------------------- Make OpenProp_Performance.txt
filename_performance = strcat(filename,'_Performance.txt');
fid = fopen(filename_performance,'wt');    % 1-Aug-2013 BEPPS: 'w' changed to 'wt'  so newline character '\n' appears properly on Windows machines

fprintf(fid,'\t\t\t\t\t %s \n\n',filename_performance);
fprintf(fid,'\t\t\t\t\t OpenProp Performance Table\n\n');
fprintf(fid,'Date and time: %s\n\n',Date_string);

fprintf(fid,' r/R\t V*\t beta\t betai\t Gamma\t CL\t Sigma\t dBetai\n');
for k = 1:Mp
    fprintf(fid,'%.3f\t %.2f\t %.2f\t %.2f\t %.4f\t %.3f\t %.3f\t %.2f\n'...
    ,RC(k),Vstar(k),Beta_c(k),BetaI_c(k),Gamma(k),CL(k),Sigma(k),dBetai(k));
end

fprintf(fid,'\nr/R \t [ ], radial position of control points / propeller radius.\n');
fprintf(fid,'V* \t [m/s], total inflow velocity.\n');
fprintf(fid,'beta \t [deg], undisturbed flow angle.\n');
fprintf(fid,'betai \t [deg], hydrodynamic Pitch angle.\n');
fprintf(fid,'Gamma \t [m^2/s], vortex sheet strength.\n');
fprintf(fid,'CL \t [ ], section lift coefficient.\n');
fprintf(fid,'Sigma \t [ ], cavitation number.\n');
fprintf(fid,'d_alpha  [deg], inflow variation bucket width.\n');

fclose(fid);


% cd(currentdir);

end