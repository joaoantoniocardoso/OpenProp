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
% =================================================== Make_ParReports Function
%
% Make graphical and text reports for a parametric propeller design study
%
% -------------------------------------------------------------------------

function [] = Make_ParReports(pt)

% ------------------------------------------------- Unpack design variables
parinput  = pt.parinput;
paroutput = pt.paroutput;

% '------ Performance inputs ------'
Zall    = parinput.Z;           % [numberZ x 1], [ ]   number of blades
Nall    = parinput.N;           % [numberN x 1], [RPM] rotation rate
Dall    = parinput.D;           % [numberD x 1], [m]   propeller diameter

% Vs      = parinput.Vs;          % [1 x 1]
% THRUST  = parinput.THRUST;       % [1 x 1]
% % '------ Geometry inputs ------'
% Mp      = parinput.Mp;          % [1 x 1]
% Rhub    = parinput.Rhub;        % [1 x 1]
% XR      = parinput.XR;          % [length(XR) x 1]
% XVA     = parinput.XVA;         % [length(XR) x 1]
% XVT     = parinput.XVT;         % [length(XR) x 1]
% XCD     = parinput.XCD;         % [length(XR) x 1]
% XCoD    = parinput.XCoD;        % [length(XR) x 1]
% % '------ Computational inputs ------'
% 
% % '------ Cavitation inputs ------'
% rho     = parinput.rho;         % [1 x 1], [kg/m^3], density of seawater
% dVs     = parinput.dVs;         % [1 x 1]
% H       = parinput.H;           % [1 x 1]
% % '------ Duct inputs ------'
% TAU     = parinput.TAU;         % [1 x 1]
% Rduct_oR= parinput.Rduct_oR;    % [1 x 1]
% CDd     = parinput.CDd;         % [1 x 1]
% CDoCL   = parinput.CDoCL;


% '------ Performance metrics ------'
% indexd by (iZ,iN,iD)
CT      = paroutput.CT;
CQ      = paroutput.CQ;
CP      = paroutput.CP;
KT      = paroutput.KT;
KQ      = paroutput.KQ;
EFFY    = paroutput.EFFY;
VMIV    = paroutput.VMIV;
TAU     = paroutput.TAU;
CTP     = paroutput.CTP;
CTH     = paroutput.CTH;
CTD     = paroutput.CTD;
KTP     = paroutput.KTP;


% ------------------------------------------------- Plot efficiency results
Fig_P = figure('units','normalized','position',[0.01 .06 .4 .3],...
               'name','Efficiency','numbertitle','off');
for iZ = 1:length(Zall) % for each number of blades
    if length(Zall)==4
       subplot(2,2,iZ);
    else
       subplot(length(Zall),1,iZ);
    end
    hold on,

    % - plot EFFY vs. D for different N
    for iN = 1:length(Nall)
        plot(Dall(:),reshape(EFFY(iZ,iN,:),length(Dall),1),'.-','MarkerSize',20);
    end
    
    for iN = 1:length(Nall)
        str_legend(iN) = {strcat(num2str(Nall(iN)),' RPM')};
    end
    legend(str_legend,'location','southwest');    
    grid on;
    axis([Dall(1) Dall(end) 0 1]),
    box on,
    xlabel('Propeller Diameter (m)');             
    ylabel('Efficiency');
    title_prefix = {'Number of Blades: '};
    
    title(strcat(title_prefix,num2str(Zall(iZ))))
end

    
    
    