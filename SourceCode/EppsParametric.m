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
% ================================================= EppsParametric Function
% Last modified: 11/5/2011 by Brenden Epps
%
% Perform a parametric design study using the EppsOptimizer function for
% the propeller/turbine design optimization.
%
% -------------------------------------------------------------------------
          
function paroutput = EppsParametric(parinput)

% -------------------------------------------------------------------------
if isfield(parinput,'Propeller_flag')
    if Propeller_flag == 0
        disp('Sorry, EppsParametric is only set up for propellers...please update the code for turbines...')
        return
    end
end
% -------------------------------------------------------------------------


% -------------------------------------------------- Unpack input variables
Zall    = parinput.Z;           % [numberZ x 1], [ ]   number of blades
Nall    = parinput.N;           % [numberN x 1], [RPM] rotation rate
Dall    = parinput.D;           % [numberD x 1], [m]   propeller diameter

if isfield(parinput,'Vs' ), Vs = parinput.Vs;   else  Vs = 1;   end % m/s

if isfield(parinput,'Duct_flag'),  Duct_flag = parinput.Duct_flag;  else   Duct_flag = 0;  end  % 0 == no duct, 1 == duct
  
                         
% ----------------------------------- Initialize output performance metrics
numberZ = length(Zall);
numberN = length(Nall);
numberD = length(Dall);

JS     = zeros(numberZ,numberN,numberD);
KT     = zeros(numberZ,numberN,numberD);
KQ     = zeros(numberZ,numberN,numberD);
CT     = zeros(numberZ,numberN,numberD);
CQ     = zeros(numberZ,numberN,numberD);
CP     = zeros(numberZ,numberN,numberD);
CTH    = zeros(numberZ,numberN,numberD);
VMIV   = zeros(numberZ,numberN,numberD);
VMWV   = zeros(numberZ,numberN,numberD);
  EFFY = zeros(numberZ,numberN,numberD);
ADEFFY = zeros(numberZ,numberN,numberD);
QF     = zeros(numberZ,numberN,numberD);


if Duct_flag == 1
TAU  = zeros(numberZ,numberN,numberD);
CTD  = zeros(numberZ,numberN,numberD);
end    
    
    

propNo = 1;                       % propeller number
Nprops = numberZ*numberN*numberD; % total number of props to analyze


for iZ = 1:numberZ
    for iN = 1:numberN
        for iD = 1:numberD
            if     propNo < 10
                disp(['Optimizing propeller   ',num2str(propNo),' of ',num2str(Nprops)])
            elseif propNo < 100
                disp(['Optimizing propeller  ' ,num2str(propNo),' of ',num2str(Nprops)])
            else
                disp(['Optimizing propeller '  ,num2str(propNo),' of ',num2str(Nprops)])
            end
             
            % ------------------------------------- Modify input parameters
            Z = Zall(iZ);
            N = Nall(iN);
            D = Dall(iD);
            
            Js = Vs/((N/60)*D);
            
            parinput.Z  = Z;
            parinput.R  = D/2;
            parinput.Js = Js;
            parinput.L  = pi/Js;
            
            % --------------------------------- Perform design optimization
            design = EppsOptimizer(parinput);
            
            
            % ----------------------------------- Record output performance
              JS(iZ,iN,iD) = Js;
              KT(iZ,iN,iD) = design.KT;
              KQ(iZ,iN,iD) = design.KQ;
              CT(iZ,iN,iD) = design.CT;
              CQ(iZ,iN,iD) = design.CQ;
              CP(iZ,iN,iD) = design.CP;
             CTH(iZ,iN,iD) = design.CTH;
             
            VMIV(iZ,iN,iD) = design.VMIV;
            VMWV(iZ,iN,iD) = design.VMWV;
             
            EFFY(iZ,iN,iD) = design.EFFY;
          ADEFFY(iZ,iN,iD) = design.ADEFFY;
              QF(iZ,iN,iD) = design.QF;
              
            if Duct_flag == 1
             TAU(iZ,iN,iD) = design.TAU;
             CTD(iZ,iN,iD) = design.CTD;
            end
            
              
            propNo = propNo + 1;
        end
    end
end
% -------------------------------------- Pack up output performance metrics
paroutput.part1 = '------   Free parameters   ------';
paroutput.Z     = Zall;
paroutput.N     = Nall;
paroutput.D     = Dall;
paroutput.part2 = '------ Performance metrics (iZ,iN,iD) ------';
paroutput.Js    = JS;
paroutput.KT    = KT;
paroutput.KQ    = KQ;
paroutput.CT    = CT;
paroutput.CQ    = CQ;
paroutput.CP    = CP;
paroutput.CTH   = CTH;
paroutput.EFFYo = CT ./ CP;              %   open-water   efficiency

paroutput.VMIV  = VMIV;
paroutput.Ja    = JS .* VMIV;
paroutput.EFFY  = EFFY;

paroutput.VMWV  = VMWV;
            
if Duct_flag == 1
paroutput.TAU   = TAU;
paroutput.CTD   = CTD;
end
% ============================================== END EppsOptimizer Function
% =========================================================================