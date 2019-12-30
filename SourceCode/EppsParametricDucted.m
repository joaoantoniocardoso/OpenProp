% =========================================================================
% ================================================= EppsParametric Function
%
% Perform a parametric design study using the EppsOptimizer function for
% the propeller/turbine design optimization.
%
% Author: Brenden Epps
% Last modified: 12/28/2010 by Brenden Epps
% -------------------------------------------------------------------------


% USE THIS CODE FOR PARAMETRIC STUDY VARYING....
%
% Z, N, TAU
%
% ... with fixed D ...

          
function paroutput = EppsParametricDucted(parinput)
         
% -------------------------------------------------- Unpack input variables
Zall    = parinput.Z;           % [numberZ x 1], [ ]   number of blades
Nall    = parinput.N;           % [numberN x 1], [RPM] rotation rate
Tall    = parinput.TAU;         % [numberT x 1], [ ]   thrust ratio == propeller thrust / total thrust



if isfield(parinput,'Duct_flag'),  Duct_flag = parinput.Duct_flag;  else   Duct_flag = 0;  end  % 0 == no duct, 1 == duct
  
if Duct_flag == 1

    D = parinput.D;           % [numberT x 1], [m]   propeller diameter
    
    parinput.R  = D/2;

    
    if isfield(parinput,'Vs' ), Vs = parinput.Vs;   else  Vs = 1;   end % m/s

else
    disp('ERROR: Only use this code with ducted propellers.')
    return    
end


% ----------------------------------- Initialize output performance metrics
numberZ = length(Zall);
numberN = length(Nall);
numberT = length(Tall);

JS   = zeros(numberZ,numberN,numberT);
CT   = zeros(numberZ,numberN,numberT);
CQ   = zeros(numberZ,numberN,numberT);
CP   = zeros(numberZ,numberN,numberT);
KT   = zeros(numberZ,numberN,numberT);
KQ   = zeros(numberZ,numberN,numberT);
EFFY = zeros(numberZ,numberN,numberT);
TAU  = zeros(numberZ,numberN,numberT);
CTP  = zeros(numberZ,numberN,numberT);
CTH  = zeros(numberZ,numberN,numberT);
CTD  = zeros(numberZ,numberN,numberT);
KTP  = zeros(numberZ,numberN,numberT);
VMIV = zeros(numberZ,numberN,numberT);

propNo = 1;                       % propeller number
Nprops = numberZ*numberN*numberT; % total number of props to analyze


for iZ = 1:numberZ
    for iN = 1:numberN
        for iT = 1:numberT
            if     propNo < 10
                disp(['Optimizing propeller   ',num2str(propNo),' of ',num2str(Nprops)])
            elseif propNo < 100
                disp(['Optimizing propeller  ' ,num2str(propNo),' of ',num2str(Nprops)])
            else
                disp(['Optimizing propeller '  ,num2str(propNo),' of ',num2str(Nprops)])
            end
             
            % ------------------------------------- Modify input parameters
            Z   = Zall(iZ);
            N   = Nall(iN);
            TAU = Tall(iT);
            
            Js = Vs/((N/60)*D);
            
            parinput.Z   = Z;
            parinput.TAU = TAU;
            parinput.Js  = Js;
            parinput.L   = pi/Js;
            
            % --------------------------------- Perform design optimization
            design = EppsOptimizer(parinput);
            
            
            % ----------------------------------- Record output performance
              JS(iZ,iN,iT) = Js;
              CT(iZ,iN,iT) = design.CT;
              CQ(iZ,iN,iT) = design.CQ;
              CP(iZ,iN,iT) = design.CP;
              KT(iZ,iN,iT) = design.KT;
              KQ(iZ,iN,iT) = design.KQ;
            EFFY(iZ,iN,iT) = design.EFFY;
            VMIV(iZ,iN,iT) = design.VMIV;
             TAU(iZ,iN,iT) = design.TAU;
             CTP(iZ,iN,iT) = design.CTP;
             CTH(iZ,iN,iT) = design.CTH;
             CTD(iZ,iN,iT) = design.CTD;
             KTP(iZ,iN,iT) = design.KTP;
            
            propNo = propNo + 1;
        end
    end
end
% -------------------------------------- Pack up output performance metrics
paroutput.part1 = '------ Performance metrics ------';
paroutput.index = '(iZ,iN,iT)';
paroutput.Z     = Zall;
paroutput.N     = Nall;
paroutput.Tall  = Tall;
paroutput.D     = D;
paroutput.Js    = JS;
paroutput.CT    = CT;
paroutput.CQ    = CQ;
paroutput.CP    = CP;
paroutput.KT    = KT;
paroutput.KQ    = KQ;
paroutput.EFFY  = EFFY;
paroutput.VMIV  = VMIV;
paroutput.TAU   = TAU;
paroutput.CTP   = CTP;
paroutput.CTH   = CTH;
paroutput.CTD   = CTD;
paroutput.KTP   = KTP;
% ============================================== END EppsOptimizer Function
% =========================================================================