% -------------------------------------------------------------------------
% Run the following script to use OpenProp via the command line.
%
% NOTE: First set Matlab's working directory to 
%           /Examples/ExP01_PropellerActuatorDisk
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
addpath ./SourceCode ../SourceCode ../../SourceCode
clr,

% Load inputs:
ActuatorDisk_input


% Perform design optimization
pt.design = EppsOptimizer(pt.input)


fig; 
    plot(pt.design.RC,pt.design.G,'.-b','MarkerSize',20)
    plotaxes,
    xlabel('r/R'), ylabel('G')

    
pause(0.01)


%%
% -------------------------------------------------------------------------
% Show several propellers designed for various Js and CT
clr,

% Load inputs:
ActuatorDisk_input


JSall    = 0.1:0.2:0.9;
CTDESall = 0.1:1  :2.1;

[Js,CT]  = meshgrid(JSall,CTDESall);
EFFY     = zeros(size(CT));

ADEFFY   = 2./(1+sqrt(1+CT)); % actuator disk efficiency


for i = 1:size(CT,1)
    for j = 1:size(CT,2)
        
        pt.input.Js    = Js(i,j);
        pt.input.CTDES = CT(i,j);
        
        design = EppsOptimizer(pt.input);
        
          EFFY(i,j) = design.EFFY;
    end
end


figure, hold on
    for i = 1:size(CT,1)
        plot(Js(i,:),  EFFY(i,:),'.-b','MarkerSize',20)
        plot([0 1],  ADEFFY(i,1)*[1 1],'--r')
    end
    axis([0 1 0 1])
    box on, grid on
    xlabel('Js'), ylabel('EFFY')