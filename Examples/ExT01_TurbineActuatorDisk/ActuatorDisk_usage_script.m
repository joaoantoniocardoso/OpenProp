% -------------------------------------------------------------------------
% Turbine "actuator disk" example:
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
addpath ./SourceCode ../SourceCode ../../SourceCode
clr,

% Load inputs:
ActuatorDisk_input,
pause(0.01),


% -------------------------------------------------------------------------
% Example grid refinement study:
i = pt.input;

i.Mp = 20; % vortex panels

     d20 = EppsOptimizer(i);

     
i.Mp = 80; % vortex panels

     d80 = EppsOptimizer(i);



fig; 
    plot(d80.RC,d80.G,'-b')
    plot(d20.RC,d20.G,'ro')
    xlabel('r/R'), ylabel('G')

disp('---------------------------------------')    
disp(['CP (20 panels)  : ',num2str(d20.CP)])
disp(['CP (80 panels)  : ',num2str(d80.CP)])
disp(['CP actuator disk: ',num2str(d80.CPBetz)])
  

% Moral: 
%   You need 80 panels to match CP with CPBetz, 
%   The circulation is the same with 20 or 80 panels
%   The code runs much faster with 20 panels.  
    

%%
% -------------------------------------------------------------------------
% Show several turbines designed for various L and CDoCL
clr,

% Load inputs:
ActuatorDisk_input    
   
input = pt.input;

input.Mp           = 80; % vortex panels
input.Viscous_flag = 1;  % viscous forces on



L       = [1 2 3 4 5 6 7 8 9 10];        
CDoCL   = [0 0.02 0.04];

[L,CDoCL]  = meshgrid(L,CDoCL);

CP     = zeros(size(L));
CPBetz = zeros(size(L)); % actuator disk with swirl power coefficient


for i = 1:size(L,1)
    for j = 1:size(L,2)
        
        input.L     =     L(i,j);
        input.XCD   = CDoCL(i,j) * input.XCLmax;
        
        design = EppsOptimizer(input);
        
            CP(i,j) = -design.CP;
        CPBetz(i,j) = -design.CPBetz;
    end
end


CPBetzLimit = 16/27; % Betz limit for infinite L and Z


fig;
    plot([0 10],CPBetzLimit*[1 1],'k--')
        
    for i = 1:size(L,1)
        plot(L(i,:),CPBetz(i,:),'-k')
        plot(L(i,:),    CP(i,:),'ro-')
    end 
    axis([0 10 0 0.6])
    
    title('Parametric study: turbines optimized with CD/CL = 0, 0.02, 0.04')

    xlabel('tip speed ratio')
    ylabel('power coefficient')    

    HB = plot([-1 -1],[-1 -2],'k');
    HD = plot([-1 -1],[-1 -2],'ro-'); 

    legend([HB,HD],'Actuator disk','Z = 100 blades','Location','SouthEast')
    
    