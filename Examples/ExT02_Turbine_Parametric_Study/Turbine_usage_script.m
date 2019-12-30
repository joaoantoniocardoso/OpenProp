% -------------------------------------------------------------------------
% Turbine design example:
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
addpath ./SourceCode ../SourceCode ../../SourceCode
clr,

% Load inputs:
% The only differences between this and the "actuator disk" example are that here:
%   Z  = 3 blades 
%   Mp = 20 panels
%
Turbine_input,
pause(0.01),

input = pt.input;

% -------------------------------------------------------------------------
% Show several turbines designed for various L
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

legend([HB,HD],'Actuator disk','Z = 3 blades','Location','SouthEast')


