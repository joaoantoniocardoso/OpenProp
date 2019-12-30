% -------------------------------------------------------------------------
% Turbine design example:
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
addpath ./SourceCode ../SourceCode ../../SourceCode
clr,

% Load inputs:
% The only differences between this and the "actuator disk" example are that here: 
%   Mp = 20 panels
%
Turbine_input,
pause(0.01),

input = pt.input;

CDoCL     = 0.01;
input.XCD = CDoCL * input.XCLmax;


% -------------------------------------------------------------------------
% Show several turbines designed for various L
L       = [0.125 0.25 0.5 1 2 3 4 5 6 7 8 9 10];        
Z       = [2 3 4 5];

[L,Z]  = meshgrid(L,Z);

CP     = zeros(size(L));
CPBetz = zeros(size(L)); % actuator disk with swirl power coefficient

for i = 1:size(L,1)
    for j = 1:size(L,2)
        
        input.L = L(i,j);
        input.Z = Z(i,j);
         
        design = EppsOptimizer(input);
        
        if design.converged == 1
                CP(i,j) = -design.CP;
            CPBetz(i,j) = -design.CPBetz;
        else
                CP(i,j) = NaN;
            CPBetz(i,j) = NaN;
        end
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
    
    xlabel('tip speed ratio')
    ylabel('power coefficient')
    
HB = plot([-1 -1],[-1 -2],'k');
HD = plot([-1 -1],[-1 -2],'ro-'); 

legend([HB,HD],'Actuator disk with CD/CL = 0.01','Turbine with 2, 3, 4, 5 blades','Location','SouthEast')


