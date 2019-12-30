% -------------------------------------------------------------------------
% Ducted turbine design example:
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
addpath ./SourceCode ../SourceCode ../../SourceCode
clr,

Ducted_Turbine_input,


%
% Perform design optimization
pt.d = EppsOptimizer(pt.i);

% fig; 
%     plot(pt.d.RC,-pt.d.G,'.-b')
%     axis([0 1 0 0.03])
%     xlabel('r/R'), ylabel('G')
% 
% 
% fig; 
%     plot(pt.d.RC,pt.d.CL)
%     xlabel('r/R'), ylabel('CL')
%%
% Create graphical and text reports
Make_Reports(pt)

%%
% close all,
pt.i.Plot_flag = 0;

% Analyze off-design states
LAMBDAall = [0 : 0.1 : 15];           % tip-speed ratio
pt.s = Analyze(pt,LAMBDAall)


fig;
    plot(pt.s.L,   -pt.s.CP  ,'.-b')
    
    plot(pt.i.L,   -pt.d.CP  ,'xr')
    
    xlabel('L'), ylabel('CP')
    axis([0 15 0 2])

    


%%
% Determine propeller geometry
pt.geometry = Geometry(pt)


