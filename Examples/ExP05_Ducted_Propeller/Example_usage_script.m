% -------------------------------------------------------------------------
% Ducted propeller design example.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
addpath ./SourceCode ../SourceCode ../../SourceCode
clr,

% Load inputs by running the Example_input script:
Example_input


% Perform design optimization
pt.design = EppsOptimizer(pt.input)

fig; 
    plot(pt.design.RC,pt.design.G,'.-b')
    xlabel('r/R'), ylabel('G')
    plotaxes,

%%
% Create graphical and text reports
Make_Reports(pt)


%%
% Determine propeller geometry
pt.geometry = Geometry(pt)


%%
close all,
pt.input.Plot_flag = 0;

% Analyze off-design states
Js_all    = [0.4:0.05:0.75];     % advance coefficient
LAMBDAall = pi./Js_all;           % tip-speed ratio
pt.states = Analyze(pt,LAMBDAall)


figure, hold on,
    plot(pt.states.Js,   pt.states.KT  ,'.-b')
    plot(pt.states.Js,10*pt.states.KQ  ,'.-r')
    plot(pt.states.Js,   pt.states.EFFY,'.-g')
    
    plot(pt.input.Js,   pt.design.KT  ,'x')
    plot(pt.input.Js,10*pt.design.KQ  ,'x')
    plot(pt.input.Js,   pt.design.EFFY,'x')
    
    xlabel('Js'), ylabel('KT, 10*KQ, EFFY')
    axis([0.4 0.8 0 0.9])

    


