% -------------------------------------------------------------------------
% Example case in (Epps et al., FAST 2011)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
addpath ./SourceCode ../SourceCode ../../SourceCode
clr,

% Load inputs:
DDGprop_input


% Perform design optimization
pt.design = EppsOptimizer(pt.input)


%%
pt.geometry = Geometry(pt)


%%
pt.states = AnalyzeAuto(pt)

s1 = pt.states;

pt.input = rmfield(pt.input,'dCLdALPHA');

s2 = AnalyzeAuto(pt);

fig;
    plot(s1.Js,s1.KT   ,'b')
    plot(s1.Js,s1.KQ*10,'b')
    plot(s1.Js,s1.EFFY ,'b')
    
    
    plot(s2.Js,s2.KT   ,'r')
    plot(s2.Js,s2.KQ*10,'r')
    plot(s2.Js,s2.EFFY ,'r')
    
    
    xlabel('Js'), ylabel('KT, 10*KQ, EFFY')
    

    

