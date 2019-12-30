% -------------------------------------------------------------------------
% Example performance analysis from given propeller geometry.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
addpath ./SourceCode ../SourceCode ../../SourceCode
clr,

Prop4119_ginput,


% -------------------------------------------------------------------------
% Analyze the propeller at its design point:
s1 = AnalyzeGeometry(pt.ginput, pi/pt.ginput.Js)
% -------------------------------------------------------------------------


% Circulation distribution
fig;
    plot(s1.RC,s1.G,'b.-','MarkerSize',20)
    axis([0.2 1 0 0.04])
    xlabel('r/R'), ylabel('G')
    plotaxes,
    
% (Jessup, 1989) - "EXPERIMENTAL DATA - Visc. Cor."
    RG4119EV = [ 
      2.257142961E-01  1.625766873E-02
      2.500000000E-01  2.009202480E-02
      2.742857337E-01  2.131901741E-02
      3.000000119E-01  2.331288338E-02
      3.500000238E-01  2.806748390E-02
      4.014285803E-01  3.036809921E-02
      4.499999881E-01  3.205521345E-02
      5.014286041E-01  3.276073694E-02
      5.514286160E-01  3.374233246E-02
      6.014285684E-01  3.358895779E-02
      6.500000358E-01  3.358895779E-02
      7.014285326E-01  3.113496780E-02
      7.514285445E-01  3.052147388E-02
      8.000000119E-01  2.898772955E-02
      8.757143021E-01  2.453987837E-02
      8.999999762E-01  2.546012402E-02
      9.228571653E-01  2.009202480E-02 
      9.499999881E-01  5.828220844E-03
      1.000000000E+00  0.000000000E+00];    
    
    plot(RG4119EV(:,1),RG4119EV(:,2),'x-k','MarkerSize',10,'LineWidth',1);     


% -------------------------------------------------------------------------
%%




% -------------------------------------------------------------------------
% Analyze states:
Js_all        = [0.25:0.05:1.1, 1.12:0.02:1.16];     % advance coefficient
LAMBDAall     = pi./Js_all;                       % tip-speed ratio
    pt.states = AnalyzeGeometry(pt.ginput,LAMBDAall);
s = pt.states
% -------------------------------------------------------------------------

%%
fig;
    plot(s.Js,s.EFFY ,'.-' ,'LineWidth',2,'MarkerSize',20,'Color',[0 0.8 0])
    plot(s.Js,s.KT   ,'.-b','LineWidth',2,'MarkerSize',20)
    plot(s.Js,10*s.KQ,'.-r','LineWidth',2,'MarkerSize',20)    
   
    % Design point
    plot(pt.ginput.Js*[1 1],[0 2],'k--','LineWidth',1);
    
    xlabel('Js','Fontsize',18), 
    ylabel('KT, 10*KQ, EFFY','Fontsize',18)
    axis([0.25 1.2 0 1.0])
    set(gca,'Fontsize',16)
    

% ------- Experimental data for propeller 4119
JKTKQ = [    
    0.5000    0.2800    0.4600    0.4844
    0.7000    0.2050    0.3600    0.6344
    0.8330    0.1500    0.2800    0.7102
    0.9000    0.1250    0.2400    0.7460
    1.0000    0.0750    0.1800    0.6631
    1.1000    0.0300    0.1100    0.4775];

 J4119   = JKTKQ(:,1);
KT4119   = JKTKQ(:,2);
KQ4119   = JKTKQ(:,3);
EFFY4119 = JKTKQ(:,4);  % == (KT4119./(KQ4119/10)).*(J4119/(2*pi));

    % Propeller 4119 data
    plot(J4119,KT4119  ,'x-k','MarkerSize',10,'LineWidth',2)%,'MarkerFaceColor','k')
    plot(J4119,KQ4119  ,'x-k','MarkerSize',10,'LineWidth',2)%,'MarkerFaceColor','k')
    plot(J4119,EFFY4119,'x-k','MarkerSize',10,'LineWidth',2)%,'MarkerFaceColor','k')
  


