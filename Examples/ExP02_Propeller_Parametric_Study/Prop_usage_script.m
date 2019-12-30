% -------------------------------------------------------------------------
% Propeller parametric study example.
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
addpath ./SourceCode ../SourceCode ../../SourceCode
clr,

Prop_input,

D  = pt.i.D;              % propeller diameter [m]
Vs = pt.i.Vs;             % ship velocity [m/s]


JsALL = [0.1 0.2:0.2:2.0];  
NJS   = length(JsALL);

for j = 1:NJS
    disp(['>>>>>>>>------------ j = ',num2str(j),' -----------------<<<<<<<<<<'])
    
    Js        = JsALL(j);    	% advance ratio == Vs/nD
    
    pt.i.N   = 60*Vs/(Js*D);    % propeller speed [RPM]
    pt.i.Js  = Js;              % advance coefficient, Js = Vs/nD = pi/L
    pt.i.L   = pi/Js;           % tip speed ratio, L = omega*R/V
    
    
    % Perform design optimization:
    
    pt.i.Viscous_flag = 0;       
       d(j) = EppsOptimizer(pt.i);         
         
    
    pt.i.Viscous_flag = 1;
       v(j) = EppsOptimizer(pt.i);
end

% save workspace1

    disp(['>>>>>>>>------------ workspace1 saved  -----------------<<<<<<<<<<'])
    
% -------------------------------------------------------------------------
% Plot results
% -------------------------------------------------------------------------
% clr,
green = [0 0.8 0];

EDISC = 0.8970;

fig;
	HA = plot([0 2.5],EDISC*[1 1],'--','LineWidth',2,'color','k');
	axis([0 2 0 1])


% ---------------------------- 
% load workspace1

clear EFFYd EFFYv

for i = 1:NJS
    EFFYd(i) = d(i).EFFY;    
    EFFYv(i) = v(i).EFFY;
end

   HI = plot(JsALL,EFFYd,'b.-','LineWidth',2,'MarkerSize',20);
   HV = plot(JsALL,EFFYv,'r.-','LineWidth',2,'MarkerSize',20); 
% ---------------------------- 


SmallFontSize = 22;
LargeFontSize = 24;

set(gca,'XTick',[0:0.2:2.4],'XtickLabel',{'0','','0.4','','0.8','','1.2','','1.6','','2.0','','2.4'})
set(gca,'YTick',[0:0.1:1.0],'YtickLabel',{'0','','0.2','','0.4','','0.6','','0.8','','1.0'})
   
HL = legend([HA,HI,HV],'Inviscid 1D (actuator disk)','Inviscid with swirl','Viscous with swirl');

set(HL,'location','SouthEast','FontSize',SmallFontSize-4,'FontName','Times')

xlabel('Js'   ,'FontSize',LargeFontSize,'FontName','Times'), 
ylabel('EFFY' ,'FontSize',LargeFontSize,'FontName','Times')
set(gca,'FontSize',SmallFontSize,'FontName','Times')
axis([0 2.0 0 1.0]), box on, grid on,    
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%%









% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%%
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Off design analysis of inviscid cases
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
clr,
addpath ../SourceCode

Prop_input,
clc,

i = pt.i;

i.ITER         = 100;

i.Viscous_flag = 0;
i.Chord_flag   = 1;
i.XCD          = 0;
i.XCLmax       = 0.2;
i.dCLdALPHA    = 2*pi;
i

D  = pt.i.D;              % propeller diameter [m]
Vs = pt.i.Vs;             % ship velocity [m/s]


deltaJ = 0.01;     % resolution of performance curves
Lmax   = pi/2;     % maximum L for off-design analysis

JsALL = [0.1 0.2:0.2:2.0];  % ConeyProp110617
NJS   = length(JsALL);

for j = 1:NJS
    j 
    
    Js        = JsALL(j);    	% advance ratio == Vs/nD
    
    i.N   = 60*Vs/(Js*D);    % propeller speed [RPM]
    i.Js  = Js;              % advance coefficient, Js = Vs/nD = pi/L
    i.L   = pi/Js;           % tip speed ratio, L = omega*R/V
    
    
    % Perform design optimization
    pt(j).filename = pt(1).filename;
    pt(j).date     = pt(1).date;
    pt(j).i        = i;
    
    pt(j).d = EppsOptimizer(i);       
    
    % Off-design analysis
    pt(j).s = AnalyzeAuto(pt(j),deltaJ,Lmax);
end



save workspace2

%%
% -------------------------------------------------------------------------
% Plot results:
% -------------------------------------------------------------------------
clr,
fig; 
load colors
green  = [0 0.8 0];
color2 = [0.75 0 0.75];        
color3 = [0 0 1];
osize  = 8;


% ----------------------------    
% Actuator disk limit
    EDISC = 0.8970;
	plot([0 2.5],EDISC*[1 1],'k--','LineWidth',2)
% ----------------------------    

% ----------------------------     
% Required KT
JSplot = 0:0.01:2.5;
CT0 = 0.512;
CT1 = 0.512 - 0;
CT2 = 0.512 + 0.3;

KT0 = pi/8 * CT0 * JSplot.^2;
KT1 = pi/8 * CT1 * JSplot.^2;
KT2 = pi/8 * CT2 * JSplot.^2;

JSPLOTall = [JSplot(1:end),JSplot(end:-1:1)];
KTPLOTall = [   KT1(1:end),   KT2(end:-1:1)];
HF = fill(JSPLOTall,KTPLOTall,[0.8,0.8,0.8]);

% HJ = plot(JSplot,KT1,'m-.','LineWidth',2);
% HJ = plot(JSplot,KT2,'m-.','LineWidth',2);
HJ = plot(JSplot,KT0,'-','LineWidth',2,'color',green);
% ----------------------------    
    

% ----------------------------    
load workspace2 

for j = 1:NJS
    plot(pt(j).i.Js,   pt(j).d.KT   ,'^' ,'LineWidth',2,'MarkerSize',4,'MarkerFaceColor',green,'color',green);
    plot(pt(j).i.Js,   pt(j).d.EFFY ,'rv','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r');
    
       Jall(j) = pt(j).d.Js;
    EFFYall(j) = pt(j).d.EFFY;
end

plot(Jall,EFFYall ,'r-','LineWidth',2);
    

    HKQ = plot([-1 -1],[-1 -2],'^-' ,'LineWidth',2,'MarkerSize',4,'MarkerFaceColor',green,'color',green);
               
    
% ------------------------------------------------------------ Color matrix
CLRS = [    0.75    0       0.75;   ... % (4) Purple
            0.25    0.25    0.75;   ... % (9) Navy blue
            0.25    0.5     0.25;   ... % (8) Forest Green
            1       0.5     0;      ... % (5) Orange
            0.75    0.75    0;      ... % (11) Burnt Yellow
            0.25    0.5     0.75];      % (10) Steel blue 
           
% ----------------------------    OFF DESIGN STATES
% for j = 2:2:NJS-2
count = 0;
for j = NJS-3: -2: 2
    count = count + 1;
    
    color = CLRS(  mod(count-1,size(CLRS,1))+1  ,:);

    if ~isempty(pt(j).s)
    
        indices = find(pt(j).s.KT > 0 & pt(j).s.KQ > 0);
        
        HT = plot(pt(j).s.Js(indices),   pt(j).s.KT(indices)  ,'-','LineWidth',2,'MarkerSize',20,'Color',color);
        HE = plot(pt(j).s.Js(indices),   pt(j).s.EFFY(indices),'-','LineWidth',2,'MarkerSize',20,'Color',color);
        
    end

    indices = find(pt(j).s.KT ~= 0 & pt(j).s.KQ ~= 0);
    
    J1  = interp1(pt(j).s.CT(indices),pt(j).s.Js(indices),CT1);
    J2  = interp1(pt(j).s.CT(indices),pt(j).s.Js(indices),CT2);
    
    KT1 = interp1(pt(j).s.Js,pt(j).s.KT,J1);
    KT2 = interp1(pt(j).s.Js,pt(j).s.KT,J2);

    EFFY1 = interp1(pt(j).s.Js,pt(j).s.EFFY,J1);
    EFFY2 = interp1(pt(j).s.Js,pt(j).s.EFFY,J2);
    
    plot(J1*[1 1],[KT1 EFFY1],'k--','Linewidth',1)
    plot(J2*[1 1],[KT2 EFFY2],'k--','Linewidth',1)
    
    
    HT2 = plot(pt(j).i.Js,   pt(j).d.KT   ,'.','LineWidth',2,'MarkerSize',20,'Color',color);
    HE2 = plot(pt(j).i.Js,   pt(j).d.EFFY ,'.','LineWidth',2,'MarkerSize',20,'Color',color);
    
    HT2 = plot(J1,   KT1   ,'.','LineWidth',2,'MarkerSize',20,'Color',color);
    HE2 = plot(J1,   EFFY1 ,'.','LineWidth',2,'MarkerSize',20,'Color',color);
    
    HT2 = plot(J2,   KT2   ,'.','LineWidth',2,'MarkerSize',20,'Color',color);
    HE2 = plot(J2,   EFFY2 ,'.','LineWidth',2,'MarkerSize',20,'Color',color);
end
% ----------------------------


axis([0 2.0 0 1.0]),


set(gca,'XTick',[0:0.2:2.0],'XtickLabel',{'0','','0.4','','0.8','','1.2','','1.6','','2.0'})
set(gca,'YTick',[0:0.1:1.0],'YtickLabel',{'0','','0.2','','0.4','','0.6','','0.8','','1.0'})

    
% HL = legend([H02,HKQ,HF],'design EFFY','design KT','KT envelope');


SmallFontSize = 22;
LargeFontSize = 24;


% set(HL,'location','SouthEast','FontSize',SmallFontSize,'FontName','Times')

xlabel('Js'       ,'FontSize',LargeFontSize,'FontName','Times'), 
ylabel('KT, EFFY' ,'FontSize',LargeFontSize,'FontName','Times')
set(gca           ,'FontSize',SmallFontSize,'FontName','Times')
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
























