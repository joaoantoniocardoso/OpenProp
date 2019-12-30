% ----------------------------------------------------------------------- %
%                                                                         %
%                              0111000                                    %
%                           100 1 100 001                                 %
%                         10    1  1  1 00                                %
%                        01  1  1  1      0                               %
%                       0 1  1  1   1  1 1 0                              %
%                       0   1   1   1  1  1 0                             %
%                       0 1     1   1  1  1 0                             %
%                       0 1  1  1   1  0  1 0                             %
%                       0 1  1  1   0  1    0                             %
%                       01 1        1  1 1 0                              %
%                        0    0  1  0 1   0                               %
%                         0         1    0                                %
%                    10010 0 1101111110111                                %
%                  10 1 1  1111111111 11 11                               %
%                 0 1 1 1 11111111101011010111                            %
%                01 11    11111111 1  1    1 110                          %
%               011    1 1 111111110011  1 1 1 110                        %
%               0   11 1 1 1 111      0  1 1 1   10                       %
%               0 1   11  1  0         1 1 1 1 1 1 0                      %
%               1  11 1 1   11          0  1 1 1 1 11                     %
%                0     1 1  0           011  1 1 1 10                     %
%                10 1   1  0             0  1 1 1  11                     %
%                 10     01               01      10                      %
%                   10001                   001 100                       %
%                                             111                         %
%                                                                         %
%             ____                   _____                                %
%            / __ \                 |  __ \                               %
%           | |  | |_ __   ___ _ __ | |__) | __ ___  _ __                 %
%           | |  | | '_ \ / _ \ '_ \|  ___/ '__/ _ \| '_ \                %
%           | |__| | |_) |  __/ | | | |   | | | (_) | |_) |               %
%            \____/| .__/ \___|_| |_|_|   |_|  \___/| .__/                %
%                  | |                              | |                   %
%                  |_|                              |_|                   %
%                                                                         %
%             An integrated rotor design and analysis tool.               %
%                                                                         %
%                                                                         %
% Copyright (C) 2011, Brenden Epps.                                       %
%                                                                         %
% This program is free software; you can redistribute it and/or modify it %
% under the terms of the GNU General Public License as published by the   %
% Free Software Foundation.                                               %
%                                                                         %
% This program is distributed in the hope that it will be useful, but     %
% WITHOUT ANY WARRANTY; without even the implied warranty of              %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    %
% See the GNU General Public License for more details.                    %
%                                                                         %
% ----------------------------------------------------------------------- %


% =========================================================================
% OpenProp_v3.3.4
% Last modified: 10/21/2013 Brenden Epps
% =========================================================================
%
% This code runs the parametric-design GUI.
%
%--------------------------------------------------------------------------

% =========================================================================
% ========================== Parametric Analysis ==========================

function OpenPropParam

clear all;

addpath ../SourceCode


%%

% ------------------------- Setup global variables ------------------------

global Param_Main;    % Main GUI figure

global Select;

global OpenPropDirectory SpecificationsValues FlagValues Filename filename RangeValues... DuctValues FoilValues
       XR_in XCoD_in XCD_in VAI_in VTI_in ri_in;
%        ... % Xt0oD_in skew0_in rake0_in...
%        Meanlinecell Thicknesscell ; % CavValues Values

global N_R0;                % Number of input radii

global Col_Label;

global XCoD_values XCLmax_values XCD_values;

% --- Set GUI element variables as global ---


% =========================================================================

%%

% =========================================================================
% ========================== Initiate Main Figure =========================


% -------------------- Declare Default Input Variables --------------------

% Variables for geometry input table
XR_def          = [.2 .3 .4 .5 .6 .7 .8 .9 .95 1];          % Radius ratio (r/R) at control point
N_R0            = length(XR_def);                           % Number of input radii

XCD_def         = ones(1,N_R0).*0.008;                      % Coefficient of drag at XR

XCoD_def        = [0.1600 0.1812 0.2024 0.2196 0.2305 0.2311 0.2173 0.1807 0.1388 0.0010];     % chord to diameter ratio at XR


XCLmax_def      = 0.5 + (0.2-0.5)*(XR_def-XR_def(1))/(XR_def(end)-XR_def(1));  % CLmax distribution

% --- Set defaults for change callbacks ---
  XCoD_values   =   XCoD_def;
XCLmax_values   = XCLmax_def;
   XCD_values   =    XCD_def;
% -----------------------------------------


% General variables
Z_def           = 3;                                        % Number of blades
N_def           = 200;                                      % Propeller speed [RPM]
D_def           = 2;                                        % Propeller diameter [m]
T_def           = 40000;                                    % Required thrust [N]
Vs_def          = 5;                                        % Ship speed [m/s]
Dhub_def        = 0.4;                                      % Hub diameter [m]


% Advanced variables
rho_def         = 1000;                                     % Water density
Mp_def          = 20;                                       % Number of vortex panels over the radius
Np_def          = 20;                                       % Number of points over the chord


n_def           = N_def/60;                                 % ** propeller speed [rev/s] = Vs/(Js*D) = N/60
lambda_def      = n_def*2*pi*(D_def/2)/Vs_def;
Js_def          = Vs_def/(n_def*D_def);                     % ** Js = Vs/(n*D) ,  advance coefficient
KT_def          = T_def/(rho_def*n_def^2*D_def^4);          % ** KT = THRUST/(rho*n^2*D^4)

CT_def          = T_def/(0.5*rho_def*Vs_def^2*pi*(D_def/2)^2);


% Ducted Propeller variables
TAU_def         = 1;                                        % Thrust ratio
CDd_def         = 0.008;                                    % Duct section drag coefficient

% Parametric analysis defaults
BladeMin        = 3;
BladeMax        = 6;

SpeedMin        = 100;
SpeedMax        = 200;
SpeedInc        = 50;

DiameterMin     = 2;
DiameterMax     = 4;
DiameterInc     = .5;




filename        = 'DefaultPropeller';                           % Filename prefix

% =========================================================================


% -------------------------- GUI Switching Check --------------------------
% 

if exist('OpenPropTempFile0307122010.mat','file')
    
    load('OpenPropTempFile0307122010.mat');
    
    Z_def       = Z;
    N_def       = N;
    D_def       = D;
    T_def       = THRUST;
    Vs_def      = Vs;
    Dhub_def    = Dhub;
    
    rho_def     = rho;
    Mp_def      = Mp;
    Np_def      = Np;
    
    XR_def      = XR;
    XCoD_def    = XCoD;
    XCD_def     = XCD;
%     'VAI',
%     'VTI',
%     'ri'
    
    
    n_def           = N_def/60;                                 % ** propeller speed [rev/s] = Vs/(Js*D) = N/60
    lambda_def      = n_def*2*pi*(D_def/2)/Vs_def;
    Js_def          = Vs_def/(n_def*D_def);                     % ** Js = Vs/(n*D) ,  advance coefficient
    KT_def          = T_def/(rho_def*n_def^2*D_def^4);          % ** KT = THRUST/(rho*n^2*D^4)
    
    CT_def          = T_def/(0.5*rho_def*Vs_def^2*pi*(D_def/2)^2);
    
    
    % Ducted Propeller variables
    TAU_def         = 1;                                        % Thrust ratio
    CDd_def         = 0.008;                                    % Duct section drag coefficient
    
    delete('OpenPropTempFile0307122010.mat');
    
end

% =========================================================================


% -------------------------- GUI Layout Constants -------------------------
% 
% This section presents the constant values for the dimensions of those
% elements used in the GUI, as well as the construction formulas for the
% different panels, including margins. The formulas are presented in a
% linear reading order based on the actual order in the GUI.

titlefontsize   = 25;
subttlfontsize  = 20;
panelfontsize   = 11;
buttonfontsize  = 12;

titleht         = 3;
editbox         = 10;
textbox         = 20;
cavtextbox      = 22;
filenametext    = 15;
pushbox         = 10;
runbox          = 20;

selectboxht     = 2;
selectbox       = 25;

textht          = 1.5;
editboxht       = 2;
pushht          = 3;

Specificationsboxht	= 1 + editboxht * 11 + 2;
Specificationsbox 	= 2 + textbox + editbox + 2;

BladeDesignboxht	= 1 + editboxht * 11 + 2;
BladeDesignbox      = 2 + editbox * 3 + 2;

Inflowboxht    	= 1 + editboxht * 11 + 2;
Inflowbox      	= 2 + editbox * 3 + 2;


Flagboxht       = BladeDesignboxht; % 2 + textht * 8 + 2;
Flagbox         = 1 + textbox + 1;

Toolboxht       = ( 1 + pushht + 1 + editboxht + 2 + 1);
Toolbox         = ( Specificationsbox + BladeDesignbox + Inflowbox...
                    + Flagbox ) / 2; % - Ductbox; % - Valuesbox;

Rangeboxht      = Toolboxht;
Rangebox        = Toolbox;

filenamebox     = Toolbox - (2 + filenametext + 2);

buttonspace     = (Toolbox-pushbox*2-runbox)/4;

Windowht        = 1 + Toolboxht + Specificationsboxht + 1 + titleht + 1; % + Ductboxht
Window          = 1 + Specificationsbox + BladeDesignbox + Inflowbox + Flagbox + 1;

GUIselectionboxht   = 1 + selectboxht + 1;
GUIselectionbox     = 1 + selectbox + 1;

% =========================================================================

%%

% ------------------------- Create figure for GUI -------------------------

close all;

Param_Main    = figure('units','characters','position',[5 55-Windowht Window Windowht],...
                     'numbertitle','off','name','OpenProp','menubar','none',...'toolbar','figure',...
                     'resize','off','color',[0.702 0.702 0.702]);

if strcmp(computer,'GLNX32') || strcmp(computer,'GLNXA64')
    
    set(Param_Main,'resize','on');
    
end

OpenPropDirectory = 'OpenProp_v3.3.4';
OpenPropVersion   = 'OpenProp v3.3.4';


Title       = uicontrol(Param_Main,'style','text','fontsize',titlefontsize,...
                        'fontweight','bold','units','characters','position',...
                        [GUIselectionbox-12 Windowht-1-titleht Window-GUIselectionbox 3],'string',OpenPropVersion);


% -------------------------------
% --------- Setup panels --------

Specifications	= uipanel('parent',Param_Main,'title','Specifications','fontsize',panelfontsize,...
                      'fontweight','bold','units','characters','position',...
                      [1 1+Toolboxht Specificationsbox Specificationsboxht],'clipping','on');

BladeDesign     = uipanel('parent',Param_Main,'title','Blade Design Values','fontsize',...
                      panelfontsize,'fontweight','bold','units','characters',...
                      'position',[1+Specificationsbox 1+Toolboxht BladeDesignbox BladeDesignboxht],...
                      'clipping','on');

Inflow          = uipanel('parent',Param_Main,'title','Inflow Profile Values','fontsize',...
                      panelfontsize,'fontweight','bold','units','characters',...
                      'position',[1+BladeDesignbox+Specificationsbox...
                      1+Toolboxht Inflowbox Inflowboxht],...
                      'clipping','on');

Flags           = uibuttongroup('parent',Param_Main,'title','Options','fontsize',...
                      panelfontsize,'fontweight','bold','units','characters',...
                      'position',[1+Specificationsbox+BladeDesignbox+Inflowbox...
                      1+Toolboxht Flagbox Flagboxht],'clipping','on','SelectionChangeFcn',@checkTurbine);

Tools           = uipanel('parent',Param_Main,'title','Tools','fontsize',...
                      panelfontsize,'fontweight','bold','units','characters',...
                      'position',[1+Rangebox 1 Toolbox Toolboxht],...
                      'clipping','on');

Range           = uipanel('parent',Param_Main,'title','Range','fontsize',...
                      panelfontsize,'fontweight','bold','units','characters',...
                      'position',[1 1 Toolbox Toolboxht],...
                      'clipping','on');

% -------------------------------
% -------------------------------

% =========================================================================

%%
% ----------------------- Package Selection Elements ----------------------

% % Select(1)     = uicontrol(Param_Main,'units','characters','style',...
% %                             'pushbutton','string','PS','position',...
% %                             [1+1 1+Rangeboxht+Specificationsboxht+2 selectbox selectboxht],...
% %                             'horizontalalignment','left','value',1,...
% %                             'tooltipstring','Parametric Study');
% % 
% % Select(2)     = uicontrol(Param_Main,'units','characters','style',...
% %                             'pushbutton','string','SP','position',...
% %                             [1+1+selectbox 1+Rangeboxht+Specificationsboxht+2 selectbox selectboxht],...
% %                             'horizontalalignment','left','callback','OpenPropSingle',...
% %                             'tooltipstring','Single Propeller Design');
% % 
% % Select(3)     = uicontrol(Param_Main,'units','characters','style',...
% %                             'pushbutton','string','A','position',...
% %                             [1+1+2*selectbox 1+Rangeboxht+Specificationsboxht+2 selectbox selectboxht],...
% %                             'horizontalalignment','left','callback','OpenPropAnalyze',...
% %                             'tooltipstring','Analyze Propeller');

% Selectcell      = {'Single Design','Parametric Study','Off-design Analysis'};
Selectcell      = {'Single Design','Parametric Study'};

Select	= uicontrol(Param_Main,'units','characters','style','popupmenu',...
                            'position',[1+1 1+Rangeboxht+Specificationsboxht+2 selectbox selectboxht],...
                            'backgroundcolor','w','string',Selectcell,'value',2,'callback',@Selectfn);

% =========================================================================

%%

% --------------------- Specifications Panel Elements ---------------------

SpecificationsStrings = {... 
                           'Required thrust (N):'...  
                           'Ship speed (m/s):'...
                           'Hub diameter (m):'...
                           'Fluid density (kg/m^3):' ...
                           '# radial panels:'...
                           '# chordwise panels:'};
                       
                           

SpecificationsValues_def      = [T_def Vs_def Dhub_def rho_def Mp_def Np_def];

for index = 1 : length(SpecificationsValues_def)
    
    SpecificationsText(index)     = uicontrol(Specifications,'units','characters','style',...
                                        'text','string',SpecificationsStrings(index),...
                                        'position',[2 Specificationsboxht-4-editboxht*(index-1)...
                                        textbox textht],'horizontalalignment','left');
    
    SpecificationsValues(index)   = uicontrol(Specifications,'units','characters','style',...
                                        'edit','string',num2str(SpecificationsValues_def(index)),...
                                        'position',[2+textbox Specificationsboxht-4-editboxht*(index-1)...
                                        editbox editboxht],'backgroundcolor',[1 1 1]);
end

% =========================================================================

%%

% ----------------------- Blade Design Panel Elements ---------------------

ColName     = {'r/R' 'c/D' 'Cd'};

for index = 1 : length(ColName)
    Col_Label(index) = uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
        'FontWeight','bold','position',[2+editbox*(index-1) BladeDesignboxht-4 editbox editboxht],...
        'string',ColName(index),'enable','inactive');
end

for index = 1 : N_R0
    
    XR_in(index)	= uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
                               'backgroundcolor','w','string',num2str(XR_def(index)),'position',...
                               [2 BladeDesignboxht-4-editboxht*index editbox editboxht]);
    
    XCoD_in(index)  = uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
                                'backgroundcolor','w','string',num2str(XCoD_def(index)),'position',...
                                [2+editbox BladeDesignboxht-4-editboxht*index editbox editboxht]);
    
    XCD_in(index)   = uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
                                'backgroundcolor','w','string',num2str(XCD_def(index)),'position',...
                                [2+editbox*2 BladeDesignboxht-4-editboxht*index editbox editboxht]);
     
end


% =========================================================================

%%

% ---------------------- Inflow Profile Panel Elements --------------------

ColName2     = {'r' 'Va/Vs' 'Vt/Vs'};

for index = 1 : length(ColName2)
    Col_Label2(index) = uicontrol(Inflow,'style','edit','units','characters','FontSize',10,...
        'FontWeight','bold','position',[2+editbox*(index-1) Inflowboxht-4 editbox editboxht],...
        'string',ColName2(index),'enable','inactive');
end

for index = 1 : N_R0
    
    ri_in(index)       = uicontrol(Inflow,'style','edit','units','characters','FontSize',10,...
                               'backgroundcolor','w','string','','position',...
                               [2 BladeDesignboxht-4-editboxht*index editbox editboxht]);
    
    VAI_in(index)   = uicontrol(Inflow,'style','edit','units','characters','FontSize',10,...
                                'backgroundcolor','w','position',...
                                [2+editbox BladeDesignboxht-4-editboxht*index editbox editboxht]);
    
    VTI_in(index)   = uicontrol(Inflow,'style','edit','units','characters','FontSize',10,...
                                'backgroundcolor','w','position',...
                                [2+editbox*2 BladeDesignboxht-4-editboxht*index editbox editboxht]);
    
end

% =========================================================================

%%

% -------------------------- Flag Panel Elements --------------------------

FlagValues(1)     = uicontrol(Flags,'units','characters','style','radiobutton',...
                              'string','Propeller','value',1,'position',...
                           [1 Flagboxht-4 textbox textht]);

% FlagValues(2)     = uicontrol(Flags,'units','characters','style','radiobutton',...
%                               'string','Turbine','position',...
%                               [1 Flagboxht-4-textht textbox textht]);

FlagValues(3)     = uicontrol(Flags,'units','characters','style','checkbox',...
                              'string','Hub','value',1,'position',...
                              [1 Flagboxht-4-textht*2 textbox textht]);

FlagValues(4)     = 0;

FlagValues(5)     = uicontrol(Flags,'units','characters','style','checkbox',...
                              'string','Chord optimization','position',...
                              [1 Flagboxht-4-textht*3 textbox textht],...
                              'callback',@changeChord);

FlagValues(6)     = uicontrol(Flags,'units','characters','style','checkbox',...
                              'string','Viscous forces','value',1,'position',...
                              [1 Flagboxht-4-textht*4 textbox textht],...
                              'callback',@changeViscous);

% FlagValues(7)     = uicontrol(Flags,'units','characters','style','checkbox',...
%                               'string','Optimization plots','position',...
%                               [1 Flagboxht-4-textht*5 textbox textht]);
% 
% FlagValues(8)     = uicontrol(Flags,'units','characters','style','checkbox',...
%                               'string','2D/3D plots','value',1,'position',...
%                               [1 Flagboxht-4-textht*6 textbox textht]);

% =========================================================================

%%

% -------------------------- Tools Panel Elements -------------------------

ToolText     = uicontrol(Tools,'units','characters','style','text',...
                            'string','Filename prefix:','horizontalalignment',...
                            'left','position',[2 Toolboxht-4 filenametext textht]);

Filename     = uicontrol(Tools,'units','characters','style','edit',...
                            'backgroundcolor','w','string',filename,'position',...
                            [2+filenametext Toolboxht-4 filenamebox editboxht],...
                            'callback',@changeDir);

LoadButton      = uicontrol(Tools,'units','characters','style','pushbutton',...
                            'string','Load','fontsize',buttonfontsize,...
                            'position',[buttonspace 1 pushbox pushht],...
                            'callback',@loaddata);

SaveButton      = uicontrol(Tools,'units','characters','style','pushbutton',...
                            'string','Save','fontsize',buttonfontsize,...
                            'position',[buttonspace*2+pushbox 1 pushbox pushht],...
                            'callback',@savedata);

RunButton       = uicontrol(Tools,'units','characters','style','pushbutton',...
                            'string','Run OpenProp','fontsize',buttonfontsize,...
                            'position',[buttonspace*3+pushbox*2 1 runbox pushht],...
                            'callback',@execute);

% =========================================================================

%%

% -------------------------- Range Panel Elements -------------------------

Rangetext(1)	= uicontrol(Range,'units','characters','style',...
                           'text','string','Number of Blades','position',...
                           [2 Rangeboxht-5 cavtextbox textht],...
                           'horizontalalignment','left');

Rangetext(2)	= uicontrol(Range,'units','characters','style',...
                           'text','string','Rotation Speed (RPM)','position',...
                           [2 Rangeboxht-5-editboxht cavtextbox textht],...
                           'horizontalalignment','left');

Rangetext(3)	= uicontrol(Range,'units','characters','style',...
                           'text','string','Rotor Diameter (m)','position',...
                           [2 Rangeboxht-5-editboxht*2 cavtextbox textht],...
                           'horizontalalignment','left');

Rangetext(4)	= uicontrol(Range,'units','characters','style',...
                           'text','string','Min','position',...
                           [2+cavtextbox Rangeboxht-4,...
                           editbox editboxht]);

Rangetext(5)	= uicontrol(Range,'units','characters','style',...
                           'text','string','Increment','position',...
                           [2+cavtextbox+editbox Rangeboxht-4,...
                           editbox editboxht]);

Rangetext(6)	= uicontrol(Range,'units','characters','style',...
                           'text','string','Max','position',...
                           [2+cavtextbox+editbox*2 Rangeboxht-4,...
                           editbox editboxht]);

% ---

% BladeMin ... look after 'string'
RangeValues(1)	= uicontrol(Range,'units','characters','style',...
                           'edit','string',BladeMin,...
                           'position',[2+cavtextbox Rangeboxht-5,...
                           editbox editboxht],'backgroundcolor',[1 1 1],...
                           'enable','on','callback',@checkBlades);
                       
RangeValues(2)	= uicontrol(Range,'units','characters','style',...
                           'edit','string',SpeedMin,...
                           'position',[2+cavtextbox Rangeboxht-5-editboxht,...
                           editbox editboxht],'backgroundcolor',[1 1 1],...
                           'enable','on');

RangeValues(3)	= uicontrol(Range,'units','characters','style',...
                           'edit','string',DiameterMin,...
                           'position',[2+cavtextbox Rangeboxht-5-editboxht*2,...
                           editbox editboxht],'backgroundcolor',[1 1 1],...
                           'enable','on');
% BladeInc
RangeValues(4)	= uicontrol(Range,'units','characters','style',...
                           'edit','string',1,...
                           'position',[2+cavtextbox+editbox Rangeboxht-5,...
                           editbox editboxht],'backgroundcolor',[1 1 1],...
                           'enable','on','callback',@checkBlades);

RangeValues(5)	= uicontrol(Range,'units','characters','style',...
                           'edit','string',SpeedInc,...
                           'position',[2+cavtextbox+editbox Rangeboxht-5-editboxht,...
                           editbox editboxht],'backgroundcolor',[1 1 1],...
                           'enable','on');

RangeValues(6)	= uicontrol(Range,'units','characters','style',...
                           'edit','string',DiameterInc,...
                           'position',[2+cavtextbox+editbox Rangeboxht-5-editboxht*2,...
                           editbox editboxht],'backgroundcolor',[1 1 1],...
                           'enable','on');

% BladeMax                       
RangeValues(7)	= uicontrol(Range,'units','characters','style',...
                           'edit','string',BladeMax,...
                           'position',[2+cavtextbox+editbox*2 Rangeboxht-5,...
                           editbox editboxht],'backgroundcolor',[1 1 1],...
                           'enable','on','callback',@checkBlades);

RangeValues(8)	= uicontrol(Range,'units','characters','style',...
                           'edit','string',SpeedMax,...
                           'position',[2+cavtextbox+editbox*2 Rangeboxht-5-editboxht,...
                           editbox editboxht],'backgroundcolor',[1 1 1],...
                           'enable','on');

RangeValues(9)	= uicontrol(Range,'units','characters','style',...
                           'edit','string',DiameterMax,...
                           'position',[2+cavtextbox+editbox*2 Rangeboxht-5-editboxht*2,...
                           editbox editboxht],'backgroundcolor',[1 1 1],...
                           'enable','on');

% =========================================================================














%%

end



function execute(hObject,ED)


global OpenPropDirectory SpecificationsValues FlagValues Filename filename RangeValues...
       XR_in XCoD_in XCD_in VAI_in VTI_in ri_in ...
       XCoD_values XCLmax_values; % CavValues
   
   global paroutput;

%%


filename   	= get(Filename,'string');                       % Filename prefix

% ------------------------------------------------------------------------- 
% Figure out what the current working directory is, 
% and change directories to OpenPropDirectory/filename
rest = pwd;

while ~isempty(rest)
    [CurrentDirectory,rest] = strtok(rest,'/');
    
    if strcmp(CurrentDirectory,OpenPropDirectory)

        if isempty(rest)
            
            % you are in /OpenPropDirectory/
            mkdir(['./',filename])
               cd(['./',filename])
               addpath ../SourceCode
               
        elseif strcmp(rest(2:end),filename)
            % already in /OpenPropDirectory/filename/
            addpath ../SourceCode
            rest = [];
            
        elseif strcmp(rest(2:end),'SourceCode')
            
            mkdir(['../',filename])
               cd(['../',filename])
               addpath ../SourceCode
            rest = [];
            
        else
            % you are in /OpenPropDirectory/wrongfolder
            disp('ERROR2: Must start OpenProp from the root directory.')
            return
        end
    end
end
% -------------------------------------------------------------------------





% --------------------------- Design parameters ---------------------------


THRUST      = str2double(get(SpecificationsValues(1),'string')); 	% required thrust [N]
Vs          = str2double(get(SpecificationsValues(2),'string'));  % ship velocity [m/s]
Dhub        = str2double(get(SpecificationsValues(3),'string'));  % hub diameter [m]
rho         = str2double(get(SpecificationsValues(4),'string')); 	% water density [kg/m^3]
Mp          = str2double(get(SpecificationsValues(5),'string')); 	% number of vortex panels over the radius
Np          = str2double(get(SpecificationsValues(6),'string')); 	% Number of points over the chord [ ]

ITER        = 40;   % number of iterations in wake alignment
Rhv         = 0.5;	% hub vortex radius / hub radius

Rhub     = Dhub/2;                               % hub radius [m]

Zmin        = str2double(get(RangeValues(1),'string'));  % min number of blades
ZIncrement	= str2double(get(RangeValues(4),'string'));  % increment in number of blades
ZMax        = str2double(get(RangeValues(7),'string'));  % max number of blades

Nmin        = str2double(get(RangeValues(2),'string'));  % min propeller speed [RPM]
NIncrement	= str2double(get(RangeValues(5),'string'));  % increment in propeller speed [RPM]
NMax        = str2double(get(RangeValues(8),'string'));  % max propeller speed [RPM]

Dmin        = str2double(get(RangeValues(3),'string'));  % min propeller diameter [m]
DIncrement	= str2double(get(RangeValues(6),'string'));  % increment in propeller diameter [m]
DMax        = str2double(get(RangeValues(9),'string'));  % max propeller diameter [m]

Propeller_flag	= get(FlagValues(1),'value');               % 0 == turbine, 1 == propeller
Hub_flag  	= get(FlagValues(3),'value');                   % 0 == no hub, 1 == hub

% Duct_flag	= get(FlagValues(4),'value');                   % 0 == no duct, 1 == duct
Duct_flag   = FlagValues(4);

Chord_flag	= get(FlagValues(5),'value');                   % ** CHORD OPTIMIZATION FLAG **

Viscous_flag	= get(FlagValues(6),'value');               % 0 == viscous forces off (CD = 0), 1 == viscous forces on

% Plot_flag       = get(FlagValues(7),'value');               % 0 == do not display plots, 1 == display plots
% Make2Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 2D plot of the results, 1 == make plot
% Make3Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 3D plot of the results, 1 == make plot

Plot_flag       = 0;               % 0 == do not display plots, 1 == display plots

Make2Dplot_flag = 0;               % 0 == do not make a 2D plot of the results, 1 == make plot
Make3Dplot_flag = 0;               % 0 == do not make a 3D plot of the results, 1 == make plot


Camber_flag     = 0;                                        % ** VERIFY USE, ADD TO GUI **


XR        = str2double(get(XR_in,'string'));                % radius / propeller radius
XCoD      = str2double(get(XCoD_in,'string'));              % chord / diameter
XCD       = str2double(get(XCD_in,'string'));               % section drag coefficient

ri        = str2double(get(ri_in, 'string'));
VAI       = str2double(get(VAI_in,'string'));               % axial      inflow velocity / ship velocity
VTI       = str2double(get(VTI_in,'string'));               % tangential inflow velocity / ship velocity

ri  = ri (~isnan(ri));
VAI = VAI(~isnan(VAI));
VTI = VTI(~isnan(VTI));

%%

Z           = Zmin:ZIncrement:ZMax;
N           = Nmin:NIncrement:NMax;
D           = Dmin:DIncrement:DMax;

dVs         = 0.3/Vs;                                     	% axial inflow variation / Vs

%%

%% ------------------------------------------------- Unpack input variables
% '------ Performance inputs ------'
parinput.Z    = Z;           % [numberZ x 1], [ ]   number of blades
parinput.N    = N;           % [numberN x 1], [RPM] rotation rate
parinput.D    = D;           % [numberD x 1], [m]   propeller diameter

parinput.Vs      = Vs;          % [1 x 1]
parinput.THRUST  = THRUST;       % [1 x 1]
% '------ Geometry inputs ------'
parinput.Mp      = Mp;          % [1 x 1]
parinput.Rhub    = Rhub;        % [1 x 1]
parinput.XR      = XR;          % [length(XR) x 1]

% parinput.VAI     = VAI;         % [length(XR) x 1]
% parinput.VTI     = VTI;         % [length(XR) x 1]

parinput.XCD     = XCD;         % [length(XR) x 1]
parinput.XCoD    = XCoD;        % [length(XR) x 1]
% '------ Computational inputs ------'
parinput.ITER    = ITER;        % [1 x 1]
parinput.Propeller_flag = Propeller_flag;
parinput.Viscous_flag = Viscous_flag; % 0 == viscous forces off (CD = 0), 1 == viscous forces on
parinput.Hub_flag     = Hub_flag;     % 0 == no hub, 1 == hub
parinput.Duct_flag    = Duct_flag;    % 0 == no duct, 1 == duct
parinput.Chord_flag    = Chord_flag;    % 0 == no duct, 1 == duct
% parinput.HUF     = HUF;         % [1 x 1]
% parinput.TUF     = TUF;         % [1 x 1]
% parinput.SCF     = SCF;         % [1 x 1]
parinput.Rhv     = Rhv;         % [1 x 1]
% '------ Cavitation inputs ------'
parinput.rho     = rho;         % [1 x 1], [kg/m^3], density of seawater
parinput.dVs     = dVs;         % [1 x 1]

% parinput.H       = H;           % [1 x 1]

% '------ Duct inputs ------'
% parinput.TAU     = TAU;         % [1 x 1]
% parinput.Rduct_oR= Rduct_oR;    % [1 x 1]
% parinput.CDd     = CDd;         % [1 x 1]
% parinput.CDoCL   = CDoCL;



if ~isempty(ri)  , parinput.ri  = ri;  end
if ~isempty(VAI) , parinput.VAI = VAI; end        % [length(XR) x 1], [ ] input axial inflow velocity  at XR
if ~isempty(VTI) , parinput.VTI = VTI; end        % [length(XR) x 1], [ ] input swirl inflow velocity  


%%

paroutput        = LerbsParametric(parinput);

% 
% %%
% 
% Fig_Param = figure('units','characters','position',[15 15 130 40],...
%                    'name','Efficiency','numbertitle','off');
%     for i = 1:length(Z)
%         if length(Z)==4
%            subplot(2,2,i);
%         else
%            subplot(length(Z),1,i);
%         end
%         plot(paroutput.D,paroutput.EFFY(:,:,i));
%         str_suffix={' RPM'};
%         for j=1:length(N)
%             str_legend(j)=strcat(num2str(N(j)),str_suffix);
%         end
%         legend(str_legend,'location','southwest');    grid on;
%         xlabel('Propeller Diameter (m)');             ylabel('Efficiency');
%         title_prefix = {'Number of Blades: '};  
%         title(strcat(title_prefix,num2str(Z(i))))
%     end
% %     set(New_Parametric,'enable','on');    
% %     figure(Fig_Main);


%%
% ------------------------------------------------------------ Color matrix
CLR = [     1       0       0;      ... % (1) Red
            0       0.9     0;      ... % (2) Green
            0       0       1;      ... % (3) Blue
            0.75    0       0.75;   ... % (4) Purple
            1       0.5     0;      ... % (5) Orange
            0       1       1;      ... % (6) Cyan
            1       0       1;      ... % (7) Magenta
            0.75    0.5     0.25;   ... % (8) Brown
            0.25    0.25    0.75;   ... % (9) Navy blue
            0.25    0.5     0.75;   ... % (10) Steel blue
            0.75    0.75    0];         % (11) Burnt Yellow

        
Fig_Param = figure('units','characters','position',[15 15 130 40],...
                   'name','Efficiency','numbertitle','off');
               
    for iZ = 1:length(Z)
        
        str_legend = {0*N};
        
        if length(Z)==4
           subplot(2,2,iZ);
        else
           subplot(length(Z),1,iZ);
        end
        hold on, box on, grid on,
        ylim([0 1]),
        
        for iN = 1:length(N)
            tempEFFY = paroutput.EFFY(iZ,iN,:);
            
            plot(D,tempEFFY(:),'-','color',CLR( (mod(iN-1,11)+1) ,:));
            
            str_legend(iN)={[num2str(N(iN)),' RPM']};
        end
        
        legend(str_legend,'location','southwest');    
        
        xlabel('Propeller Diameter (m)');             
        ylabel('Efficiency');
        title(['Number of Blades: ',num2str(Z(iZ))])
    end
%     set(New_Parametric,'enable','on');    
%     figure(Fig_Main);

%%

% ========================= pt overwrite sequence =========================

pt.filename     = filename;
pt.date         = date;
pt.parinput     = parinput;
pt.paroutput	= paroutput;

pt.input.GUI_flag     = 1;

if exist([filename '.mat'],'file')
    
    temp            = pt;
    
    warning off
    load(filename)
    warning on
    
    pt.filename     = temp.filename;
    pt.date         = temp.date;
    pt.parinput     = temp.parinput;
    pt.paroutput	= temp.paroutput;
    
end

% =========================================================================

save(filename,'pt');


save([filename '_GUIps'],'Zmin','ZIncrement','ZMax','Nmin','NIncrement','NMax',...
    'Dmin','DIncrement','DMax','THRUST','Vs','Dhub','rho','Mp','Np',...
    'Propeller_flag','Hub_flag','Duct_flag',...
    'Chord_flag','Viscous_flag','Plot_flag','Make2Dplot_flag','Make3Dplot_flag',...
    'Camber_flag','filename','XR','XCoD','XCD','VAI','VTI','ri');



end



function savedata(hObject,ED)

%%

% =========================================================================
% ================================ savedata ===============================
% 
% This subfunction saves all values presented in the OpenProp GUI.
%

global OpenPropDirectory SpecificationsValues FlagValues Filename filename RangeValues...
       XR_in XCoD_in XCD_in VAI_in VTI_in ri_in ...
       XCoD_values XCLmax_values; % CavValues
   
   global paroutput;

%%


filename   	= get(Filename,'string');                       % Filename prefix


% ------------------------------------------------------------------------- 
% Figure out what the current working directory is, 
% and change directories to OpenPropDirectory/filename
rest = pwd;

while ~isempty(rest)
    [CurrentDirectory,rest] = strtok(rest,'/');
    
    if strcmp(CurrentDirectory,OpenPropDirectory)

        if strcmp(rest(2:end),filename)
            % already in /OpenPropDirectory/filename/
            rest = [];
            
            addpath ../SourceCode
            
        elseif strcmp(rest(2:end),'SourceCode')
            
            mkdir(['../',filename])
               cd(['../',filename])
            rest = [];
            
            addpath ../SourceCode
            
        elseif isempty(rest)
            
            % you are in /OpenPropDirectory/
            mkdir(['./',filename])
               cd(['./',filename])
               
            addpath ../SourceCode
               
        else
            % you are in /OpenPropDirectory/wrongfolder
            disp('ERROR1: Must start OpenProp from the root directory.')
            return
        end
    end
end
% -------------------------------------------------------------------------


THRUST      = str2double(get(SpecificationsValues(1),'string')); 	% required thrust [N]
Vs          = str2double(get(SpecificationsValues(2),'string'));  % ship velocity [m/s]
Dhub        = str2double(get(SpecificationsValues(3),'string'));  % hub diameter [m]
rho         = str2double(get(SpecificationsValues(4),'string')); 	% water density [kg/m^3]
Mp          = str2double(get(SpecificationsValues(5),'string')); 	% number of vortex panels over the radius
Np          = str2double(get(SpecificationsValues(6),'string')); 	% Number of points over the chord [ ]
Rhub     = Dhub/2;                               % hub radius [m]

Zmin        = str2double(get(RangeValues(1),'string'));  % min number of blades
ZIncrement	= str2double(get(RangeValues(4),'string'));  % increment in number of blades
ZMax        = str2double(get(RangeValues(7),'string'));  % max number of blades

Nmin        = str2double(get(RangeValues(2),'string'));  % min propeller speed [RPM]
NIncrement	= str2double(get(RangeValues(5),'string'));  % increment in propeller speed [RPM]
NMax        = str2double(get(RangeValues(8),'string'));  % max propeller speed [RPM]

Dmin        = str2double(get(RangeValues(3),'string'));  % min propeller diameter [m]
DIncrement	= str2double(get(RangeValues(6),'string'));  % increment in propeller diameter [m]
DMax        = str2double(get(RangeValues(9),'string'));  % max propeller diameter [m]

Propeller_flag	= get(FlagValues(1),'value');               % 0 == turbine, 1 == propeller
Hub_flag  	= get(FlagValues(3),'value');                   % 0 == no hub, 1 == hub

% Duct_flag	= get(FlagValues(4),'value');                   % 0 == no duct, 1 == duct
Duct_flag   = FlagValues(4);

Chord_flag	= get(FlagValues(5),'value');                   % ** CHORD OPTIMIZATION FLAG **

Viscous_flag	= get(FlagValues(6),'value');               % 0 == viscous forces off (CD = 0), 1 == viscous forces on

% Plot_flag       = get(FlagValues(7),'value');               % 0 == do not display plots, 1 == display plots
% Make2Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 2D plot of the results, 1 == make plot
% Make3Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 3D plot of the results, 1 == make plot

Plot_flag       = 0;
Make2Dplot_flag = 0;               % 0 == do not make a 2D plot of the results, 1 == make plot
Make3Dplot_flag = 0;               % 0 == do not make a 3D plot of the results, 1 == make plot

Camber_flag     = 0;                                        % ** VERIFY USE, ADD TO GUI **


XR        = str2double(get(XR_in,'string'));                % radius / propeller radius
XCoD      = str2double(get(XCoD_in,'string'));              % chord / diameter
XCD       = str2double(get(XCD_in,'string'));               % section drag coefficient

ri        = str2double(get(ri_in, 'string'));
VAI       = str2double(get(VAI_in,'string'));               % axial      inflow velocity / ship velocity
VTI       = str2double(get(VTI_in,'string'));               % tangential inflow velocity / ship velocity

%%

Z           = Zmin:ZIncrement:ZMax;
N           = Nmin:NIncrement:NMax;
D           = Dmin:DIncrement:DMax;

dVs         = 0.3/Vs;                                     	% axial inflow variation / Vs

%%

save([filename '_GUIps'],'Zmin','ZIncrement','ZMax','Nmin','NIncrement','NMax',...
    'Dmin','DIncrement','DMax','THRUST','Vs','Dhub','rho','Mp','Np',...
    'Propeller_flag','Hub_flag','Duct_flag',...
    'Chord_flag','Viscous_flag','Plot_flag','Make2Dplot_flag','Make3Dplot_flag',...
    'Camber_flag','filename','XR','XCoD','XCD','VAI','VTI','ri'); % 'Cav_flag','H','dV',

%%


end



function loaddata(hObject,ED)

%%

% =========================================================================
% ================================ loaddata ===============================
% 
% This subfunction loads all values presented in the OpenProp GUI.
%

global OpenPropDirectory SpecificationsValues FlagValues Filename filename RangeValues...
       XR_in XCoD_in XCD_in VAI_in VTI_in;


%%

% ------------------------------------------------------------------------- 
% Figure out what the current working directory is, 
% and change directories to /OpenPropDirectory/
rest = pwd;

while ~isempty(rest)
    [CurrentDirectory,rest] = strtok(rest,'/');
end

if     strcmp(CurrentDirectory,OpenPropDirectory)    % OpenPropDirectory == 'OpenProp_v3.3.4';
    
    % stay in /OpenPropDirectory/
    uiload;
    cd(['./',filename])
    
elseif strcmp(CurrentDirectory,'SourceCode')
    
    cd('../')
    uiload;
    cd(['./',filename])
        
else
    % already in /OpenPropDirectory/filename
    uiload;
end
% -------------------------------------------------------------------------


%%

% % ------------------------------------------------------------------------
% % --- Change directory to Saved Files folder to load, then change back ---
% 
% cd('../OpenProp_saved_files/');
% uiload;
% cd('../SourceCode/');

% ------------------------------------------------------------------------
% ------------------------------------------------------------------------


set(SpecificationsValues(1),'string',num2str()); 	% required thrust [N]
set(SpecificationsValues(2),'string',num2str(Vs));  % ship velocity [m/s]
set(SpecificationsValues(3),'string',num2str(Dhub));  % hub diameter [m]
set(SpecificationsValues(4),'string',num2str(rho)); 	% water density [kg/m^3]
set(SpecificationsValues(5),'string',num2str(Mp)); 	% number of vortex panels over the radius
set(SpecificationsValues(6),'string',num2str(Np)); 	% Number of points over the chord [ ]


Rhub     = Dhub/2;                               % hub radius [m]

set(RangeValues(1),'string',num2str(Zmin));  % min number of blades
set(RangeValues(4),'string',num2str(ZIncrement));  % increment in number of blades
set(RangeValues(7),'string',num2str(ZMax));  % max number of blades

set(RangeValues(2),'string',num2str(Nmin));  % min propeller speed [RPM]
set(RangeValues(5),'string',num2str(NIncrement));  % increment in propeller speed [RPM]
set(RangeValues(8),'string',num2str(NMax));  % max propeller speed [RPM]

set(RangeValues(3),'string',num2str(Dmin));  % min propeller diameter [m]
set(RangeValues(6),'string',num2str(DIncrement));  % increment in propeller diameter [m]
set(RangeValues(9),'string',num2str(DMax));  % max propeller diameter [m]

set(FlagValues(1),'value',Propeller_flag);               % 0 == turbine, 1 == propeller
set(FlagValues(3),'value',Hub_flag);                   % 0 == no hub, 1 == hub

% Duct_flag	= get(FlagValues(4),'value');                   % 0 == no duct, 1 == duct

set(FlagValues(5),'value',Chord_flag);                   % ** CHORD OPTIMIZATION FLAG **

set(FlagValues(6),'value',Viscous_flag);               % 0 == viscous forces off (CD = 0), 1 == viscous forces on
set(FlagValues(7),'value',Plot_flag);               % 0 == do not display plots, 1 == display plots

set(FlagValues(8),'value',Make2Dplot_flag);               % 0 == do not make a 2D plot of the results, 1 == make plot
set(FlagValues(8),'value',Make3Dplot_flag);               % 0 == do not make a 3D plot of the results, 1 == make plot

% Camber_flag     = 0;                                        % ** VERIFY USE, ADD TO GUI **

for index = 1 : length(XR);
    
    set(XR_in(index),'string',num2str(XR(index)));                    % radius / propeller radius
    set(XCoD_in(index),'string',num2str(XCoD(index)));                % chord / diameter
    set(XCD_in(index),'string',num2str(XCD(index)));                  % section drag coefficient
    
    set(VAI_in(index),'string',num2str(VAI(index)));                  % axial      inflow velocity / ship velocity
    set(VTI_in(index),'string',num2str(VTI(index)));                  % tangential inflow velocity / ship velocity
    
end

for index = 1 : length(ri)
    if isnan(ri(index))
        set(ri_in(index),'string','');
    else
        set(ri_in(index),'string',num2str(ri(index)));
    end
    
    if isnan(VAI(index))
        set(VAI_in(index),'string','');
    else
        set(VAI_in(index),'string',num2str(VAI(index)));
    end
    
    if isnan(VTI(index))
        set(VTI_in(index),'string','');
    else
        set(VTI_in(index),'string',num2str(VTI(index)));
    end
end

%%

Z           = Zmin:ZIncrement:ZMax;
N           = Nmin:NIncrement:NMax;
D           = Dmin:DIncrement:DMax;

dVs         = 0.3/Vs;                                     	% axial inflow variation / Vs

%%

set(Filename,'string',filename);                	% Filename prefix

% ----------------------------------------------
% ----------------------------------------------

%%

end


% -------------------------------------------------------------------------
function changeChord(hObject,ED)

%%

% global SpecificationsValues DuctValues CavValues FlagValues FoilValues Filename filename...
%        XR_in XCoD_in XCD_in VAI_in VTI_in Xt0oD_in skew0_in rake0_in...
%        Meanlinecell Thicknesscell;

global FlagValues XCoD_in;

global Col_Label XCoD_values XCLmax_values;

if get(FlagValues(5),'value')
    
    XCoD_values     = str2double(get(XCoD_in,'string'));
    
    set(Col_Label(2),'string','XCLmax');
    
    for i = 1:length(XCoD_in)
        set(XCoD_in(i),'string',num2str(XCLmax_values(i)));
    end
else
    
    XCLmax_values  = str2double(get(XCoD_in,'string'));

    
    set(Col_Label(2),'string','c/D');
    
    for i = 1:length(XCoD_in)
        set(XCoD_in(i),'string',num2str(XCoD_values(i)));
    end        
end

end
% -------------------------------------------------------------------------




function Selectfn(hObject,ED)

%%

global Select;

global  SpecificationsValues FlagValues Filename filename RangeValues...
    XR_in XCoD_in XCD_in VAI_in VTI_in ri_in; % ...

if get(Select,'value')==1
    
    % =========================================================================
    % ================================ savedata ===============================
    %
    % This subfunction saves all values presented in the OpenProp GUI.
    %
    
    filename   	= get(Filename,'string');                       % Filename prefix
    
    THRUST      = str2double(get(SpecificationsValues(1),'string')); 	% required thrust [N]
    Vs          = str2double(get(SpecificationsValues(2),'string'));  % ship velocity [m/s]
    Dhub        = str2double(get(SpecificationsValues(3),'string'));  % hub diameter [m]
    rho         = str2double(get(SpecificationsValues(4),'string')); 	% water density [kg/m^3]
    Mp          = str2double(get(SpecificationsValues(5),'string')); 	% number of vortex panels over the radius
    Np          = str2double(get(SpecificationsValues(6),'string')); 	% Number of points over the chord [ ]
    
    
    Rhub     = Dhub/2;                               % hub radius [m]
    
    Zmin        = str2double(get(RangeValues(1),'string'));  % min number of blades
    ZIncrement	= str2double(get(RangeValues(4),'string'));  % increment in number of blades
    ZMax        = str2double(get(RangeValues(7),'string'));  % max number of blades
    
    Nmin        = str2double(get(RangeValues(2),'string'));  % min propeller speed [RPM]
    NIncrement	= str2double(get(RangeValues(5),'string'));  % increment in propeller speed [RPM]
    NMax        = str2double(get(RangeValues(8),'string'));  % max propeller speed [RPM]
    
    Dmin        = str2double(get(RangeValues(3),'string'));  % min propeller diameter [m]
    DIncrement	= str2double(get(RangeValues(6),'string'));  % increment in propeller diameter [m]
    DMax        = str2double(get(RangeValues(9),'string'));  % max propeller diameter [m]
    
    Propeller_flag	= get(FlagValues(1),'value');               % 0 == turbine, 1 == propeller
    Hub_flag  	= get(FlagValues(3),'value');                   % 0 == no hub, 1 == hub
    
    % Duct_flag	= get(FlagValues(4),'value');                   % 0 == no duct, 1 == duct
    Duct_flag   = FlagValues(4);
    
    Chord_flag	= get(FlagValues(5),'value');                   % ** CHORD OPTIMIZATION FLAG **
    
    Viscous_flag	= get(FlagValues(6),'value');               % 0 == viscous forces off (CD = 0), 1 == viscous forces on

%     Plot_flag       = get(FlagValues(7),'value');               % 0 == do not display plots, 1 == display plots
%     Make2Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 2D plot of the results, 1 == make plot
%     Make3Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 3D plot of the results, 1 == make plot
    
    Plot_flag       = 0;
    Make2Dplot_flag = 0;
    Make3Dplot_flag = 0;
    
    Camber_flag     = 0;                                        % ** VERIFY USE, ADD TO GUI **
    
    XR        = str2double(get(XR_in,'string'));                % radius / propeller radius
    XCoD      = str2double(get(XCoD_in,'string'));              % chord / diameter
    XCLmax    = str2double(get(XCoD_in,'string'));              % CL max

    XCD       = str2double(get(XCD_in,'string'));               % section drag coefficient
    
    ri        = str2double(get(ri_in, 'string'));
    VAI       = str2double(get(VAI_in,'string'));               % axial      inflow velocity / ship velocity
    VTI       = str2double(get(VTI_in,'string'));               % tangential inflow velocity / ship velocity
    
    %%
    
    Z           = Zmin:ZIncrement:ZMax;
    N           = Nmin:NIncrement:NMax;
    D           = Dmin:DIncrement:DMax;
    
    dVs         = 0.3/Vs;                                     	% axial inflow variation / Vs
    
    %%
    
    save('OpenPropTempFile0307122010','Zmin','ZIncrement','ZMax','Nmin','NIncrement','NMax',...
        'Dmin','DIncrement','DMax','THRUST','Vs','Dhub','rho','Mp','Np',...
        'Propeller_flag','Hub_flag','Duct_flag',...
        'Chord_flag','Viscous_flag','Plot_flag','Make2Dplot_flag','Make3Dplot_flag',...
        'Camber_flag','filename','XR','XCoD','XCD','VAI','VTI','ri'); % 'Cav_flag','H','dV',
    
    %%
    
    OpenPropSingle;
    
    
elseif get(Select,'value')==2
    % Do nothing - already in Single Design
elseif get(Select,'value')==3
    OpenPropAnalyze;
end

end



function checkTurbine(hObject,ED)

global FlagValues SpecificationsValues;

if get(FlagValues(1),'value')
    set(SpecificationsValues(1),'enable','on')
else
    % Grays out THRUST field if turbine is selected
    set(SpecificationsValues(1),'string','','enable','off')
end


end



function checkBlades(hObject,ED)

%%

global RangeValues;

Z = floor(str2double(get(RangeValues(1),'string')));

set(RangeValues(1),'string',num2str(Z));

Z = floor(str2double(get(RangeValues(4),'string')));

set(RangeValues(4),'string',num2str(Z));

Z = floor(str2double(get(RangeValues(7),'string')));

set(RangeValues(7),'string',num2str(Z));


end



function changeDir(hObject,ED)

global Filename filename OpenPropDirectory;

%%

filename    = get(Filename,'string');                       % Filename prefix

% ------------------------------------------------------------------------- 
% Figure out what the current working directory is, 
% and change directories to OpenPropDirectory/filename
rest    = pwd;
root  = '';

while ~isempty(rest)
    [CurrentDirectory,rest] = strtok(rest,'/');
    
    if strcmp(CurrentDirectory,OpenPropDirectory)
        
        if isempty(rest)
            
            % you are in /OpenPropDirectory/
            %             addpath ../SourceCode
            
        elseif strcmp(rest(2:end),filename)
            % already in /OpenPropDirectory/filename/
            %             addpath ../SourceCode
            rest = [];
            
        elseif strcmp(rest(2:end),'SourceCode')
            
            %             addpath ../SourceCode
            rest = [];
            
        else
            % you are in /OpenPropDirectory/wrongfolder
%             disp('NOTICE: Must start OpenProp from the root directory.')
            
            cd([root '/' OpenPropDirectory]);
            
%             disp('Changed to root directory.')
            
        end
    end
    
    root  = [root '/' CurrentDirectory];
    
end
% -------------------------------------------------------------------------

end
%%



