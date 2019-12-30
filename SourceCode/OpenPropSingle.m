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
% This code runs the single-design GUI.
%
%--------------------------------------------------------------------------

% =========================================================================
% ========================= Initiate OpenProp ======================
function OpenPropSingle

    clear all;
    clear global;

    % warning off
    addpath ../SourceCode
    % addpath  ./SourceCode


    % =========================================================================
    % ------------------------- Setup global variables ------------------------

    global Fig_Main;    % Main GUI figure

    global Select;

    global OpenPropDirectory SpecificationsValues DuctValues FlagValues FoilValues Filename...
           XR_in XCoD_in XCD_in VAI_in VTI_in ri_in Xt0oD_in skew0_in rake0_in...
           Meanline_cell Thickness_cell Values; % CavValues


    % ------ Variable definitions ------
    %
    % Variables marked with ** are not taken from user input through the GUI.
    % Instead, these variables are computed through formulas or remain unused
    % in this version of the GUI but are still present to allocate space for
    % future use.
    % 
    % % 
    % % global Filename             % Filename prefix
    % % 
    % % 
    % % global SpecificationsValues;      % Contains the following variables:
    % % % General variables
    % % global Z_in;                % Number of blades
    % % global N_in;                % Propeller speed [RPM]
    % % global D_in;                % Propeller diameter [m]
    % % global T_in;                % Required thrust [N]
    % % global Vs_in;               % Ship speed [m/s]
    % % global Dhub_in;             % Hub diameter [m]
    % % 
    % % global Js_in;               % ** Js = Vs/(n*D) ,  advance coefficient
    % % global KT_in;               % ** KT = THRUST/(rho*n^2*D^4)
    % % global n_in;                % ** propeller speed [rev/s] = Vs/(Js*D) = N/60
    % % % Advanced variables
    % % global rho_in;              % Water density
    % % global Mp_in;               % Number of vortex panels over the radius
    % % global Np_in;               % Number of points over the chord
    % % global ITER_in;             % Max iterations
    % % global Rhv_in;              % Hub vortex radius/hub radius
    % %     global SpecificationsText;
    % % 
    % % 
    % % global DuctValues;          % Contains the following variables:
    % % % Ducted Propeller variables
    % % global TAU;                 % Thrust ratio
    % % global CDd;                 % Duct section drag coefficient
    % % 
    % % 
    % % global CavValues;           % Contains the following variables:
    % % % Cavitation Analysis variables
    % % global H_in;                % Shaft centerline depth [m]
    % % global dV_in;               % Inflow variation [m/s]
    % % 
    % % 
    % % global FlagValues;          % Contains the following variables:
    % % % Flags
    % % global Propeller_flag;      % 0 == turbine, 1 == propeller
    % % global Hub_flag;            % 0 == no hub, 1 == hub
    % % global Duct_flag;           % 0 == no duct, 1 == duct
    % % 
    % % global Chord_flag;          % ** CHORD OPTIMIZATION FLAG **
    % % 
    % % global Viscous_flag;        % 0 == viscous forces off (CD = 0), 1 == viscous forces on
    % % global Plot_flag;           % 0 == do not display plots, 1 == display plots
    % % 
    % % global FoilValues;          % Contains the following variables:
    % % % Airfoil type
    % % global Meanline;            % Meanline form
    % % global Thickness;           % Thickness form
    % % 
    % % % Variables for geometry input table
    % % global XR_in;               % Radius ratio (r/R) at control point
    % % global XCoD_in;             % Camber to diameter ratio at XR
    % % global XCD_in;              % Coefficient of drag at XR
    % % global VAI_in;              % Va/Vs (wake profile at XR)
    % % global VTI_in;              % Vt/Vs (wake profile at XR)
    % % global Xt0oD_in;            % Thickness to camber ratio at XR
    % % global skew0_in;            % Skew at XR
    % % global rake0_in;            % Rake at XR

    global N_R0;                % Number of input radii

    global Col_Label;

    global XCoD_values XCLmax_values XCD_values;


    % --- Set GUI element variables as global ---

    % =========================================================================



    % =========================================================================
    % ========================== Initiate Main Figure =========================
     Meanline_cell   = {'NACA a=0.8' 'NACA a=0.8 (modified)' 'Parabolic'};
    Thickness_cell   = {'NACA 65A010' 'NACA 65A010 (modified)' 'Elliptical' 'Parabolic' 'NACA 66 (DTRC modified)'};

    % -------------------- Declare Default Input Variables --------------------

    % Variables for geometry input table
    XR_def          = [.2 .3 .4 .5 .6 .7 .8 .9 .95 1];          % Radius ratio (r/R)
    N_R0            = length(XR_def);                           % Number of input radii

    
    
    XCoD_def        = [0.1600 0.1812 0.2024 0.2196 0.2305 0.2311 0.2173 0.1807 0.1388 0.0010];     % chord to diameter ratio at XR
    XCD_def         = ones(1,N_R0).*0.008;                      % Coefficient of drag at XR

    % VAI_def         = ones(1,N_R0);                             % Va/Vs (wake profile at XR)
    % VTI_def         = zeros(1,N_R0);                            % Vt/Vs (wake profile at XR)

    % t0oc0_def       = [0.2055 0.1553 0.1180 0.09016 0.06960 0.05418 0.04206 0.03321 0.03228 0.03160];     % Thickness to chord ratio at XR

    % Xt0oD_def       = t0oc0_def .* XCoD_def;  % thickness / diameter at XR

      Xt0oD_def       = [0.0329    0.0281    0.0239    0.0198    0.0160    0.0125 0.0091    0.0060    0.0045  0];  % thickness / diameter at XR

    skew0_def       = zeros(N_R0);                              % Skew at XR
    rake0_def       = zeros(N_R0);                              % Rake at XR

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
    T_def           = 25000;                                    % Required thrust [N]
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


    % Cavitation Analysis variables
    % H_def           = 3;                                        % Shaft centerline depth [m]
    % dV_def          = 0.3;                                      % Inflow variation [m/s]


    % % Flags
    % Propeller_flag  = 1;                                        % 0 == turbine, 1 == propeller
    % Hub_flag        = 1;                                        % 0 == no hub, 1 == hub
    % Duct_flag       = 0;                                        % 0 == no duct, 1 == duct
    % 
    % Chord_flag      = 0;                                        % ** CHORD OPTIMIZATION FLAG **
    % 
    % Viscous_flag    = 1;                                        % 0 == viscous forces off (CD = 0), 1 == viscous forces on
    % Plot_flag       = 0;                                        % 0 == do not display plots, 1 == display plots
    % 

    % % Airfoil type
    % Thickness       = 'NACA66 (DTRC Modified)';                 % Thickness form
    % Meanline        = 'NACA a=0.8';                             % Meanline form

    filename        = 'DefaultPropeller';                           % Filename prefix

    % =========================================================================


    % =========================================================================
    % -------------------------- GUI Switching Check --------------------------
    if exist('OpenPropTempFile0307122010.mat','file')

        load('OpenPropTempFile0307122010.mat');

        Z_def       = ceil(mean([Zmin ZMax]));
        N_def       = mean([Nmin NMax]);
        D_def       = mean([Dmin DMax]);
        T_def       = THRUST;
        Vs_def      = Vs;
        Dhub_def    = Dhub;

        rho_def     = rho;
        Mp_def      = Mp;
        Np_def      = Np;

        XR_def      = XR;
        XCoD_def    = XCoD;
        XCD_def     = XCD;

        XCLmax_def  = 0.5 + (0.2-0.5)*(XR_def-XR_def(1))/(XR_def(end)-XR_def(1));  % CLmax distribution

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
    % -------------------------------------------------------------------------
    % =========================================================================


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

    Ductboxht           = 1 + editboxht * 3 + 2;
    Ductbox             = 2 + textbox + editbox + 2;

    BladeDesignboxht	= 1 + editboxht * 11 + 2;
    BladeDesignbox      = 2 + editbox * 6 + 2;

    Inflowboxht    	= 1 + editboxht * 11 + 2;
    Inflowbox      	= 2 + editbox * 3 + 2;

    % Cavboxht        = 1 + editboxht * 3 + 2;
    % Cavbox          = 2 + cavtextbox + editbox + 2;

    Valuesboxht        = 1 + editboxht * 3 + 2;
    Valuesbox          = 2 + cavtextbox + editbox + 2 + cavtextbox + 4 + editbox + 2;

    Flagboxht       = 2 + textht * 8 + 2;
    Flagbox         = 1 + textbox + 1;

    Foilboxht       = BladeDesignboxht - Flagboxht;
    Foilbox         = 1 + textbox + 1;

    Toolboxht       = 1 + pushht + 1 + editboxht + 2;
    Toolbox         = Specificationsbox + BladeDesignbox + Inflowbox + Flagbox - Ductbox - Valuesbox;

    filenamebox     = Toolbox - (2 + filenametext + 2);

    buttonspace     = (Toolbox-pushbox*2-runbox)/4;

    Windowht        = 1 + Ductboxht + Specificationsboxht + 1 + titleht + 1;
    Window          = 1 + Specificationsbox + BladeDesignbox + Inflowbox + Flagbox + 1;

    GUIselectionboxht   = 1 + selectboxht + 1;
    GUIselectionbox     = 1 + selectbox   + 1;
    % =========================================================================


    % =========================================================================
    % ------------------------- Create figure for GUI -------------------------
    close all;

    Fig_Main    = figure('units','characters','position',[5 55-Windowht Window Windowht],...
                         'numbertitle','off','name','OpenProp','menubar','none',...'toolbar','figure',...
                         'resize','off','color',[0.702 0.702 0.702]);

    % % -------------------------------                 
    % if strcmp(computer,'GLNX32') || strcmp(computer,'GLNXA64')

        set(Fig_Main,'resize','on');

    % end
    % % -------------------------------


    % -------------------------------------------------------------------------
    OpenPropDirectory = 'OpenProp_v3.3.4';
    OpenPropVersion   = 'OpenProp v3.3.4';

    Title       = uicontrol(Fig_Main,'style','text','fontsize',titlefontsize,...
                            'fontweight','bold','units','characters','position',...
                            [GUIselectionbox-12 Windowht-1-titleht Window-GUIselectionbox 3],'string',{OpenPropVersion});


    % -------------------------------------------------------------------------
    % --------- Setup panels --------

    % % GUIselection    = uibuttongroup('parent',Fig_Main,'fontsize',panelfontsize,...      %,'title',''
    % %                       'fontweight','bold','units','characters','position',...
    % %                       [1 1+Ductboxht+Specificationsboxht GUIselectionbox GUIselectionboxht],'clipping','on');

    Specifications	= uipanel('parent',Fig_Main,'title','Specifications','fontsize',panelfontsize,...
                          'fontweight','bold','units','characters','position',...
                          [1 1+Ductboxht Specificationsbox Specificationsboxht],'clipping','on');

    Duct            = uipanel('parent',Fig_Main,'title','    Ducted Propeller','fontsize',...
                          panelfontsize,'fontweight','bold','units','characters',...
                          'position',[1 1 Ductbox Ductboxht],'clipping','on');

    BladeDesign     = uipanel('parent',Fig_Main,'title','Blade Design Values','fontsize',...
                          panelfontsize,'fontweight','bold','units','characters',...
                          'position',[1+Ductbox 1+Ductboxht BladeDesignbox BladeDesignboxht],...
                          'clipping','on');

    Inflow          = uipanel('parent',Fig_Main,'title','Inflow Profile Values','fontsize',...
                          panelfontsize,'fontweight','bold','units','characters',...
                          'position',[1+Ductbox+BladeDesignbox 1+Ductboxht Inflowbox Inflowboxht],...
                          'clipping','on');

    % Cav             = uipanel('parent',Fig_Main,'title','    Cavitation Analysis','fontsize',...
    %                       panelfontsize,'fontweight','bold','units','characters',...
    %                       'position',[1+Ductbox 1 Cavbox Cavboxht],'clipping','on');

    Calculator      = uipanel('parent',Fig_Main,'title','Non-dimensional Parameters','fontsize',...
                          panelfontsize,'fontweight','bold','units','characters',...
                          'position',[1+Ductbox 1 Valuesbox Valuesboxht],'clipping','on');

    Flags           = uibuttongroup('parent',Fig_Main,'title','Options','fontsize',...
                          panelfontsize,'fontweight','bold','units','characters',...
                          'position',[1+Specificationsbox+BladeDesignbox+Inflowbox...
                          1+Ductboxht+Foilboxht Flagbox Flagboxht],'clipping','on',...
                          'SelectionChangeFcn',@checkTurbine);

    Foil            = uipanel('parent',Fig_Main,'title','Airfoil type','fontsize',...
                          panelfontsize,'fontweight','bold','units','characters',...
                          'position',[1+Specificationsbox+BladeDesignbox+Inflowbox...
                          1+Ductboxht Foilbox Foilboxht],'clipping','on');

    Tools           = uipanel('parent',Fig_Main,'title','Tools','fontsize',...
                          panelfontsize,'fontweight','bold','units','characters',...
                          'position',[1+Ductbox+Valuesbox 1 Toolbox Toolboxht],...
                          'clipping','on');
    % -------------------------------------------------------------------------


    % -------------------------------------------------------------------------
    % ----------------------- Package Selection Elements ----------------------

    % % Select(1)     = uicontrol(Fig_Main,'units','characters','style',...
    % %                             'pushbutton','string','PS','position',...
    % %                             [1+1 1+Ductboxht+Specificationsboxht+2 selectbox selectboxht],...
    % %                             'horizontalalignment','left','callback','OpenPropParam',...
    % %                             'tooltipstring','Parametric Study');
    % % 
    % % Select(2)     = uicontrol(Fig_Main,'units','characters','style',...
    % %                             'pushbutton','string','SP','position',...
    % %                             [1+1+selectbox 1+Ductboxht+Specificationsboxht+2 selectbox selectboxht],...
    % %                             'horizontalalignment','left','value',1,...
    % %                             'tooltipstring','Single Propeller Design');
    % % 
    % % Select(3)     = uicontrol(Fig_Main,'units','characters','style',...
    % %                             'pushbutton','string','A','position',...
    % %                             [1+1+2*selectbox 1+Ductboxht+Specificationsboxht+2 selectbox selectboxht],...
    % %                             'horizontalalignment','left','callback','OpenPropAnalyze',...
    % %                             'tooltipstring','Analyze Propeller');
    % % 

%     Selectcell      = {'Single Design','Parametric Study','Off-design Analysis'};
    Selectcell      = {'Single Design','Parametric Study'};

    Select          = uicontrol(Fig_Main,'units','characters','style','popupmenu',...
                                'position',[1+1 1+Ductboxht+Specificationsboxht+2 selectbox selectboxht],...
                                'backgroundcolor','w','string',Selectcell,'value',1,'callback',@Selectfn);
    % -------------------------------------------------------------------------


    % -------------------------------------------------------------------------
    % --------------------- Specifications Panel Elements ---------------------

    SpecificationsStrings   = {'Number of blades:'...
                               'Rotation speed (RPM):'...
                               'Rotor diameter (m):'...
                               'Required thrust (N):'...
                               'Ship speed (m/s):'...
                               'Hub diameter (m):'...
                               'Fluid density (kg/m^3):' ...
                               '# radial panels:'...
                               '# chordwise panels:'};

    SpecificationsValues_def      = [Z_def N_def D_def T_def Vs_def Dhub_def rho_def Mp_def Np_def];

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

    % Set callback for those Specs that affect Js and lambda

    set(SpecificationsValues(1),'callback',@checkBlades);
    set(SpecificationsValues(2),'callback',@updateValues);
    set(SpecificationsValues(3),'callback',@updateValues);
    set(SpecificationsValues(4),'callback',@updateValues);
    set(SpecificationsValues(5),'callback',@updateValues);
    set(SpecificationsValues(6),'callback',@updateValues);
    set(SpecificationsValues(7),'callback',@updateValues);
    % -------------------------------------------------------------------------


    % -------------------------- Duct Panel Elements --------------------------

    Ducttext(1)     = uicontrol(Duct,'units','characters','style',...
                                'text','string','Thrust Ratio:','position',...
                                [2 Ductboxht-4 textbox textht],...
                                'horizontalalignment','left');

    Ducttext(2)     = uicontrol(Duct,'units','characters','style',...
                                'text','string','Duct section drag (Cd):','position',...
                                [2 Ductboxht-4-editboxht textbox textht],...
                                'horizontalalignment','left');

    Ducttext(3)     = uicontrol(Duct,'units','characters','style',...
                                'text','string','duct D / prop D:','position',...
                                [2 Ductboxht-4-editboxht*2 textbox textht],...
                                'horizontalalignment','left');

    DuctValues(1)   = uicontrol(Duct,'units','characters','style',...
                                'edit','string',num2str(TAU_def),...
                                'position',[2+textbox Ductboxht-4,...
                                editbox editboxht],'backgroundcolor',[1 1 1],...
                               'enable','off');

    DuctValues(2)   = uicontrol(Duct,'units','characters','style',...
                                'edit','string',num2str(CDd_def),...
                                'position',[2+textbox Ductboxht-4-editboxht,...
                                editbox editboxht],'backgroundcolor',[1 1 1],...
                               'enable','off');

    DuctValuesoff(1)= uicontrol(Duct,'units','characters','style',...
                                'edit','string','1',...
                                'position',[2+textbox Ductboxht-4-editboxht*2,...
                                editbox editboxht],'backgroundcolor',[1 1 1],...
                                'enable','off');
    % -------------------------------------------------------------------------


    % ----------------------- Blade Design Panel Elements ---------------------

    ColName     = {'r/R' 'c/D' 'Cd' 't0/D' 'Skew' 'Xs/D'};

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

    % 	f0oc_in(index)  = uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
    %                                 'backgroundcolor','w','string',num2str(f0oc_def(index)),'position',...
    %                                 [2+editbox*5 BladeDesignboxht-4-editboxht*index editbox editboxht]);

        Xt0oD_in(index)	= uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
                                     'backgroundcolor','w','string',num2str(Xt0oD_def(index)),'position',...
                                     [2+editbox*3 BladeDesignboxht-4-editboxht*index editbox editboxht]);

        skew0_in(index)	= uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
                                     'backgroundcolor','w','string',num2str(skew0_def(index)),'position',...
                                     [2+editbox*4 BladeDesignboxht-4-editboxht*index editbox editboxht]);

        rake0_in(index)	= uicontrol(BladeDesign,'style','edit','units','characters','FontSize',10,...
                                     'backgroundcolor','w','string',num2str(rake0_def(index)),'position',...
                                     [2+editbox*5 BladeDesignboxht-4-editboxht*index editbox editboxht]);

    end
    % -------------------------------------------------------------------------


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

    % set(VAI_in(1),'string','1');
    % set(VTI_in(1),'string','0');
    % -------------------------------------------------------------------------


    % ------------------- Cavitation Analysis Panel Elements ------------------

    % Cavtext(1)     = uicontrol(Cav,'units','characters','style',...
    %                            'text','string','Shaft centerline depth (m):','position',...
    %                            [2 Cavboxht-4 cavtextbox textht],...
    %                            'horizontalalignment','left');
    % 
    % Cavtext(2)     = uicontrol(Cav,'units','characters','style',...
    %                            'text','string','Inflow variation (m/s):','position',...
    %                            [2 Cavboxht-4-editboxht cavtextbox textht],...
    %                            'horizontalalignment','left');
    % 
    % Cavtext(3)     = uicontrol(Cav,'units','characters','style',...
    %                            'text','string','Ideal angle of attackii:','position',...
    %                            [2 Cavboxht-4-editboxht*2 cavtextbox textht],...
    %                            'horizontalalignment','left');
    % 
    % CavValues(1)   = uicontrol(Cav,'units','characters','style',...
    %                            'edit','string',num2str(TAU_def),...
    %                            'position',[2+cavtextbox Cavboxht-4,...
    %                            editbox editboxht],'backgroundcolor',[1 1 1],...
    %                            'enable','off');
    % 
    % CavValues(2)   = uicontrol(Cav,'units','characters','style',...
    %                            'edit','string',num2str(CDd_def),...
    %                            'position',[2+cavtextbox Cavboxht-4-editboxht,...
    %                            editbox editboxht],'backgroundcolor',[1 1 1],...
    %                            'enable','off');
    % 
    % CavValuesoff(1)= uicontrol(Cav,'units','characters','style',...
    %                            'edit','string','1',...
    %                            'position',[2+cavtextbox Cavboxht-4-editboxht*2,...
    %                            editbox editboxht],'backgroundcolor',[1 1 1],...
    %                            'enable','off');
    % -------------------------------------------------------------------------


    % -------------------------- Values Panel Elements ------------------------

    Valuestext(1)     = uicontrol(Calculator,'units','characters','style',...
                               'text','string','J = V/nD =','position',...
                               [2 Valuesboxht-4 cavtextbox textht],...
                               'horizontalalignment','left');

    Valuestext(2)     = uicontrol(Calculator,'units','characters','style',...
                               'text','string','L = omega*R/V =','position',...
                               [2 Valuesboxht-4-editboxht cavtextbox textht],...
                               'horizontalalignment','left');

    Valuestext(3)     = uicontrol(Calculator,'units','characters','style',...
                               'text','string','KT = T/(rho*n^2*D^4) =','position',...
                               [2+cavtextbox+editbox+2 Valuesboxht-4-editboxht cavtextbox+4 textht],...
                               'horizontalalignment','left');

    Valuestext(4)     = uicontrol(Calculator,'units','characters','style',...
                               'text','string','CT = T/(1/2*rho*V^2*pi*R^2) =','position',...
                               [2+cavtextbox+editbox+2 Valuesboxht-4 cavtextbox+4 textht],...
                               'horizontalalignment','left');

    Values(1)   = uicontrol(Calculator,'units','characters','style',...
                               'edit','string',Js_def,...
                               'position',[cavtextbox-4, Valuesboxht-4, editbox, editboxht],'backgroundcolor',[1 1 1],...
                               'enable','off');

    Values(2)   = uicontrol(Calculator,'units','characters','style',...
                               'edit','string',lambda_def,...
                               'position',[cavtextbox-4, Valuesboxht-4-editboxht, editbox, editboxht],'backgroundcolor',[1 1 1],...
                               'enable','off');

    Values(3)   = uicontrol(Calculator,'units','characters','style',...
                               'edit','string',KT_def,...
                               'position',[2+cavtextbox*2+editbox+2+4 Valuesboxht-4-editboxht,...
                               editbox editboxht],'backgroundcolor',[1 1 1],...
                               'enable','off');

    Values(4)   = uicontrol(Calculator,'units','characters','style',...
                               'edit','string',CT_def,...
                               'position',[2+cavtextbox*2+editbox+2+4 Valuesboxht-4,...
                               editbox editboxht],'backgroundcolor',[1 1 1],...
                               'enable','off');

    %Values(3)= uicontrol(Calculator,'units','characters','style',...
    %                            'edit','string','1',...
    %                            'position',[2+cavtextbox Valuesboxht-4-editboxht*2,...
    %                            editbox editboxht],'backgroundcolor',[1 1 1],...
    %                            'enable','off');

    % -------------------------------------------------------------------------


    % -------------------------------------------------------------------------
    % -------------------------- Flags Panel Elements --------------------------

    FlagValues(1)     = uicontrol(Flags,'units','characters','style','radiobutton',...
                                  'string','Propeller','value',1,'position',[1 Flagboxht-4 textbox textht]);

    FlagValues(2)     = uicontrol(Flags,'units','characters','style','radiobutton',...
                                  'string','Turbine','position',[1 Flagboxht-4-textht textbox textht]);

    FlagValues(3)     = uicontrol(Flags,'units','characters','style','checkbox',...
                                  'string','Hub','value',1,'position',[1 Flagboxht-4-textht*2 textbox textht]);

    % ---------------------------------------------------
    % ------ Ducted flag relocated to Ducted panel ------

    % FlagValues(4)     = uicontrol(Flags,'units','characters','style','checkbox',...
    %                               'string','Ducted','position',...
    %                               [1 Flagboxht-4-textht*3 textbox textht]);

    FlagValues(4)     = uicontrol(Duct,'units','characters','style','checkbox',...
                                  'string','','position',[1 Ductboxht-textht...
                                  textbox textht],'callback',@changeDuct);
    % ---------------------------------------------------

    FlagValues(5)     = uicontrol(Flags,'units','characters','style','checkbox',...
                                  'string','Chord optimization','position',[1 Flagboxht-4-textht*3 textbox textht],...
                                  'callback',@changeChord);

    FlagValues(6)     = uicontrol(Flags,'units','characters','style','checkbox',...
                                  'string','Viscous forces','value',1,'position',[1 Flagboxht-4-textht*4 textbox textht],...
                                  'callback',@changeViscous);

    FlagValues(7)     = uicontrol(Flags,'units','characters','style','checkbox',...
                                  'string','Optimization plots','position',[1 Flagboxht-4-textht*6 textbox textht]);

    FlagValues(8)     = uicontrol(Flags,'units','characters','style','checkbox',...
                                  'string','Geometry plots','value',1,'position',[1 Flagboxht-4-textht*7 textbox textht]);


    FlagValues(9)     = uicontrol(Flags,'units','characters','style','checkbox',...
                                  'string','Performance curve','value',0,'position',[1 Flagboxht-4-textht*8 textbox textht]);

    % ------ Added Cav_flag ------

    % FlagValues(10)     = uicontrol(Cav,'units','characters','style','checkbox',...
    %                               'string','','position',[1 Cavboxht-textht...
    %                               textbox textht],'callback',@changeCav);

    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------

    % ---------------------- Airfoil Type Panel Elements ----------------------
    FoilText(1)     = uicontrol(Foil,'units','characters','style','text',...
                                'string','Meanline type:','position',...
                                [1 Foilboxht-3.5 textbox textht]);

    FoilValues(1)	= uicontrol(Foil,'units','characters','style','popupmenu',...
                                'position',[1 Foilboxht-3.5-textht textbox textht],...
                                'backgroundcolor','w','string',...
                                Meanline_cell);

    FoilText(2)     = uicontrol(Foil,'units','characters','style','text',...
                                'string','Thickness type:','position',...
                                [1 Foilboxht-3.5-textht*2 textbox textht]);

    FoilValues(2)	= uicontrol(Foil,'units','characters','style','popupmenu',...
                                'position',[1 Foilboxht-3.5-textht*3 textbox textht],...
                                'backgroundcolor','w','string',...
                                Thickness_cell);
    % -------------------------------------------------------------------------


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
    % -------------------------------------------------------------------------

end
%%
% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
% -------------------------------------------------------------------------
%               And the meek shall inherit the earth...                   %         
% -------------------------------------------------------------------------
% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
% =========================================================================
%%
% ================================ execute ================================
function execute(hObject,ED)
    % 
    % This subfunction runs OpenProp.
    %
    newPlots;

    global Plots PlotPanels Toggle OnDesignValues ConversionValues systemToggle;

    global OpenPropDirectory SpecificationsValues DuctValues FlagValues FoilValues Filename...
           XR_in XCoD_in XCD_in VAI_in ri_in VTI_in Xt0oD_in skew0_in rake0_in...
           Meanline_cell Thickness_cell; % CavValues

    global pt

    % ------------------------------------------------------------------------- 
    % Figure out what the current working directory is, 
    % and change directories to OpenPropDirectory/filename
    filename = get(Filename,'string');                       % Filename prefix

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

    Z           = str2double(get(SpecificationsValues(1),'string'));  % number of blades
    N           = str2double(get(SpecificationsValues(2),'string'));  % propeller speed [RPM]
    D           = str2double(get(SpecificationsValues(3),'string'));	% propeller diameter [m]
    THRUST      = str2double(get(SpecificationsValues(4),'string')); 	% required thrust [N]
    Vs          = str2double(get(SpecificationsValues(5),'string'));  % ship velocity [m/s]
    Dhub        = str2double(get(SpecificationsValues(6),'string'));  % hub diameter [m]
    rho         = str2double(get(SpecificationsValues(7),'string')); 	% water density [kg/m^3]
    Mp          = str2double(get(SpecificationsValues(8),'string')); 	% number of vortex panels over the radius
    Np          = str2double(get(SpecificationsValues(9),'string')); 	% Number of points over the chord [ ]
    
    ITER        = 40;   % number of iterations in analysis
    Rhv         = 0.5;	% hub vortex radius / hub radius

    TAU         = str2double(get(DuctValues(1),'string'));      % Thrust ratio
    CDd         = str2double(get(DuctValues(2),'string'));      % Duct section drag coefficient

    % H           = str2double(get(CavValues(1),'string'));       % Shaft centerline depth [m]
    % dV          = str2double(get(CavValues(2),'string'));       % Inflow variation [m/s]


    % --------------------------------- Flags ---------------------------------

    Propeller_flag	= get(FlagValues(1),'value');               % 0 == turbine, 1 == propeller
    Hub_flag  	    = get(FlagValues(3),'value');                   % 0 == no hub, 1 == hub
    Duct_flag	    = get(FlagValues(4),'value');                   % 0 == no duct, 1 == duct

    Chord_flag	    = get(FlagValues(5),'value');                   % ** CHORD OPTIMIZATION FLAG **

    Viscous_flag	= get(FlagValues(6),'value');               % 0 == viscous forces off (CD = 0), 1 == viscous forces on
    Plot_flag       = get(FlagValues(7),'value');               % 0 == do not display plots, 1 == display plots



    Make2Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 2D plot of the results, 1 == make plot
    Make3Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 3D plot of the results, 1 == make plot

    Analyze_flag	= get(FlagValues(9),'value');

    % Cav_flag	= get(FlagValues(10),'value');                   % 0 == do not run cavitation mapping, 1 == run cavitation mapping


    % -------------------------------------------------------------------------
    % ---------------------- Blade 2D section properties ----------------------

    Meanline_index  = get(FoilValues(1),'value');
    Meanline        = char(Meanline_cell(Meanline_index));              % Meanline form

    Thickness_index	= get(FoilValues(2),'value');
    Thickness       = char(Thickness_cell(Thickness_index));            % Thickness form


    XR          = str2double(get(XR_in,  'string'));             	% radius / propeller radius
    
    XCoD        = str2double(get(XCoD_in,'string'));           	% chord / diameter
    XCLmax      = str2double(get(XCoD_in,'string'));            % maximum lift coefficient (for chord optimization)
    
    XCD     	= str2double(get(XCD_in, 'string'));            	% section drag coefficient

    ri          = str2double(get(ri_in, 'string'));
    VAI         = str2double(get(VAI_in,'string'));             % axial      inflow velocity / ship velocity
    VTI         = str2double(get(VTI_in,'string'));          	% tangential inflow velocity / ship velocity

    ri  = ri (~isnan(ri));
    VAI = VAI(~isnan(VAI));
    VTI = VTI(~isnan(VTI));


    % f0oc0     = str2double(get(f0oc_in,'String'));          	% max section camber    / chord
    Xt0oD       = str2double(get(Xt0oD_in,'String'));           % max section thickness / chord
    skew0       = str2double(get(skew0_in,'String'));          	% skew
    rake0       = str2double(get(rake0_in,'String'));         	% rake


    % ----------------------- Compute derived quantities ----------------------

    n           = N/60;                                         % ** propeller speed [rev/s] = Vs/(Js*D) = N/60
    Js          = Vs/(n*D);                                     % ** Js = Vs/(n*D) ,  advance coefficient
    KT          = THRUST/(rho*n^2*D^4);                        	% ** KT = THRUST/(rho*n^2*D^4)
    L           = pi/Js;                                        % tip speed ratio

    R           = D/2;                                          % propeller radius [m]
    Rhub        = Dhub/2;                                       % hub radius [m]
    Rhub_oR     = Rhub/R;

    CTDES       = THRUST/(0.5*rho*Vs^2*pi*R^2);                 % CT thrust coefficient required          

    % dVs         = dV/Vs;                                        % axial inflow variation / Vs

    
    
    % *************************************************************************
    % *************************************************************************
    input.part1      = '------ Performance inputs ------';
    input.Z          = Z;           % [1 x 1], [ ] number of blades
    input.N          = N;           % propeller speed [RPM]
    input.D          = D;           % propeller diameter [m]  
    input.Vs         = Vs;          % [1 x 1], [m/s] ship speed
    input.Js         = Js;          % [1 x 1], [ ] advance coefficient, Js = Vs/nD = pi/L
    input.L          = L;           % [1 x 1], [ ] tip speed ratio, L = omega*R/V
    input.THRUST     = THRUST;      % required thrust [N]
    input.CTDES      = CTDES;       % [1 x 1], [ ] desired thrust coefficient
    input.TAU        = TAU;          % Thrust ratio
    
    input.part2      = '------ Geometry inputs ------';
    input.Mp         = Mp;          % [1 x 1], [ ] number of blade sections
    input.Np         = Np;          % [1 x 1], [ ] number of points along the chord
    input.R          = R;           % [1 x 1], [m] propeller radius
    input.Rhub       = Rhub;        % [1 x 1], [m] hub radius
    input.XR         = XR;          % [length(XR) x 1], [ ] input radius/propeller radiusat XR
    input.XCD        = XCD;         % [length(XR) x 1], [ ] input drag coefficient       at XR
    input.XCoD       = XCoD;        % [length(XR) x 1], [ ] input chord / diameter       at XR
    input.Xt0oD      = Xt0oD;       % [length(XR) x 1], [ ] input thickness / chord      at XR 
    input.skew0      = skew0;       % [length(XR) x 1], [ ] input skew  [deg]      at XR 
    input.rake0      = rake0;       % [length(XR) x 1], [ ] input rake X/D       at XR 
    input.Meanline   = Meanline;    % 2D section meanline  flag
    input.Thickness  = Thickness;   % 2D section thickness flag 
    input.XCLmax     = XCLmax;

    if ~isempty(ri)  , input.ri  = ri;  end
    if ~isempty(VAI) , input.VAI = VAI; end        % [length(XR) x 1], [ ] input axial inflow velocity  at XR
    if ~isempty(VTI) , input.VTI = VTI; end        % [length(XR) x 1], [ ] input swirl inflow velocity  


    input.Rduct     = R;
    input.Cduct     = D/2;
    input.CDd       = CDd;

    input.part3      = '------ Computational inputs ------';
    input.Propeller_flag  = Propeller_flag; % 0 == turbine, 1 == propeller
    input.Viscous_flag    = Viscous_flag;   % 0 == viscous forces off (CD = 0), 1 == viscous forces on
    input.Hub_flag        = Hub_flag;       % 0 == no hub, 1 == hub
    input.Duct_flag       = Duct_flag;      % 0 == no duct, 1 == duct
    input.Plot_flag       = Plot_flag;      % 0 == do not display plots, 1 == display plots
    input.Chord_flag      = Chord_flag;     % 0 == do not optimize chord lengths, 1 == optimize chord lengths

    input.Make2Dplot_flag = Make2Dplot_flag;
    input.Make3Dplot_flag = Make3Dplot_flag;
    % input.Make_Rhino_flag = Make_Rhino_flag;
    input.ITER            = ITER;           % [ ] number of iterations
    input.Rhv              = Rhv;         % [1 x 1], [ ] hub vortex radius / hub radius

    input.part4      = '------ Cavitation inputs ------';
    input.rho        = rho;         % [1 x 1], [kg/m^3] fluid density
    % input.dVs        = dVs;         % [1 x 1], [ ] ship speed variation / ship speed
    % input.H          = H;           % [1 x 1]

    input.part5      = '------ Duct inputs ------';


    % ---------------------------- Pack up propeller/turbine data structure, pt
    pt.filename = filename; % (string) propeller/turbine name
    pt.date     = date;     % (string) date created
  % pt.notes    = ' ';    % (string or cell matrix)   notes
    pt.input    = input;    % (struct) input parameters
    pt.design   = [];       % (struct) design conditions
    pt.geometry = [];       % (struct) design geometry
    pt.states	= [];       % (struct) off-design state analysis

    pt.input.GUI_flag = 1;

    % *************************************************************************
    % *************************************************************************

    
    
    % =========================================================================
    % ============================ execution script ===========================
    
    % ---------------------------------------------------------------------
    % Plot from input:

    % Expanded Blade, only if Chord_flag = 0
    if Chord_flag == 0

        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotPanels(1));
        axes(h);
        hold on;

        XXR   = XR(1) + (XR(end)-XR(1))*(sin((0:60)*pi/(2*60))); 
        XXCoD = InterpolateChord(XR,XCoD,XXR);

        plot(XXR, XXCoD,'b','LineWidth',2);
        plot(XXR,-XXCoD,'b','LineWidth',2);

        plot(XR, XCoD,'.b','MarkerSize',16)
        plot(XR,-XCoD,'.b','MarkerSize',16)

        xlabel('r/R','Fontsize',16,'FontName','Times');   
        ylabel('c/R','Fontsize',16,'FontName','Times');       
        set(gca,'Fontsize',14,'FontName','Times')
        grid on, box on,

    else
        set(Toggle(1),'enable','off');
    end

    % Thickness profile:
        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotPanels(2));
        axes(h);
        hold on; 

        XXR    = XR(1) + (XR(end)-XR(1))*(sin((0:60)*pi/(2*60))); 
        XXt0oD = pchip(XR,Xt0oD,XXR);

        plot(XXR, XXt0oD,'b','LineWidth',2);
        plot(XXR,-XXt0oD,'b','LineWidth',2);

        plot(XR, Xt0oD,'.b','MarkerSize',16)
        plot(XR,-Xt0oD,'.b','MarkerSize',16)

        xlabel('r/R','FontSize',16,'FontName','Times');   
        ylabel('t0/D','FontSize',16,'FontName','Times');       
        set(gca,     'FontSize',14,'FontName','Times');
        grid on, box on,

    % Inflow profile:
        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotPanels(3));
        axes(h);
        hold on;

        if ~isempty(VAI)

            plot(VAI,ri,'b','LineWidth',2);
            plot(VTI,ri,'r','LineWidth',2);

            xlim([-0.1 1.1*max(VAI)])
        else
            plot( ones(size(XR)),XR,'b','LineWidth',2);
            plot(zeros(size(XR)),XR,'r','LineWidth',2);

            xlim([-0.1 1.1])
        end

        xlabel('VA / Vs (blue),   VT / Vs (red)','FontSize',16,'FontName','Times');   
        ylabel('r / R','FontSize',16,'FontName','Times');       
        set(gca,     'FontSize',14,'FontName','Times');
        grid on, box on,

    
    % ---------------------------------------------------------------------
    % Perform design optimization
    pt.design   = EppsOptimizer(input);
    % ---------------------------------------------------------------------

   
    % ---------------------------------------------------------------------
    % Set On Design Performance values

    pt.design.Q = pt.design.CQ * 0.5 * rho * Vs^2 * pi*D^2/4 * D/2; % [Nm]  torque

    omega = 2*pi*n; % [rad/s]

    pt.design.P = pt.design.Q * omega;

    if Propeller_flag == 1
        set(OnDesignValues(1),'string',num2str(pt.design.Js));
        set(OnDesignValues(2),'string',num2str(pt.design.KT));
        set(OnDesignValues(3),'string',num2str(pt.design.KQ));
        set(OnDesignValues(4),'string',num2str(pt.design.EFFY));
    else
        set(OnDesignValues(1),'string',num2str(pi/pt.design.L));
        set(OnDesignValues(2),'string',' ');
        set(OnDesignValues(3),'string',' ');
        set(OnDesignValues(4),'string',' ');
    end
%     set(OnDesignValues(2),'string',num2str(pt.design.KT));
%     set(OnDesignValues(3),'string',num2str(pt.design.KQ));
%     set(OnDesignValues(4),'string',num2str(pt.design.EFFY));

    if Propeller_flag == 1
        set(OnDesignValues(5),'string',num2str(pt.design.ADEFFY));
    end

    set(OnDesignValues(6),'string',num2str(pt.design.CT));
    set(OnDesignValues(7),'string',num2str(pt.design.CQ));
    set(OnDesignValues(8),'string',num2str(pt.design.CP));

    set(ConversionValues(1),'string',num2str(pt.input.Vs));
    set(ConversionValues(2),'string',num2str(pt.input.N));
    set(ConversionValues(3),'string',num2str(pt.input.D));
    set(ConversionValues(4),'string',num2str(pt.input.THRUST));
    set(ConversionValues(5),'string',num2str(pt.design.Q));
    set(ConversionValues(6),'string',num2str(pt.design.P));

    set(systemToggle,'enable','on');
    % ---------------------------------------------------------------------

    
    
    % ---------------------------------------------------------------------    
    disp(' ')
    disp('Creating graphical and text reports')
    disp(' ')

    % Create graphical and text reports
    Make_Reports(pt);
    % ---------------------------------------------------------------------
    
   
    % ---------------------------------------------------------------------
    % Plot t0/D vs RC
        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotPanels(11));
        axes(h);
        hold on;

        plot([Rhub,pt.design.RC,1],interp1(pt.design.RC, pt.design.t0oD,[Rhub,pt.design.RC,1],'spline','extrap'),'b','LineWidth',2);
        plot([Rhub,pt.design.RC,1],interp1(pt.design.RC,-pt.design.t0oD,[Rhub,pt.design.RC,1],'spline','extrap'),'b','LineWidth',2);

        plot(pt.design.RC, pt.design.t0oD,'b.','MarkerSize',16);
        plot(pt.design.RC,-pt.design.t0oD,'b.','MarkerSize',16);

        xlabel('r/R','FontSize',16,'FontName','Times');   
        ylabel('t0/D','FontSize',16,'FontName','Times');            
        set(gca,     'FontSize',14,'FontName','Times');
        grid on; box on,


    % ---------------------------------------------------------------------
    % Plot CL vs RC
        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotPanels(12));
        axes(h);
        hold on;

        plot([Rhub,pt.design.RC,1],interp1(pt.design.RC, pt.design.CL,[Rhub,pt.design.RC,1],'spline','extrap'),'b','LineWidth',2);

        plot(pt.design.RC,pt.design.CL,'b.','MarkerSize',16);

        xlabel('r/R','FontSize',16,'FontName','Times');   
        ylabel('CL','FontSize',16,'FontName','Times');            
        set(gca,     'FontSize',14,'FontName','Times');
        grid on; box on,


            % Set lower y limit to 0
            ylimits = get(gca,'Ylim');

            if abs(ylimits(2)) > abs(ylimits(1))

                set(gca,'Ylim',[0 ylimits(2)]);
            else
                set(gca,'Ylim',[ylimits(1) 0]);
            end   
    % ---------------------------------------------------------------------

    
    % ---------------------------------------------------------------------
    % Determine propeller geometry
    if Make2Dplot_flag == 1 | Make3Dplot_flag == 1
        pt.geometry = Geometry(pt);
    end
    % ---------------------------------------------------------------------
    
    % ---------------------------------------------------------------------
    % if Cav_flag
    %     
    %     % Pefrorm cavitation analysis
    %     Cav_CavitationMap(pt);
    %     
    %     VLMbucket
    %     
    % end
    % ---------------------------------------------------------------------
    
     
    % ---------------------------------------------------------------------
    if Analyze_flag == 1
        % % Analyze off-design states
        % Js_all      = [1.05:-0.1:0.55];     % advance coefficient
        % LAMBDAall   = pi./Js_all;           % tip-speed ratio
        % pt.states	= Analyze(pt,LAMBDAall)

        pt.states      = AnalyzeAuto(pt);

        % ---------------------------------------------------------------------
        set(0,'CurrentFigure',Plots);
        h = axes('parent',PlotPanels(15));
        axes(h);
        hold on;

        if Propeller_flag == 1

            VMIV = pt.design.VMIV;

            Js_curve = linspace(min(pt.states.Js),max(pt.states.Js),100);

            EFFY_curve = (Js_curve/(2*pi)) * VMIV .* pchip(pt.states.Js,pt.states.KT,Js_curve)./pchip(pt.states.Js,pt.states.KQ,Js_curve); 

            % Efficiency (green squares)
                    plot(Js_curve,EFFY_curve,'-','LineWidth',2,'Color',[0 0.8 0])
            Heffy = plot(pt.states.Js,pt.states.EFFY,'sk','MarkerSize',5,'LineWidth',1,'MarkerFaceColor',[0 0.8 0]); 

            % Thrust coefficient (blue diamonds)
                    plot(pt.states.Js,pt.states.KT,'b-','LineWidth',2)
            Hkt   = plot(pt.states.Js,pt.states.KT,'dk','MarkerSize',5,'LineWidth',1,'MarkerFaceColor','b');

            % Torque coefficient (red circles)
                  plot(pt.states.Js,10*pt.states.KQ,'r-','LineWidth',2)    
            Hkq = plot(pt.states.Js,10*pt.states.KQ,'ok','MarkerSize',5,'LineWidth',1,'MarkerFaceColor','r');

            % Design point
            plot(pt.design.Js*[1 1],[0 2],'k--','LineWidth',1);

            xlabel('Js','FontSize',16,'FontName','Times'), 
            ylabel('KT, 10*KQ, EFFY','FontSize',16,'FontName','Times')
            axis([min(pt.states.Js) max(pt.states.Js) 0 0.9])
            set(gca,     'FontSize',14,'FontName','Times');
            box on, grid on,


            ylimits = get(gca,'Ylim');
            set(gca,'Ylim', [0  max(1,ylimits(2)) ] );

        else

            % Power coefficient (blue dots)
            plot(pt.states.L,-pt.states.CP,'b.-','LineWidth',2,'MarkerSize',12)

            % Design point
            plot(pt.design.L*[1 1],[0 0.6],'k--','LineWidth',1);

            % % Betz limit
            % plot([0 ceil(max(pt.states.L))],(16/27)*[1 1],'k--','LineWidth',1);

            xlabel('L','FontSize',16,'FontName','Times'), 
            ylabel('CP','FontSize',16,'FontName','Times'),

            set(gca,'Ytick',[0:0.1:0.6])

            axis([0 ceil(max(pt.states.L))  0 0.6])
            set(gca,     'FontSize',14,'FontName','Times');
            box on, grid on,   

        end

    end
    % ---------------------------------------------------------------------
    
  

    % ---------------------------------------------------------------------
    % Plot parametric study results
    if exist([filename '.mat'],'file')
        
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

        temp    = pt;

        load([filename '.mat']);

        if isfield(pt,'paroutput')

            paroutput   = pt.paroutput;
            parinput    = pt.parinput;

            % --- For EFFY vs N ---

            set(0,'CurrentFigure',Plots);
            h = axes('parent',PlotPanels(4));
            axes(h);

            hold on, box on, grid on,
            ylim([0 1]),

            tempEFFY = paroutput.EFFY(1,:,1);

            plot(paroutput.N,tempEFFY(:),'-','color',CLR( (mod(0,11)+1) ,:));

            str_legend ={[num2str(parinput.D),' m  ']};


            legend(str_legend,'location','southwest');

            xlabel('Rotation Speed (RPM)  ');
            ylabel('Efficiency');
            title(['Number of Blades: ',num2str(Z),'  '])

            % --- For EFFY vs D ---

            set(0,'CurrentFigure',Plots);
            h = axes('parent',PlotPanels(5));
            axes(h);

            hold on, box on, grid on,
            ylim([0 1]),

            tempEFFY = paroutput.EFFY(1,1,:);

            plot(paroutput.D,tempEFFY(:),'-','color',CLR( (mod(0,11)+1) ,:));

            str_legendD ={[num2str(parinput.N),' RPM  ']};

            legend(str_legendD,'location','southwest');

            xlabel('Propeller Diameter (m)  ');
            ylabel('Efficiency');
            title(['Number of Blades: ',num2str(Z),'  '])

            % =========
            % =========

        else

            set(Toggle(4),'enable','off');
            set(Toggle(5),'enable','off');

        end

        pt      = temp;

    else

        set(Toggle(4),'enable','off');
        set(Toggle(5),'enable','off');

    end
    % ---------------------------------------------------------------------
    
    
    % ---------------------------------------------------------------------
    % pt overwrite sequence:

    if exist([filename '.mat'],'file')

        disp(['Found original file:  ',filename,'.mat']);

        temp  = pt;

        load(filename)

        disp(['Overwriting file:  ',filename,'.mat']);

        pt.filename = temp.filename;
        pt.date     = temp.date;
        pt.input	= temp.input;
        pt.design   = temp.design;
        pt.geometry = temp.geometry;
        pt.states   = temp.states;

    end
    
    pt

    save(filename,'pt');

    save([filename,'_GUIsd'],'Z','N','D','THRUST','Vs','Dhub','rho','Mp','Np',...
        'TAU','CDd','Propeller_flag','Hub_flag','Duct_flag',...
        'Chord_flag','Viscous_flag','Plot_flag','Make2Dplot_flag','Make3Dplot_flag','Analyze_flag',...
        'Meanline','Meanline_index','Thickness',...
        'Thickness_index','filename','XR','XCoD','XCD','ri','VAI','VTI','Xt0oD','skew0',...
        'rake0');
    % ---------------------------------------------------------------------
    
end
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% ================================ savedata ===============================
function savedata(hObject,ED)
    % 
    % This subfunction saves all values presented in the OpenProp GUI.
    %

    global OpenPropDirectory SpecificationsValues DuctValues FlagValues FoilValues Filename...
           XR_in XCoD_in XCD_in VAI_in VTI_in ri_in Xt0oD_in skew0_in rake0_in...
           Meanline_cell Thickness_cell XCoD_values XCLmax_values; % CavValues

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


    Z           = str2double(get(SpecificationsValues(1),'string'));  % number of blades
    N           = str2double(get(SpecificationsValues(2),'string'));  % propeller speed [RPM]
    D           = str2double(get(SpecificationsValues(3),'string'));	% propeller diameter [m]
    THRUST      = str2double(get(SpecificationsValues(4),'string')); 	% required thrust [N]
    Vs          = str2double(get(SpecificationsValues(5),'string'));  % ship velocity [m/s]
    Dhub        = str2double(get(SpecificationsValues(6),'string'));  % hub diameter [m]
    rho         = str2double(get(SpecificationsValues(7),'string')); 	% water density [kg/m^3]
    Mp          = str2double(get(SpecificationsValues(8),'string')); 	% number of vortex panels over the radius
    Np          = str2double(get(SpecificationsValues(9),'string')); 	% Number of points over the chord [ ]

    TAU         = str2double(get(DuctValues(1),'string'));      % Thrust ratio
    CDd         = str2double(get(DuctValues(2),'string'));      % Duct section drag coefficient

    % H           = str2double(get(CavValues(1),'string'));       % Shaft centerline depth [m]
    % dV          = str2double(get(CavValues(2),'string'));       % Inflow variation [m/s]


    Propeller_flag	= get(FlagValues(1),'value');               % 0 == turbine, 1 == propeller
    Hub_flag  	= get(FlagValues(3),'value');                   % 0 == no hub, 1 == hub
    Duct_flag	= get(FlagValues(4),'value');                   % 0 == no duct, 1 == duct

    Chord_flag	= get(FlagValues(5),'value');                   % ** CHORD OPTIMIZATION FLAG **

    Viscous_flag	= get(FlagValues(6),'value');               % 0 == viscous forces off (CD = 0), 1 == viscous forces on
    Plot_flag       = get(FlagValues(7),'value');               % 0 == do not display plots, 1 == display plots

    Make2Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 2D plot of the results, 1 == make plot
    Make3Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 3D plot of the results, 1 == make plot

    Analyze_flag	= get(FlagValues(9),'value'); 

    % Cav_flag	= get(FlagValues(10),'value');                   % 0 == do not run cavitation mapping, 1 == run cavitation mapping



    Meanline_index	= get(FoilValues(1),'value');
    Meanline        = char(Meanline_cell(Meanline_index));       	% Meanline form

    Thickness_index	= get(FoilValues(2),'value');
    Thickness       = char(Thickness_cell(Thickness_index));    	% Thickness form


    XR        = str2double(get(XR_in,'string'));                % radius / propeller radius
    XCoD      = str2double(get(XCoD_in,'string'));              % chord / diameter
    XCD       = str2double(get(XCD_in,'string'));               % section drag coefficient

    ri        = str2double(get(ri_in, 'string'));
    VAI       = str2double(get(VAI_in,'string'));               % axial      inflow velocity / ship velocity
    VTI       = str2double(get(VTI_in,'string'));               % tangential inflow velocity / ship velocity

    % f0oc0     = str2double(get(f0oc_in,'String'));              % max section camber    / chord
    Xt0oD     = str2double(get(Xt0oD_in,'String'));             % max section thickness / chord
    skew0     = str2double(get(skew0_in,'String'));             % skew
    rake0     = str2double(get(rake0_in,'String'));             % rake


    save([filename,'_GUIsd'],'Z','N','D','THRUST','Vs','Dhub','rho','Mp','Np',...
        'TAU','CDd','Propeller_flag','Hub_flag','Duct_flag',...
        'Chord_flag','Viscous_flag','Plot_flag','Make2Dplot_flag','Make3Dplot_flag','Analyze_flag',...
        'Meanline','Meanline_index','Thickness',...
        'Thickness_index','filename','XR','XCoD','XCD','ri','VAI','VTI','Xt0oD','skew0',...
        'rake0','XCoD_values','XCLmax_values'); % 'Cav_flag','H','dV',

end
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% ================================ loaddata ===============================
function loaddata(hObject,ED)
    % 
    % This subfunction loads all values presented in the OpenProp GUI.
    %

    global OpenPropDirectory SpecificationsValues DuctValues FlagValues FoilValues Filename...
           XR_in XCoD_in XCD_in VAI_in VTI_in ri_in Xt0oD_in skew0_in rake0_in; % CavValues

    %%

    % ------------------------------------------------------------------------- 
    % Figure out what the current working directory is, 
    % and change directories to /OpenPropDirectory/
    rest = pwd;

    while ~isempty(rest)
        [CurrentDirectory,rest] = strtok(rest,'/');
    end

    if     strcmp(CurrentDirectory,OpenPropDirectory)    % OpenPropDirectory == 'OpenProp_vX.X.X';

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


    % ------------------------------------------------------------------------
    % ------------------------------------------------------------------------

    set(SpecificationsValues(1),'string',num2str(Z));         % number of blades
    set(SpecificationsValues(2),'string',num2str(N));         % propeller speed [RPM]
    set(SpecificationsValues(3),'string',num2str(D));         % propeller diameter [m]
    set(SpecificationsValues(4),'string',num2str(THRUST)); 	% required thrust [N]
    set(SpecificationsValues(5),'string',num2str(Vs));        % ship velocity [m/s]
    set(SpecificationsValues(6),'string',num2str(Dhub));      % hub diameter [m]
    set(SpecificationsValues(7),'string',num2str(rho));       % water density [kg/m^3]
    set(SpecificationsValues(8),'string',num2str(Mp));        % number of vortex panels over the radius
    set(SpecificationsValues(9),'string',num2str(Np));        % Number of points over the chord [ ]

    set(DuctValues(1),'string',num2str(TAU));           % Thrust ratio
    set(DuctValues(2),'string',num2str(CDd));           % Duct section drag coefficient

    % set(CavValues(1),'string',num2str(H));              % Shaft centerline depth [m]
    % set(CavValues(2),'string',num2str(dV));             % Inflow variation [m/s]

    set(FlagValues(1),'value',Propeller_flag);       	% 0 == turbine, 1 == propeller
    set(FlagValues(3),'value',Hub_flag);              	% 0 == no hub, 1 == hub
    set(FlagValues(4),'value',Duct_flag);              	% 0 == no duct, 1 == duct


    set(FlagValues(5),'value',Chord_flag);              % ** CHORD OPTIMIZATION FLAG **
    changeChord;       % Testing whether running the callback function on loading updates fields and avoids loading issue


    set(FlagValues(6),'value',Viscous_flag);          	% 0 == viscous forces off (CD = 0), 1 == viscous forces on
    set(FlagValues(7),'value',Plot_flag);         % 0 == do not display plots, 1 == display plots
    set(FlagValues(8),'value',Make3Dplot_flag);        	% 0 == do not display plots, 1 == display plots

    set(FlagValues(9),'value',Analyze_flag);

    % set(FlagValues(10),'value',Cav_flag);               	% 0 == do not run cavitation mapping, 1 == run cavitation mapping


    set(FoilValues(1),'value',Meanline_index);          	% Meanline form
    set(FoilValues(2),'value',Thickness_index);        	% Thickness form

    set(Filename,'string',filename);                	% Filename prefix

    % ----------------------------------------------
    % --- Loop to set new values for input table ---

    for index = 1 : length(XR);

        set(XR_in(index),'string',num2str(XR(index)));                    % radius / propeller radius
        set(XCoD_in(index),'string',num2str(XCoD(index)));                % chord / diameter
        set(XCD_in(index),'string',num2str(XCD(index)));                  % section drag coefficient

        % set(f0oc_in(index),'String',num2str(f0oc0(index)));             % max section camber    / chord
        set(Xt0oD_in(index),'String',num2str(Xt0oD(index)));              % max section thickness / chord
        set(skew0_in(index),'String',num2str(skew0(index)));              % skew
        set(rake0_in(index),'String',num2str(rake0(index)));              % rake

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
    % ----------------------------------------------

end
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
% function changeCav(hObject,ED)
% 
%     global CavValues FlagValues;
% 
%     if get(FlagValues(10),'value')
% 
%         set(CavValues,'enable','on');
% 
%     else
% 
%         set(CavValues,'enable','off');
% 
%     end
% 
% end
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function changeDir(hObject,ED)

    global Filename OpenPropDirectory;

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
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function changeDuct(hObject,ED)

    global DuctValues FlagValues;

    if get(FlagValues(4),'value')

        set(DuctValues,'enable','on');
    else
        set(DuctValues,'enable','off');
    end
end
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function checkBlades(hObject,ED)

    global SpecificationsValues;

    Z = floor(str2double(get(SpecificationsValues(1),'string')));

    set(SpecificationsValues(1),'string',num2str(Z));

end
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function checkTurbine(hObject,ED)

    global FlagValues SpecificationsValues;

    if get(FlagValues(1),'value')
        set(SpecificationsValues(4),'enable','on')
    else
        % Grays out THRUST field if turbine is selected
        set(SpecificationsValues(4),'string','','enable','off')
    end


end
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function changeChord(hObject,ED)
    % global SpecificationsValues DuctValues CavValues FlagValues FoilValues Filename...
    %        XR_in XCoD_in XCD_in VAI_in VTI_in Xt0oD_in skew0_in rake0_in...
    %        Meanline_cell Thickness_cell;

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
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function changeViscous(hObject,ED)

    global FlagValues XCD_in XCD_values;

    if get(FlagValues(6),'value')

        set(XCD_in,'enable','on');

        for index = 1 : length(XCD_values)

            set(XCD_in(index),'string',num2str(XCD_values(index)));

        end

    else

        XCD_values  = str2double(get(XCD_in,'string'));

        set(XCD_in,'enable','off','string','0');

    end

end
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function updateValues(hObject,ED)

%%

% Update values displayed in the Calculator panel

global SpecificationsValues Values;

N           = str2double(get(SpecificationsValues(2),'string'));  % propeller speed [RPM]
D           = str2double(get(SpecificationsValues(3),'string'));	% propeller diameter [m]
THRUST      = str2double(get(SpecificationsValues(4),'string')); 	% required thrust [N]
Vs          = str2double(get(SpecificationsValues(5),'string'));  % ship velocity [m/s]
Dhub        = str2double(get(SpecificationsValues(6),'string'));  % hub diameter [m]
rho         = str2double(get(SpecificationsValues(7),'string')); 	% water density [kg/m^3]

n       = N/60;    % ** propeller speed [rev/s] = Vs/(Js*D) = N/60
lambda 	= n*2*pi*(D/2)/Vs;
Js      = Vs/(n*D);              	% ** Js = Vs/(n*D) ,  advance coefficient
KT      = THRUST/(rho*n^2*D^4);          % ** KT = THRUST/(rho*n^2*D^4)

CT      = THRUST/(0.5*rho*Vs^2*pi*(D/2)^2);

set(Values(1),'string',num2str(Js));
set(Values(2),'string',num2str(lambda));
set(Values(3),'string',num2str(KT));
set(Values(4),'string',num2str(CT));

end
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function newPlots

    global PlotPanels Toggle OnDesignValues systemToggle ConversionValues UnitsText;

    % ---------------------------------------------------------------------
    % -------------------------- GUI Layout Constants -------------------------
    %
    % This section presents the constant values for the dimensions of those
    % elements used in the GUI, as well as the construction formulas for the
    % different panels, including margins. The formulas are presented in a
    % linear reading order based on the actual order in the GUI.

    titlefontsize   = 25;
    panelfontsize   = 11;

    titleht         = 3;
    textht          = 1.5;
    editboxht       = 2;
    pushht          = 2; % 1.75;

    textbox         = 20;
    ondesignbox     = 15;
    editbox         = 10;

    Windowht        = 44; % 40;           % 1 + Ductboxht + Specificationsboxht + 1 + titleht + 1;
    Window          = 142;          % 1 + Specificationsbox + Inputbox + Flagbox + 1;

    Togglebox       = textbox + 2;
    Toggleboxht     = Windowht - 1 - titleht - 1;

    Displaybox      = Window - 1 - Togglebox - 1;
    Displayboxht    = Windowht - 1 - titleht - 1;

    Conversionbox   = 2 + textbox + ondesignbox + editbox + 2;
    Conversionboxht	= pushht * 8 * 1.25 + 1;

    % =========================================================================
    % ---------------------- Create figure for Plots GUI ----------------------

    Plots       = figure('units','characters','position',[25 20 Window Windowht],...
                         'name','Plots','numbertitle','off','toolbar','figure',... 'menubar','none'
                         'color',[0.702 0.702 0.702],'resize','off');
    % --- --- ---

    Title       = uicontrol(Plots,'style','text','fontsize',titlefontsize,...
                            'fontweight','bold','units','characters','position',...
                            [0 Windowht-1-titleht Window 3],'string',{'OpenProp v3.3.4'});

    % =========================================================================
    % --------- Setup panels --------

    Togglebar           = uibuttongroup('parent',Plots,'title','Figures and plots',...
                                        'fontsize',panelfontsize,'fontweight','bold',...
                                        'units','characters','position',[1 1 Togglebox...
                                        Toggleboxht],'clipping','on','SelectionChangeFcn',...
                                        @togglefn);

    % ---
    % 
    % The panels defined from here on are set up using a (chronological) order
    % criteria based on the different design stages, as presented in the following
    % list:
    % 
    % Inputs:
    % 1 Expanded Blade
    % 2 Thickness profile
    % 3 Inflow profile
    % 
    % Parametric:
    % 4 Efficiency vs Diameter
    % 5 Efficiency vs Rotational Speed (N)
    % 
    % 
    % Design results:
    % 
    % 6 On Design Perf. (default)
    % 
    % 7 Circulation distribution
    % 8 Induced Velocity
    % 9 Beta & Beta i
    % 
    % 10 Expanded Blade (after)
    % 11 t0/D vs Radius at control points
    % 12 CL vs Radius at control points
    % 
    % 13 2D Geometry
    % 14 3D Geometry
    % 
    % 15 Performance curves
    % 
    % ---


    PlotPanels(1)       = uipanel('parent',Plots,'title','Expanded Blade (input blade)',  	'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(2)       = uipanel('parent',Plots,'title','Blade Thickness (input blade)', 	'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(3)       = uipanel('parent',Plots,'title','Inflow Profile',                	'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(4)       = uipanel('parent',Plots,'title','Efficiency vs Diameter',          'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(5)       = uipanel('parent',Plots,'title','Efficiency vs Rotation Speed',    'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(6)       = uipanel('parent',Plots,'title','On-design Performance',           'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','on');

    PlotPanels(7)       = uipanel('parent',Plots,'title','Circulation Distribution',        'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(8)       = uipanel('parent',Plots,'title','Induced Velocity',                'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(9)       = uipanel('parent',Plots,'title','Inflow Angle',                    'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(10)       = uipanel('parent',Plots,'title','Expanded Blade (as designed)',   'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(11)       = uipanel('parent',Plots,'title','Blade Thickness (as designed)',  'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(12)       = uipanel('parent',Plots,'title','Lift Coefficient',               'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(13)       = uipanel('parent',Plots,'title','2D Geometry',                    'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(14)       = uipanel('parent',Plots,'title','3D Geometry',                    'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    PlotPanels(15)       = uipanel('parent',Plots,'title','Performance Curves',             'fontsize',panelfontsize,'fontweight','bold','units','characters','position',[1+Togglebox 1 Displaybox Displayboxht],'clipping','on','visible','off');

    % ---------------------------------------------------------------------
    
    % =========================================================================
    % ------------------------ Togglebar Panel Elements -----------------------

    ToggleText(1)	= uicontrol(Togglebar,'units','characters','style',...
                                'text','string','From Inputs:','position',[1 Toggleboxht-2-pushht textbox textht]);

    Toggle(1)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Expanded Blade','position',[1 Toggleboxht-2-pushht*2 textbox pushht]);

    Toggle(2)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Blade Thickness','position',[1 Toggleboxht-2-pushht*3 textbox pushht]);

    Toggle(3)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Inflow Profile','position',[1 Toggleboxht-2-pushht*4 textbox pushht]);

    ToggleText(2)	= uicontrol(Togglebar,'units','characters','style',...
                                'text','string','Parametric:','position',[1 Toggleboxht-2-pushht*5 textbox textht]);

    Toggle(4)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Effy vs D','position',[1 Toggleboxht-2-pushht*6 textbox pushht]);

    Toggle(5)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Effy vs N','position',[1 Toggleboxht-2-pushht*7 textbox pushht]);

    ToggleText(3)	= uicontrol(Togglebar,'units','characters','style',...
                                'text','string','Design results:','position',[1 Toggleboxht-2-pushht*8 textbox textht]);

    Toggle(6)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Design Performance','position',[1 Toggleboxht-2-pushht*9 textbox pushht],'value',1);

    Toggle(7)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Circulation Distribution','position',[1 Toggleboxht-2-pushht*10 textbox pushht]);

    Toggle(8)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Induced Velocity','position',[1 Toggleboxht-2-pushht*11 textbox pushht]);

    Toggle(9)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Inflow Angle','position',[1 Toggleboxht-2-pushht*12 textbox pushht]);

    Toggle(10)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Expanded Blade','position',[1 Toggleboxht-2-pushht*13 textbox pushht]);

    Toggle(11)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Blade Thickness','position',[1 Toggleboxht-2-pushht*14 textbox pushht]);

    Toggle(12)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Lift Coefficient','position',[1 Toggleboxht-2-pushht*15 textbox pushht]);

    Toggle(13)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','2D Geometry','position',[1 Toggleboxht-2-pushht*16 textbox pushht]);

    Toggle(14)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','3D Geometry','position',[1 Toggleboxht-2-pushht*17 textbox pushht]);

    Toggle(15)       = uicontrol(Togglebar,'units','characters','style',...
                                'togglebutton','string','Performance Curves','position',[1 Toggleboxht-2-pushht*18 textbox pushht]);

    % ---------------------------------------------------------------------
    
    
    % =========================================================================
    % ------------------ On Design Performance Panel Elements -----------------

    ODtextsrc           = {'J =' 'KT =' 'KQ =' 'EFFY =' 'ADEFFY =' 'CT =' 'CQ =' 'CP ='};

    for index = 1 : length(ODtextsrc)

        OnDesignText(index)     = uicontrol(PlotPanels(6),'units','characters','style',...
                                'text','string',ODtextsrc(index),'position',...
                                [2 Displayboxht-4-pushht*index*1.25 ondesignbox textht],...
                                'horizontalalignment','left');

        OnDesignValues(index)	= uicontrol(PlotPanels(6),'units','characters','style',...
                                'edit','string','','position',...
                                [2+ondesignbox Displayboxht-4-pushht*index*1.25 editbox editboxht],...
                                'enable','off');

    end


    % Create a subpanel to group those values that may be converted between SI and English

    ConversionPanel     = uibuttongroup(PlotPanels(6),'units','characters','position',...
                                        [2+ondesignbox+editbox+6 Displayboxht-4-1-pushht*8*1.25...
                                        Conversionbox Conversionboxht],'clipping','on',...
                                        'SelectionChangeFcn',@convertfn);

    systemToggle(1)     = uicontrol(ConversionPanel,'style','radiobutton','units','characters',...
                                    'position',[2 Conversionboxht-2.5 editbox textht],...
                                    'string','SI','horizontalalignment','left','enable','off');

    systemToggle(2)     = uicontrol(ConversionPanel,'style','radiobutton','units','characters',...
                                    'position',[2+editbox Conversionboxht-2.5 editbox textht],...
                                    'string','English','horizontalalignment','left','enable','off');


    Conversiontextsrc   = {'Ship speed (Vs) =' 'Rotation speed (N) =' 'Diameter (D) ='...
                           'Thrust (T) =' 'Torque (Q) =' 'Power (P) ='};

    unitssrc            = {'m/s' 'RPM' 'm' 'N' 'Nm' 'W'};

    for index = 1 : length(Conversiontextsrc)

        ConversionText(index)	= uicontrol(ConversionPanel,'units','characters','style',...
                                'text','string',Conversiontextsrc(index),'position',...
                                [2 Conversionboxht-2.5-pushht*(index)*1.25 textbox textht],...
                                'horizontalalignment','left');

        ConversionValues(index)	= uicontrol(ConversionPanel,'units','characters','style',...
                                'edit','string','','position',...
                                [2+textbox Conversionboxht-2.5-pushht*(index)*1.25 ondesignbox editboxht],...
                                'enable','off');

        UnitsText(index)	= uicontrol(ConversionPanel,'units','characters','style',...
                                'text','string',unitssrc(index),'position',...
                                [2+textbox+ondesignbox+1 Conversionboxht-2.5-pushht*(index)*1.25 editbox textht],...
                                'horizontalalignment','left');



    end


end
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function togglefn(hObject,ED)

    global PlotPanels Toggle;

    for index = 1 : length(Toggle)

        if get(Toggle(index),'value')
            set(PlotPanels(index),'visible','on');
        else
            set(PlotPanels(index),'visible','off');
        end
    end

end
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function convertfn(hObject,ED)

    global systemToggle ConversionValues UnitsText;

    global pt;

    if get(systemToggle(1),'value')

        set(ConversionValues(1),'string',num2str(pt.input.Vs));
        set(ConversionValues(2),'string',num2str(pt.input.N));
        set(ConversionValues(3),'string',num2str(pt.input.D));
        set(ConversionValues(4),'string',num2str(pt.input.THRUST));
        set(ConversionValues(5),'string',num2str(pt.design.Q));
        set(ConversionValues(6),'string',num2str(pt.design.P));

        unitssrc            = {'m/s' 'RPM' 'm' 'N' 'Nm' 'W'};

        for index = 1 : length(unitssrc)

            set(UnitsText(index),'string',unitssrc(index));

        end

    else

        set(ConversionValues(1),'string',num2str(pt.input.Vs*1.94384449));     % Convert to knots
        set(ConversionValues(2),'string',num2str(pt.input.N));
        set(ConversionValues(3),'string',num2str(pt.input.D*3.2808399));       % Convert to feet
        set(ConversionValues(4),'string',num2str(pt.input.THRUST*0.224808943));     % Convert to lbf
        set(ConversionValues(5),'string',num2str(pt.design.Q*0.737562149277));  % Convert to lb ft
        set(ConversionValues(6),'string',num2str(pt.design.P*0.00134102209));   % Convert to Hp

        unitssrc            = {'knots' 'RPM' 'ft' 'lb' 'lb/ft' 'HP'};

        for index = 1 : length(unitssrc)

            set(UnitsText(index),'string',unitssrc(index));

        end

    end

end
% =========================================================================
% =========================================================================
% =========================================================================
%%
% =========================================================================
% =========================================================================
% =========================================================================
function Selectfn(hObject,ED)

    global Select;

    global OpenPropDirectory SpecificationsValues DuctValues FlagValues FoilValues Filename...
        XR_in XCoD_in XCD_in VAI_in VTI_in ri_in Xt0oD_in skew0_in rake0_in...
        Meanline_cell Thickness_cell XCoD_values XCLmax_values; % CavValues 

    if get(Select,'value')==1
        % Do nothing - already in Single Design
    elseif get(Select,'value')==2

        % =========================================================================
        % ================================ savedata ===============================
        %
        % This subfunction saves all values presented in the OpenProp GUI.
        %

        filename   	= get(Filename,'string');                       % Filename prefix

        Z           = str2double(get(SpecificationsValues(1),'string'));  % number of blades
        N           = str2double(get(SpecificationsValues(2),'string'));  % propeller speed [RPM]
        D           = str2double(get(SpecificationsValues(3),'string'));	% propeller diameter [m]
        THRUST      = str2double(get(SpecificationsValues(4),'string')); 	% required thrust [N]
        Vs          = str2double(get(SpecificationsValues(5),'string'));  % ship velocity [m/s]
        Dhub        = str2double(get(SpecificationsValues(6),'string'));  % hub diameter [m]
        rho         = str2double(get(SpecificationsValues(7),'string')); 	% water density [kg/m^3]
        Mp          = str2double(get(SpecificationsValues(8),'string')); 	% number of vortex panels over the radius
        Np          = str2double(get(SpecificationsValues(9),'string')); 	% Number of points over the chord [ ]

        TAU         = str2double(get(DuctValues(1),'string'));      % Thrust ratio
        CDd         = str2double(get(DuctValues(2),'string'));      % Duct section drag coefficient

         Propeller_flag = get(FlagValues(1),'value');               % 0 == turbine, 1 == propeller
               Hub_flag = get(FlagValues(3),'value');                   % 0 == no hub, 1 == hub
              Duct_flag = get(FlagValues(4),'value');                   % 0 == no duct, 1 == duct  
             Chord_flag	= get(FlagValues(5),'value');                   % ** CHORD OPTIMIZATION FLAG **
           Viscous_flag = get(FlagValues(6),'value');               % 0 == viscous forces off (CD = 0), 1 == viscous forces on
              Plot_flag = get(FlagValues(7),'value');               % 0 == do not display plots, 1 == display plots

        Make2Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 2D plot of the results, 1 == make plot
        Make3Dplot_flag = get(FlagValues(8),'value');               % 0 == do not make a 3D plot of the results, 1 == make plot
           Analyze_flag = get(FlagValues(9),'value'); 


        Meanline_index  	= get(FoilValues(1),'value');
        Meanline        = char(Meanline_cell(Meanline_index));       	% Meanline form

        Thickness_index	= get(FoilValues(2),'value');
        Thickness       = char(Thickness_cell(Thickness_index));    	% Thickness form

        XR        = str2double(get(XR_in,'string'));                % radius / propeller radius
        XCoD      = str2double(get(XCoD_in,'string'));              % chord / diameter
        XCD       = str2double(get(XCD_in,'string'));               % section drag coefficient

        ri        = str2double(get(ri_in, 'string'));
        VAI       = str2double(get(VAI_in,'string'));               % axial      inflow velocity / ship velocity
        VTI       = str2double(get(VTI_in,'string'));               % tangential inflow velocity / ship velocity

        Xt0oD     = str2double(get(Xt0oD_in,'String'));             % max section thickness / chord
        skew0     = str2double(get(skew0_in,'String'));             % skew
        rake0     = str2double(get(rake0_in,'String'));             % rake

        save('OpenPropTempFile0307122010','Z','N','D','THRUST','Vs','Dhub','rho','Mp','Np',...
            'TAU','CDd','Propeller_flag','Hub_flag','Duct_flag',...
            'Chord_flag','Viscous_flag','Plot_flag','Make2Dplot_flag','Make3Dplot_flag','Analyze_flag',...
            'Meanline','Meanline_index','Thickness',...
            'Thickness_index','filename','XR','XCoD','XCD','ri','VAI','VTI','Xt0oD','skew0',...
            'rake0','XCoD_values','XCLmax_values'); % 'Cav_flag','H','dV',




        OpenPropParam;

    elseif get(Select,'value')==3
        OpenPropAnalyze;
    end

end
% =========================================================================
% =========================================================================
% =========================================================================


