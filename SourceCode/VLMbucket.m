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


% -------------------------------------------------------------------------
% Last update: 11/2/2011 Brenden Epps
%
% Modified from XBucket.m Code by Chris Peterson and 
% VLMBucket.m by Richard Kimball
% -------------------------------------------------------------------------
%
% Create cavitation bucket diagram for inviscid flow.
%
% -------------------------------------------------------------------------


function [] = VLMbucket()
%%
warning off

panels      = 80;          %Sets number of panels if resetting in XFOIL

% foc_rng     = [0.0093 0.0093];  %Camber ratio range (Must be <0.1 for 4-digit NACA)
foc_rng     = [0.01 0.01];  %Camber ratio range (Must be <0.1 for 4-digit NACA)
foc_step    = 0.01;         %Camber ratio increment
foc_all     = foc_rng(1):foc_step:foc_rng(2);
% 
toc_rng     = [0.02 0.18];    % Thicness ratio range
% toc_rng     = [0.14 0.14];    % Thicness ratio range
toc_step    = 0.04;          % Thickness ratio increment
toc_all     = toc_rng(1):toc_step:toc_rng(2);
 
Alpha_rng   = [-5 8];       %Angle of attack range
Alpha_delta = 0.1;          %Angle of attack increment
% Alpha_rng   = [-0.5 0.3];       %Angle of attack range
% Alpha_delta = 0.01;          %Angle of attack increment
Alpha_all   = Alpha_rng(1):Alpha_delta:Alpha_rng(2);
 
Nt = length(  toc_all);
Na = length(Alpha_all);
Nf = length(  foc_all);

 cpmni = zeros(Nt,Na,Nf);  % Data array for minimum Cp
xcpmni = zeros(Nt,Na,Nf);  % Data array for location of CPmin

% fig(1);
% fig(2);
load colors

%Start of main calculation loops
for k = 1:Nf                                  % Calculate over range of f/c
    
    fo_c = foc_all(k);

    for j = 1:Nt                              % Calculate over range of t/c

        to_c = toc_all(j);

        
        for kk = 1:Na                        % Calculate over range of alfa
            
            alpha = Alpha_all(kk);
            
            % call VLM
            [xt, CPU, CPL] = VLM2D(panels,fo_c/0.0679,alpha*pi/180,to_c);
            
            % extract  xcp and CPmin and generate array
            [CPMINU,indU]=max(CPU);
            [CPMINL,indL]=max(CPL);
            
            if CPMINU > CPMINL
                  cpmni(j,kk,k) = CPMINU;
                xcpmini(j,kk,k) = xt(indU);
            else
                  cpmni(j,kk,k) = CPMINL;
                xcpmini(j,kk,k) = xt(indL);
            end
            
%             ind = mod(kk-1,11)+1;            
%             fig(1); plot(xt,CPU,'-','Color',CLR(ind,:)), plot(xt,CPL,'-','Color',CLR(ind,:)), 
%             fig(2); plot(cpmni(j,kk,k),alpha,'.','Color',CLR(ind,:)),% axis([0 3 -2 2]),  pause(0.1),
        end

    end
end

% Generates Bucket diagrams, new plot for each Fo/C  

handle = fig;
    cmap = colormap(hsv(Nt)); % Generates color distibution
    set(gca,'ColorOrder',cmap);
    
    % title('NACA 65A010, a=0.8, f0/c = 0.02, Lighthill')
load colors

%%
for k = 1:Nf                                  % Calculate over range of f/c
    
    fo_c = foc_all(k);

    for j = 1:Nt                              % Calculate over range of t/c

        to_c = toc_all(j);        
        
  
    
    plot(cpmni(j,:,k), Alpha_all,'LineWidth',2,'Color',cmap(j,:)); %Plots Alpha vs. -Cpmin
    
    xlim([0 3]);
    if fo_c > 0
        ylim(Alpha_rng);                         %Set plot X/Y limits
    else
        ylim([0 8]);
    end
    
    box on,
    xlabel('-CPmin  [  ]'              ,'FontName','Times','FontSize',FontSize_label); 
    ylabel('net angle of attack  [deg]','FontName','Times','FontSize',FontSize_label);
    
    
    title({['INVISCID Brockett Diagram',10, 'Cavitation Bucket', 10,...
        ' Fo/c = ', num2str(fo_c), '   ', '\Delta\alpha = ', num2str(Alpha_delta)]});

    tau = toc_rng(1):toc_step:toc_rng(2);           %Used for legend
    leg_st = cell(1,length(tau));                   %Initializes cells
    for i = 1:length(tau);                          %Set vales to cells
        leg_st(i) = {num2str(tau(i))};
    end
    legend(leg_st, 'Location', 'SouthEast')
    
    %visc_tog=0
    %
    %if visc_tog == 1
    %    figure();
    %    hold on; grid on;
    %    cmap = colormap(hsv(toc_rng(2)/toc_step+1));
%				%Generates color distibution
    %    set(gca,'ColorOrder',cmap);
    %    plot(-cpmnv(:,:,k), Alpha_all(1:length(cpmnv)));        %Plots Alpha vs. -Cpmin
    %    xlim([0 3]);    
    %    ylim(Alpha_rng);                %Set plot X/Y limits
    %    xlabel('-CP_m_i_n'); ylabel('Angle of Attack (\alpha)');
    %    if foil_type == 'NACA'
    %        title_name = [foil_type, ' ', foil_name];
    %        elseif foil_type == 'LOAD'
    %        title_name = ['Meanline: ', load_mean, '.  Thickness: ', load_thck];
    %        else
    %        title_name = 'UNKNOWN TYPE';
    %        
    %    end
    %    title({['VISCOUS Brockett Diagram',10, title_name,10,...
    %        ' Fo/c = ', num2str(fo_c), '    \Delta\alpha = ', num2str(Alpha_delta)]});
    %   tau = toc_rng(1):toc_step:toc_rng(2);           %Used for legend
    %    leg_st = cell(1,length(tau));                   %Initializes cells
    %    for i = 1:length(tau);                          %Set vales to cells
    %        leg_st(i) = {num2str(tau(i))};
    %   end
    %    legend(leg_st, 'Location', 'SouthEast')
    %end
   
    end
end

set(gca,'FontSize',FontSize_axis,'FontName','Times')
%
text(2.3,0.5,'t/c = 0.02','FontSize',FontSize_axis-2,'FontName','Times')
text(1.2,5.2,'t/c = 0.2','FontSize',FontSize_axis-2,'FontName','Times','Rotation',25)

%%

filename = 'fig_bucket_diagram';

% File types: -dpdf, -dtiff, -djpeg, -depsc
% Resolution: -r72, -r144
print(gcf,filename,'-r72','-depsc'),   

