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
        

% ------------------------------------------------------------ Rainbow Color matrix
CLR = [     1       0       0;      ... % (1) Red
            1       0.5     0;      ... % (2) Orange
            0.8     0.8     0;      ... % (3) Burnt Yellow  
            0       0.8     0;      ... % (4) Green
            0       0       1;      ... % (5) Blue
            0.75    0       0.75;   ... % (6) Purple
            0.3     0.3     0.3];   ... % (7) Gray   

        
% Cycle through the colors
figure, hold on,

for i = 1:size(CLR,1)
    
    ind = mod(i-1,size(CLR,1))+1;
    
    plot([0 1],[i i],'Linewidth',2,'Color',CLR(ind,:)),
end           
    axis([0 1 0 12])    

% ----------------------------------------------------------------------------
%% Plot a rainbow colormap:
%
figure, hold on,
plot([0 1],[1 1]    ,'Linewidth',2,'Color',[1 0 0]),          % Red
plot([0 1],[2 2]    ,'Linewidth',2,'Color',[1 0.5 0]),        % Orange
plot([0 1],[3 3]    ,'Linewidth',2,'Color',[0.95 0.95 0]),    % Yellow
plot([0 1],[4 4]    ,'Linewidth',2,'Color',[0 0.9 0]),        % Green
plot([0 1],[5 5]    ,'Linewidth',2,'Color',[0 0 1]),          % Blue
plot([0 1],[6 6]    ,'Linewidth',2,'Color',[0.75 0 0.75]),    % Purple
plot([0 1],[7 7]    ,'Linewidth',2,'Color',[0.5 0.5 0.5]),    % Gray
plot([0 1],[8 8]    ,'Linewidth',2,'Color',[0 0 0]),          % Black
plot([0 1],[9 9]    ,'Linewidth',2,'Color',[0.75 0.5 0.25]),  % Brown
plot([0 0.5],[10 10],'Linewidth',2,'Color',[0 1 1]),          % Cyan
plot([0.5 1],[10 10],'Linewidth',2,'Color','c'),              % Cyan
plot([0 0.5],[11 11],'Linewidth',2,'Color',[1 0 1]),          % Magenta
plot([0.5 1],[11 11],'Linewidth',2,'Color','m'),              % Magenta
plot([0 1],[12 12]  ,'Linewidth',2,'Color',[0.25 0.5 0.75]),  % Steel blue
plot([0 1],[13 13]  ,'Linewidth',2,'Color',[0.25 0.25 0.75]), % Navy blue
%%
%close all,

% Primary colors:
rcmap = [1    0    0; ...
         1    0.5  0; ...
         0.85  0.85 0; ...
         0    0.5  0;...
         0    0    1;...
         0.75 0    0.75;...
         0    0    0];

% Dark colors:
rcmap = [1    0    0; ...   % Red
         0.25 0.5 0.75;...  % Steel blue
         1    0.5  0; ...   % Orange
         0.75 0.75 0; ...   % Burnt Yellow
         0.25 0.25 0.75;... % Navy blue
         0    0.5  0;...    % Green
         0    0    1;...    % Blue
         0.75 0.5 0.25;...  % Brown
         0.75 0    0.75;... % Purple
         0    0    0];      % Black
    
  
     
% % Interpolate between colors    
% figure, hold on,
%     plot([0 1],[0 0],'k')
% 
%     N = 21;
% for i = 1:N
%     i_cmap = 1 + (i-1)*(size(rcmap,1)-1)/(N-1);
%     clr    = rcmap(floor(i_cmap),:)+(i_cmap-floor(i_cmap))*(rcmap(ceil(i_cmap),:)-rcmap(floor(i_cmap),:));
%     
%     plot([0 1],[i i],'Linewidth',2,'Color',clr),
% end 

% % Round off to the nearest color    
% figure, hold on,
%     plot([0 1],[0 0],'k')
% 
%     N = 21;
% for i = 1:N
%     i_cmap = floor(1 + (i-1)*(size(rcmap,1)-1)/(N-1));
%     clr    = rcmap(i_cmap,:);
%     
%     plot([0 1],[i i],'Linewidth',2,'Color',clr),
% end  


% % Cycle through the colors
% figure, hold on,
%     plot([0 1],[0 0],'k')
% 
% for i = 1:size(rcmap,1)
%     clr    = rcmap(mod(i-1,size(rcmap,1))+1,:);
%     
%     plot([0 1],[i i],'Linewidth',2,'Color',clr),
% end    



% % Cycle through the colors in the 'lines' colormap:
% rcmap = lines;
% 
% figure, hold on,
%     plot([0 1],[0 0],'k')
% 
%     N = 7;
% for i = 1:N
%     clr    = rcmap(i,:);
%     
%     plot([0 1],[i i],'Linewidth',2,'Color',clr),
% end   

    


 
