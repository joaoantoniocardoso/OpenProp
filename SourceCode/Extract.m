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


% --------------------------------------------------------------- Extract.m
% Created: 7/7/2011, Brenden Epps
%
% This function extracts field 'x' from a size (M,N) data structure D
% and assigns it to X, such that
%
%   X(m,n,p,q)) = D(m,n).x(p,q)
%
% Inputs:
%   D   size (M,N) data structure with field 'x'
%   x   string, where field D.x is size (P,Q) 
%
% Outputs:
%   X   [M,N,P,Q] matrix of values of x
%
%   Note: 
%       This implementation can handle the case where x is 
%       buried in sub-structures: 
%                                     X(m,n) = D(m,n).struct1.struct2.x
%
%       and so on, where there can be as many sub-structures as desired...,
%       To do so, execute:
%                           X = Extract(D,'struct1.struct2.x')
%
%   Note:
%       To handle the case:     X(m,n) = D.struct1.struct2(m,n).x
%       call:
%              X = Extract(D.struct1.struct2,'x')
%
% -------------------------------------------------------------------------

function [X] = Extract(D,x)

[M,N] = size(D);


% If x is buried in a sub-structure, need to extract sub-structure first
Done_flag = 0;

while Done_flag == 0

    [struct, x] = strtok(x,'.');

    if ~isempty(x)

        x = x(2:end);

        for m = 1:M
            for n = 1:N

                temp(m,n) = getfield(D(m,n),struct);

            end
        end

        D = temp; 

        clear temp

    else
        Done_flag = 1;

        x = struct;
    end
end


% Extract field x from data structure
[P,Q] = size( getfield(D(1,1),x) );

X = zeros(M,N,P,Q);


for m = 1:M
    for n = 1:N    
        X(m,n,:,:) = getfield(D(m,n),x);        
    end
end




