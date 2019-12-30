% -------------------------------------------------------------------------
% 17-Feb-2012 Blade thickness profile suitable for 3D printing
%
% Can change hub thickness "t0hub" to keep t0/c < 0.2 at the hub.
%
% Can change D and XCoD to keep t0/c < 0.2 for rest of blade.
%
% Keep t0tip and t0tpm as is, so 3D printer can resolve tip blade section.
% -------------------------------------------------------------------------
clear, close all, clc,

% -------------------------------------------------------------------------
XR      = [0.2000    0.3000    0.4000    0.5000    0.6000    0.7000    0.8000    0.9000    0.9500    0.9800    0.9900    1.0000];  % radius / propeller radius
XCoD    = [0.1740    0.2280    0.2750    0.3130    0.3380    0.3480    0.3340    0.2810    0.2190    0.1530    0.1150    0.0010];  % chord / diameter  (12/6/2011 NOTE: finite chord at tip: 0.0010)

D       = 0.254;  % m, propeller diameter 
Rhub_oR = 0.2; 

% ---------------------------------------------
% Blade thickness profile
t0hub = 0.400*0.0254; % [m] == 0.400 inch (for model), max thickness at hub section
t0tip = 0.150*0.0254; % [m] == 0.150 inch (for model), max thickness at tip section
t0tpm = 0.080*0.0254; % [m] == 0.080 inch (for model), modified tip thickness
XRmax = 0.80;              % maximum XR for which thickness reduction is less than 1%
TTRF  = t0tpm/t0tip;       % Tip Thickness Reduction Factor == modified thickness at tip / baseline thickness at tip
HTTR  = t0hub/t0tip;       % Hub-Tip Thickness Ratio        == t0(hub) / t0(tip)

t0    = t0tip*(HTTR - (HTTR-1).*(XR-Rhub_oR)/(1-Rhub_oR))  .* (1-(1-TTRF)*exp(-4.6*(1-XR)/(1-XRmax)));
t0oD0 = t0/D;
t0oc0 = t0oD0 ./ XCoD;


fig; plot(XR,t0/0.0254)

fig; plot(XR,XCoD)

fig; plot(XR,t0oc0), axis([0 1 0 0.5])

