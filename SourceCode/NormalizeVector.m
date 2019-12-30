% This function scales the given vector to have a norm of 1.

function [vx vy vz] = NormalizeVector(vx,vy,vz)


veclength = sqrt(vx.^2 + vy.^2 + vz.^2);

vx = vx./veclength;
vy = vy./veclength;
vz = vz./veclength;