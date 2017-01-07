% Returns the optimal altitude for the given network parameters
% and decreasing coverage quality
function zopt = z_optimal_decreasing(zmin, zmax, a, b)

zopt = (2*zmin)/3 + (3*zmax^2 - 6*zmax*zmin + 4*zmin^2)^(1/2)/3;
