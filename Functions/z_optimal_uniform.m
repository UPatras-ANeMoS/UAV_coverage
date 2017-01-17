% Returns the optimal altitude for the given network parameters
% and uniform coverage quality (third argument is a placeholder for a)
function zopt = z_optimal_uniform(zmin, zmax)

zopt = (2*zmin)/3 + (3*zmax^2 - 6*zmax*zmin + 4*zmin^2)^(1/2)/3;
