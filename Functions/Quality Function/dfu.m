% Uniform coverage quality function derivative with respect to z
% Return a vector with the coverage quality for each altitude in z
function f = dfu(z, zmin, zmax)
    
    f = zeros(size(z));
    
    for i=1:length(z)
        if z(i) < zmin || z(i) > zmax
            f(i) = 0;
        else
            f(i) = (4*(z - zmin)*((z - zmin)^2 - (zmax - zmin)^2))/(zmax - zmin)^4;
        end
    end
