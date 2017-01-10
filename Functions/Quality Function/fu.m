% Uniform coverage quality function
% Return a vector with the coverage quality for each altitude in z
function f = fu(z, zmin, zmax)
    
    f = zeros(size(z));
    
    for i=1:length(z)
        if z(i) < zmin
            f(i) = 1;
        elseif z(i) > zmax
            f(i) = 0;
        else
            f(i) = ( (z(i)-zmin)^2 - (zmax-zmin)^2 )^2 / (zmax-zmin)^4;
        end
    end
