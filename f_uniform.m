function c = f_uniform(zi, zmin, zmax)
    
    c = zeros(size(zi));
    
    for i=1:length(zi)
        if zi(i) < zmin
            c(i) = 1;
        elseif zi(i)> zmax
            c(i) = 0;
        else
            c(i) = ( (zi(i)-zmin)^2 - (zmax-zmin)^2 )^2 / (zmax-zmin)^4;
        end
    end
