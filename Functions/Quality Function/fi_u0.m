function c = fi_u0(x, y, xi, yi, zi, zmin, zmax, a)
    
    if ((x-xi)^2 + (y-yi)^2) <= (zi*tan(a))^2
        if zi < zmin
            c = 1;
        elseif zi> zmax
            c = 0;
        else
            c = ( (zi-zmin)^2 - (zmax-zmin)^2 )^2 / (zmax-zmin)^4;
        end
    else
        c = 0;
    end
