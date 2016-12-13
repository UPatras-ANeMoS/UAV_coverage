function c = fi_u(zi, zmin, zmax)
    
    if zi < zmin
        c = 1;
    elseif zi> zmax
        c = 0;
    else
        c = ( (zi-zmin)^2 - (zmax-zmin)^2 )^2 / (zmax-zmin)^4;
    end
