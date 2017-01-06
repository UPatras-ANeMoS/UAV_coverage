function c = dfi_p_dzi(x, y, xi, yi, zi, zmin, zmax, a, b )
    rmax = zi*tan(a);
    Ki = (1-b) / (tan(a)*zi)^2;
    
    if ((x-xi)^2 + (y-yi)^2) <= (rmax)^2
        c = (2*((((x - xi)^2 + (y - yi)^2)*(b - 1))/(zi^2*tan(a)^2) + 1)*...
            (2*zi - 2*zmin)*((zi - zmin)^2 - (zmax - zmin)^2))/(zmax - zmin)^4 - ...
            (2*((x - xi)^2 + (y - yi)^2)*(b - 1)*((zi - zmin)^2 - ...
            (zmax - zmin)^2)^2)/(zi^3*tan(a)^2*(zmax - zmin)^4);
    else
        c = 0;
    end