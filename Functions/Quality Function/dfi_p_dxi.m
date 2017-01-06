function c = dfi_p_dxi(x, y, xi, yi, zi, zmin, zmax, a, b )
    % for dfi_p_dyi just swap x with y and xi with yi
    rmax = zi*tan(a);
    Ki = (1-b) / (tan(a)*zi)^2;
    
    if ((x-xi)^2 + (y-yi)^2) <= (rmax)^2
        c = -((2*x - 2*xi)*(b - 1)*((zi - zmin)^2 - ...
            (zmax - zmin)^2)^2)/(zi^2*tan(a)^2*(zmax - zmin)^4);
    else
        c = 0;
    end