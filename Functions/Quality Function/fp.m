function f = fp(x, y, xi, yi, zi, zmin, zmax, a, b )
    fi = fu(zi, zmin, zmax);
    Ki = (1-b) / (tan(a)*zi)^2;
    
    % Zero if it is outside the sensing circle
    if ((x-xi)^2 + (y-yi)^2) <= (zi*tan(a))^2
        % Paraboloid
        f = (1 - Ki*( (x-xi)^2 + (y-yi)^2 ) ) * fi;
        % Same as fi(z)
    else
        f = 0;
    end