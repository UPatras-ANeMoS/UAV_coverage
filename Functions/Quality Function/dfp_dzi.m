% Decreasing coverage quality function (inverted paraboloid) derivative
% with respect to zi
% return the coverage quality derivative of point [x y] with respect to node 
% [xi yi zi]
function f = dfp_dzi(x, y, xi, yi, zi, zmin, zmax, a, b)
    % Sensing disk radius
    r = zi*tan(a);
    
    % Check if the point is inside the sensing circle
    if norm([x ; y] - [xi ; yi]) <= r
        % Calculate uniform quality and its derivative
        fiu = fu(zi, zmin, zmax);
        dfiu = dfu(zi, zmin, zmax);
        
        f = (((x-xi)^2 + (y-yi)^2)*(b-1)/(zi*tan(a))^2 + 1) * dfiu + ...
            2*((x-xi)^2 + (y-yi)^2)*(1-b) / (zi^3*tan(a)^2) * fiu;
    else
        f = 0;
    end
