% Decreasing coverage quality function (inverted paraboloid) derivative
% with respect to xi (swap x with y and xi with yi for derivative with
% respect to yi)
% return the coverage quality derivative of point [x y] with respect to node 
% [xi yi zi]
function f = dfp_dxi(x, y, xi, yi, zi, zmin, zmax, a, b )
    % Sensing disk radius
    r = zi*tan(a);
    
    % Check if the point is inside the sensing circle
    if norm([x ; y] - [xi ; yi]) <= r
        % Calculate uniform quality used in derivative
        fiu = fu(zi, zmin, zmax);
        
        f = 2*(1-b) * (x-xi) / r^2* fiu;
    else
        f = 0;
    end
