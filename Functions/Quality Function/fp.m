% Decreasing coverage quality function (inverted paraboloid)
% return the coverage quality of point [x y] with respect to node 
% [xi yi zi]
function f = fp(x, y, xi, yi, zi, zmin, zmax, a, b )
    % Sensing disk radius
    r = zi*tan(a);
    
    % Check if the point is inside the sensing circle
    if norm([x ; y] - [xi ; yi]) <= r
        % Calculate uniform quality for paraboloid maximum value
        fiu = fu(zi, zmin, zmax);
        
        % Value of paraboloid at point [x y]
        f = (1 - (1-b) * ((x-xi)^2 + (y-yi)^2) / r^2) * fiu;
    else
        f = 0;
		fprintf('out\n')
    end
