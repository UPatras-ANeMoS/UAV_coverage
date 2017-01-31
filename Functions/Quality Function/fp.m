% Decreasing coverage quality function (inverted paraboloid)
% return the coverage quality of point [x y] with respect to node 
% [xi yi zi]
function f = fp(x, y, xi, yi, zi, zmin, zmax, a, b )

% Sensing disk radius
r = zi*tan(a);

% Initialize result
f = zeros(size(x));

% DOES NOT WORK IF size(x) ~= size(y)
s = size(x);
s = s(1)*s(2);

for k=1:s
    % Check if the point is inside the sensing circle (MAYBE REMOVE THIS?)
    if norm([x(k) ; y(k)] - [xi ; yi]) <= r
        % Calculate uniform quality for paraboloid maximum value
        fiu = fu(zi, zmin, zmax);
        
        % Value of paraboloid at point [x y]
        f(k) = (1 - (1-b) * ((x(k)-xi)^2 + (y(k)-yi)^2) / r^2) * fiu;
    end
end

