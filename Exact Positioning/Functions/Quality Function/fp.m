% Copyright 2017 Sotiris Papatheodorou
% 
% Licensed under the Apache License, Version 2.0 (the \"License\");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%    http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an \"AS IS\" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

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

