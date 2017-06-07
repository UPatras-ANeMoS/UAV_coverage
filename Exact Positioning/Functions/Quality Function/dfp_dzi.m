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
