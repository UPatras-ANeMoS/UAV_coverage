% Copyright 2016 Sotiris Papatheodorou
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

% offsets point x a distance r perpendicular to line
% the sign of r defines the direction of the offset
% line = [linex ; liney]
function x_new = point_perp(line, x, r)

th = atan2( (line(2,2)-line(2,1)) , (line(1,2)-line(1,1)) );

if th >= 0
    x_new(1,1) = x(1) + r*sin(th);
    x_new(2,1) = x(2) - r*cos(th);   
else   
    x_new(1,1) = x(1) + r*cos(pi/2 - th);
    x_new(2,1) = x(2) - r*sin(pi/2 - th);
end

