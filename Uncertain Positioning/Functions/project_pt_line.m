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

% Returns the projection of x on line. May also check if the projection
% lies on the line segment defined by the two points in line
function [pr, on] = project_pt_line( x, line )

% Line parameters ax + by + c = 0
a = line(2,2) - line(2,1);
b = line(1,1) - line(1,2);
c = -a * line(1,1) - b * line(2,1);

pr = [ ( b*(b*x(1)-a*x(2))-a*c ) / (a*a+b*b) ;
    ( a*(-b*x(1)+a*x(2))-b*c ) / (a*a+b*b) ];

% Also check if the projection is on the line segment
if nargout == 2
    xmin = min( line(1,:) );
    xmax = max( line(1,:) );
    ymin = min( line(2,:) );
    ymax = max( line(2,:) );
    
    if (pr(1) >= xmin) && (pr(1) <= xmax) && ...
       (pr(2) >= ymin) && (pr(2) <= ymax)
        on = 1;
    else
        on = 0;
    end
end
