% Copyright 2017 Sotiris Papatheodorou
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%    http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% Given a point on the hyperbola with known foci and semi-major axis,
% the value of the parameter t of the hyperbola at that point is returned.
% If the given point is not on the hyperbola, NaN is returned.
function t = hyperbola_point2parameter(x, q1, q2, a)

% Get hyperbola parameters
c = norm(q1-q2)/2;
b = sqrt(c^2-a^2);

% Get hyperbola center and translate point
center = (q1+q2)/2;
xx = x - center;

% Get hyperbola rotation angle and rotate point
theta = atan2( q1(2)-q2(2), q1(1)-q2(1) );
xx = rot(xx, -theta);


% Error tolerance for point on hyperbola test
e = 10^-10;
if (abs(xx(1)^2/a^2 - xx(2)^2/b^2 - 1) <= e) && isreal(b)
    t = asinh( xx(2)/b );
else
    % The point is not on the hyperbola or the hyperbola is degenerate
    t = NaN;
end
