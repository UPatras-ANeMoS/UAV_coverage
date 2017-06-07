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

% Returns the Guaranteed Voronoi edges between two disks
function [E, a, c] = GV_edge(A , B , uncertA, uncertB , p , max_line)
% Voronoi edges are hyperbolas (polygons)
% Returns the Voronoi edges between points A and B
% p defines how many points the polygon will have (MAKE SURE ITS ODD)
% l is the largest x on the polygon

% Calculate parameters
c = eucl_dist(A,B) / 2;
a = (uncertA + uncertB) / 2;
b = sqrt(c^2 - a^2);


% if c <= a there are no cells
if c <= a
    E = [];
    return
end

% if a = 0 then the edge is the classic Voronoi edge, the perpendicular bisector of A and B
if a == 0
	x = [0 0 max_line max_line];
	y = [max_line -max_line -max_line max_line];
else
    % if a > 0 then the edge is a hyperbola branch
    % Use the parametric equation of the hyperbola
    t = linspace( -acosh(max_line/a), acosh(max_line/a), p);
%     t = logspace( -10, ceil(log(acosh(max_line/a))), (p-1)/2 );
%     t = [-fliplr(t) 0 t];
    x = a .* cosh(t);
    y = b .* sinh(t);
end

% Edge of point in right hand plane
EB = [ x;
        y];
% Edge of point in left hand plane
EA = [ -x;
        y];

% EA: edge of point A wrt B
% EB: edge of point B wrt A
E = [EA; EB];


% Rotation
theta = atan2( B(2) - A(2) , B(1) - A(1) );
for i = 1:length( E(1,:) ) % p points per edge
	E(1:2,i) = rot( E(1:2,i) , theta);
	E(3:4,i) = rot( E(3:4,i) , theta);
end


% Translation vector
f = (A+B) / 2;

% Translation
E = [   (E(1,:) + f(1)) ;
        (E(2,:) + f(2)) ;
        (E(3,:) + f(1)) ;
        (E(4,:) + f(2)) ];
