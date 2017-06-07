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

function b = isvisible(side, pt, pol)
% checks if side of polygon pol is visible from point pt
% side contains the two indices of the sides vertices

% create the triangle defined by side and pt
tri = [ pt pol(:,side(1)) pol(:,side(2)) ];
% faster without poly2cw
% [tri(1,:), tri(2,:)] = poly2cw(tri(1,:), tri(2,:));

% DEBUG
d=fill(tri(1,:), tri(2,:), 'g');
delete(d);

% find the intersection of the triangle and the polygon
% if it is not empty, the side is invisible
% xb = polybool('and', tri(1,:), tri(2,:), pol(1,:), pol(2,:));
xb = polybooland(tri(1,:), tri(2,:), pol(1,:), pol(2,:));

if isempty(xb)
	% the vertices of side might be colinear with pt
	% if this happens, side is invisible
	if (pt(2)-pol(2,side(1)))/(pt(1)-pol(1,side(1))) == (pt(2)-pol(2,side(2)))/(pt(1)-pol(1,side(2)))
		b = 0;
	else
    	b = 1;
	end
else
    b = 0;
end
