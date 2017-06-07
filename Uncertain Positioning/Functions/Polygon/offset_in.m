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

% Offset the convex polygon P inwards by r
function Poffset = offset_in( P, r )

N = length(P(1,:));

Poffset = P;
for i=1:N
	v1 = P(:,i);
	v2 = P(:,mod(i,N)+1);
    edge = [v1 v2];
    % Check the direction of offset with the midpoint
    midpoint = (edge(:,1) + edge(:,2)) / 2;
    dir = 1;
    midpoint = point_perp(edge, midpoint, dir*10*eps);
    if ~inpolygon(midpoint(1), midpoint(2), P(1,:), P(2,:))
        dir = -1;
    end
    
    % Create the polygon to remove from P
    edge_offset = zeros(2,4);
    edge_offset(:,1) = point_perp(edge, edge(:,1), dir*r);
    edge_offset(:,2) = point_perp(edge, edge(:,2), dir*r);
    edge_offset(:,3) = point_perp(edge, edge(:,2), -dir*r);
    edge_offset(:,4) = point_perp(edge, edge(:,1), -dir*r);
    [xr , yr] = poly2cw(edge_offset(1,:), edge_offset(2,:));
    edge_offset = [xr ; yr];
    
    [x, y] = polybool('minus', Poffset(1,:), Poffset(2,:), ...
                            edge_offset(1,:), edge_offset(2,:));
    Poffset = [x ; y];
end
