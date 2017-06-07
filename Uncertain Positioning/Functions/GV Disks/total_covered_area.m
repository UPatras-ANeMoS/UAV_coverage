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

% Returns the total area the nodes cover, not just the GV cells
function covered_area = total_covered_area(region_vert, x, guarsrad)

N = length(x(1,:));

sensor_circles = cell( [N 1] );
circle_union = zeros(2,2);

for i=1:N
    sensor_circles{i} = create_circle(x(:,i), guarsrad(i));
    
    [polyboolx, polybooly] = polybool('or',...
        circle_union(1,:), circle_union(2,:),...
        sensor_circles{i}(1,:), sensor_circles{i}(2,:) );
    
    circle_union = [polyboolx ; polybooly];
end

% intersect circle union with region Omega
[polyboolx, polybooly] = polybool('and',...
    circle_union(1,:), circle_union(2,:),...
    region_vert(1,:), region_vert(2,:) );

region_intersection = [polyboolx ; polybooly];


covered_area = polyarea_nan(region_intersection(1,:), region_intersection(2,:));

