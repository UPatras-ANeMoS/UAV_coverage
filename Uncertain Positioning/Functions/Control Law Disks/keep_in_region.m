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

% Projects the point x on the boundary of region taking into account the
% dimensions r.
function x_new = keep_in_region( region, x, r )

% offset the region inwards so that the disk of node x will be inside
region_offset = offset_in(region, r);


if inpolygon(x(1), x(2), region_offset(1,:), region_offset(2,:) )
    % Movement inside region
    x_new = x;
else
    % Movement outside region
    % Project x on the offset region
    x_new = project_pt_poly( x, region_offset );
end

