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

function h = plot_AABB(limits, linespec)

if nargin == 1
    linespec = 'b';
end
xmin = limits(1);
xmax = limits(2);
ymin = limits(3);
ymax = limits(4);

x = [xmin xmax xmax xmin xmin];
y = [ymin ymin ymax ymax ymin];

h = plot(x, y, linespec, 'MarkerEdgeColor', 'none');
