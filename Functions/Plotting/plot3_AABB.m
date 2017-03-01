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

function h = plot3_AABB(limits, linespec)

hh = ishold;
hold on
h = [];

if nargin == 1
    linespec = 'b';
end

xmin = limits(1);
xmax = limits(2);
ymin = limits(3);
ymax = limits(4);
zmin = limits(5);
zmax = limits(6);

% Bottom
x = [xmin xmax xmax xmin xmin];
y = [ymin ymin ymax ymax ymin];
z = [zmin zmin zmin zmin zmin];
htmp = plot3(x, y, z, linespec, 'MarkerEdgeColor', 'none');
h = [h htmp];

% Top
x = [xmin xmax xmax xmin xmin];
y = [ymin ymin ymax ymax ymin];
z = [zmax zmax zmax zmax zmin];
htmp = plot3(x, y, z, linespec, 'MarkerEdgeColor', 'none');
h = [h htmp];

% Sides
x = [xmin xmin];
y = [ymin ymin];
z = [zmin zmax];
htmp = plot3(x, y, z, linespec, 'MarkerEdgeColor', 'none');
h = [h htmp];
x = [xmax xmax];
y = [ymin ymin];
z = [zmin zmax];
htmp = plot3(x, y, z, linespec, 'MarkerEdgeColor', 'none');
h = [h htmp];
x = [xmax xmax];
y = [ymax ymax];
z = [zmin zmax];
htmp = plot3(x, y, z, linespec, 'MarkerEdgeColor', 'none');
h = [h htmp];
x = [xmin xmin];
y = [ymax ymax];
z = [zmin zmax];
htmp = plot3(x, y, z, linespec, 'MarkerEdgeColor', 'none');
h = [h htmp];

if ~hh
    hold off;
end