% MIT License
% 
% Copyright (c) 2016 Sotiris Papatheodorou
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

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