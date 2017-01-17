% MIT License
% 
% Copyright (c) 2017 Sotiris Papatheodorou
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

% Plot NaN delimited polygons
% The first line of P contains the x coordinates of the polygon's vertices
% and the second line their y cordinate
% 2016/3/4
function plot_phi( phi , region, gridsize )

hh = ishold;
hold on

if nargin < 3
	gridsize = 50;
end

% Scale region so that no jagged edges appear
region_sc = offset_hom(region, 1.1);
% Create meshgrid
xl = linspace(min(region_sc(1,:)), max(region_sc(1,:)), gridsize);
yl = linspace(min(region_sc(2,:)), max(region_sc(2,:)), gridsize);
[xm, ym] = meshgrid(xl, yl);
f = zeros(size(xm));
for i=1:gridsize^2
	if inpolygon(xm(i), ym(i), region_sc(1,:), region_sc(2,:))
		f(i) = phi(xm(i), ym(i));
	end
end

% Create contour plot
f = f-max(max(f));
surf(xm, ym, f, 'EdgeColor', 'none');



% Create a large scaled version of the region
region_sc = offset_hom(region, 10);
[pbx, pby] = polybool('minus', region_sc(1,:), region_sc(2,:), ...
	region(1,:), region(2,:));
[F, V] = poly2fv(pbx, pby);
patch('Faces', F, 'Vertices', V, 'FaceColor', 'w', ...
  'EdgeColor', 'none')

if ~hh
	hold off;
end
