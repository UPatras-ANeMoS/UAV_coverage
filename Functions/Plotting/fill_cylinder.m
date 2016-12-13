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

function h = fill_cylinder(center_top, height, radius, colorspec, opac, color)

if nargin < 4
	colorspec = 'b';
    opac = 1;
elseif nargin < 5
	opac = 1;
end

h = [];
% Set hold property
hh = ishold;
hold on

% Create circle
t = linspace(0,2*pi,60);
Cx = radius*cos(t) + center_top(1);
Cy = radius*sin(t) + center_top(2);

% Set minimum and maximum z
zmin = (center_top(3)-height);
zmax = center_top(3);

% Fill top
if nargin < 6
	htmp = fill3(Cx, Cy, zmax*ones(size(Cx)), ...	
		colorspec ,'facealpha',opac);
else
	htmp = fill3(Cx, Cy, zmax*ones(size(Cx)), ...
		colorspec, 'FaceColor', color,'facealpha',opac);
end
h = [h htmp];

% Fill bottom
if nargin < 6
	htmp = fill3(Cx, Cy, zmin*ones(size(Cx)), ...
		colorspec ,'facealpha',opac);
else
	htmp = fill3(Cx, Cy, zmin*ones(size(Cx)), ...
		colorspec, 'FaceColor', color,'facealpha',opac);
end
h = [h htmp];

% Fill sides
if nargin < 6
	for i=1:length(t)-1
		Px = [Cx(i:i+1) fliplr(Cx(i:i+1))];
		Py = [Cy(i:i+1) fliplr(Cy(i:i+1))];
		Pz = [zmax zmax zmin zmin];
		htmp = fill3(Px, Py, Pz, ...
		colorspec ,'EdgeColor','None','facealpha',opac);
		h = [h htmp];
	end
else
	for i=1:length(t)-1
		Px = [Cx(i:i+1) fliplr(Cx(i:i+1))];
		Py = [Cy(i:i+1) fliplr(Cy(i:i+1))];
		Pz = [zmax zmax zmin zmin];
		htmp = fill3(Px, Py, Pz, ...
		colorspec ,'FaceColor', color, 'EdgeColor','None','facealpha',opac);
		h = [h htmp];
	end
end






% Reset hold property
if ~hh
	hold off
end