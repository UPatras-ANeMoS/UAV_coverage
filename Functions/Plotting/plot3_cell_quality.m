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

function h = plot3_cell_quality(cell, quality, colorspec, opac, color)

if nargin < 3
	colorspec = 'b';
    opac = 1;
elseif nargin < 4
	opac = 1;
end

h = [];
% Set hold property
hh = ishold;
hold on

% Cell
Cx = cell(1,:);
Cy = cell(2,:);

% Set minimum and maximum z
zmin = 0;
zmax = quality;

% Fill top
if nargin < 5
	htmp = fill3(Cx, Cy, zmax*ones(size(Cx)), ...	
		colorspec ,'facealpha',opac);
else
	htmp = fill3(Cx, Cy, zmax*ones(size(Cx)), ...
		colorspec, 'FaceColor', color,'facealpha',opac);
end
h = [h htmp];

% Fill bottom
if nargin < 5
	htmp = fill3(Cx, Cy, zmin*ones(size(Cx)), ...
		colorspec ,'facealpha',opac);
else
	htmp = fill3(Cx, Cy, zmin*ones(size(Cx)), ...
		colorspec, 'FaceColor', color,'facealpha',opac);
end
h = [h htmp];

% Fill sides
if nargin < 5
	for i=1:length(Cx)-1
		Px = [Cx(i:i+1) fliplr(Cx(i:i+1))];
		Py = [Cy(i:i+1) fliplr(Cy(i:i+1))];
		Pz = [zmax zmax zmin zmin];
		htmp = fill3(Px, Py, Pz, ...
		colorspec ,'EdgeColor','None','facealpha',opac);
		h = [h htmp];
	end
else
	for i=1:length(Cx)-1
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