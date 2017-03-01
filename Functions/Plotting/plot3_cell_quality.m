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