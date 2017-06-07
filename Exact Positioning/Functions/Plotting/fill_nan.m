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

% Fill NaN delimited polygons
% 2017/1/30
function H = fill_nan( x , y , colorspec, opac, color )

if nargin < 5
    color = 'b';
end
if nargin < 4
    opac = 1;
end
if nargin < 3
    colorspec = 'b';
end

nan_indices = find( isnan(x) ); % The indices are the same for x and y
indices = [ 0 nan_indices length(x)+1 ];

% Collect the handles of all the plots
H = zeros(1, length( nan_indices )+1);
for i=1:length( nan_indices )+1
    tmpx = x( indices(i) + 1 : indices(i+1) - 1 );
    tmpy = y( indices(i) + 1 : indices(i+1) - 1 );

    h = fill(   tmpx , tmpy , ...
            colorspec, 'FaceColor', color, 'EdgeColor','None','facealpha',opac);

    hold on;
	H(i) = h;
end

