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

% Calculate the radius r constrained cell of point x (a column vector)
function rcell = rad_cell( x , cell , r)

if ~isempty(cell)
    % Define certain sensing circle
    t = linspace(0, 2*pi, 120); % Number of points for circles
    circx = r * cos(t) + x(1);
    circy = r * sin(t) + x(2);
    [circx, circy] = poly2cw(circx, circy);

    [ rcellx , rcelly ] = ...
    polybool( 'and', circx, circy, ...
    cell(1,:), cell(2,:) );

    rcell = [rcellx ; rcelly];
else
    rcell = [];
end
