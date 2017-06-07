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

% offsets polygon P inward so that the distance between the edges is r
% P = [Px ; Py]
function G = offset(P, r, dir)

if nargin < 3
    dir = 'in';
end

% if strcmp(dir, 'out')
% 	dir = 'in';
% else
% 	dir = 'out';
% end

% Check if the last vertex is the same as the first
% If not add it
if ~isequal( P(:,1), P(:,end) )
    P = [P P(:,1)];
end

[Gx , Gy] = bufferm(P(1,:), P(2,:), r, dir);

% the points until the nan are the new polygon ?
nan_index = find( isnan(Gx) );

G = [Gx(1:nan_index-1)' ; Gy(1:nan_index-1)'];

