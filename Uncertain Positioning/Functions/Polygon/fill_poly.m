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
function h = fill_poly(P , colorspec )

if nargin < 2
    colorspec = 'b';
end

x = P(1,:);
y = P(2,:);

nan_indices = find( isnan(x) ); % The indices are the same for x and y
indices = [ 0 nan_indices length(x)+1 ];

h = [];

for i=1:length( nan_indices )+1

    ht = fill(   x( indices(i) + 1 : indices(i+1) - 1 ) , ...
            y( indices(i) + 1 : indices(i+1) - 1 ) , ...
            colorspec ,'EdgeColor','None');
    
    h = [h ht];
    hold on;

end

