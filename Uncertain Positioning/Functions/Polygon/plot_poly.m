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

% plot NaN delimited polygons
function h = plot_poly( P , colorspec )
if isempty(P)
    return
end

a = ishold;

if nargin == 1
    colorspec = 'b';
end
% if the last vertex is not the same as the last, add the first vertex to
% the end
% if ~isequal( P(:,1), P(:,end) )
%     P = [P P(:,1)];
% end

x = P(1,:);
y = P(2,:);

nan_indices = find( isnan(P(1,:)) ); % The indices are the same for x and y
indices = [ 0 nan_indices length(x)+1 ];

h = [];

for i=1:length( nan_indices )+1

    ht = plot(   [x( indices(i) + 1 : indices(i+1) - 1 ) x(indices(i) + 1)] , ...
            [y( indices(i) + 1 : indices(i+1) - 1 ) y(indices(i) + 1)] , ...
            colorspec );

    h = [h ht];
    hold on;

end

% Keep the previous hold state
if a
    hold on
else
    hold off
end
