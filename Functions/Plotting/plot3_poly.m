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

% Plot NaN delimited polygons
% The first line of P contains the x coordinates of the polygon's vertices
% and the second line their y cordinate
% 2016/3/4
function H = plot3_poly( P , colorspec )

hh = ishold;
hold on

if nargin == 1
    colorspec = 'b';
end

if ~isempty(P)
    x = P(1,:);
    y = P(2,:);
	z = P(3,:);

    nan_indices = find( isnan(P(1,:)) ); % The indices are the same for x and y
    indices = [ 0 nan_indices length(x)+1 ];

    % Collect the handles of all the plots
    H = zeros(1, length( nan_indices )+1);
    for i=1:length( nan_indices )+1

        tmpP = [x( indices(i) + 1 : indices(i+1) - 1 ) ;
                y( indices(i) + 1 : indices(i+1) - 1 ) ;
				z( indices(i) + 1 : indices(i+1) - 1 )];
        % if the last vertex is not the same as the last, add the first vertex to
        % the end
        if ~isequal( tmpP(:,1), tmpP(:,end) )
            tmpP = [tmpP tmpP(:,1)];
        end

        h = plot3( tmpP(1,:), tmpP(2,:), tmpP(3,:), colorspec );

        H(i) = h;

	end

end

if ~hh
    hold off;
end
