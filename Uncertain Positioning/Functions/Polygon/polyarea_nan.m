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

% Finds the area of a polygon with holes
% polygon contours are Nan delimited
function area = polyarea_nan(x, y)
% The outer contour is the last NaN delimited segment

nan_indices = find(isnan( x ));
indices = [ 0 nan_indices length(x)+1 ];
area = 0;

if ~isempty(nan_indices)
    for i=1:length(nan_indices)+1

        tempx = x( indices(i)+1 : indices(i+1)-1 );
        tempy = y( indices(i)+1 : indices(i+1)-1 );
        if ispolycw(tempx, tempy)
        	% external contour
        	area = area + polyarea(tempx, tempy);
    	else
        	% internal contour
        	area = area - polyarea(tempx, tempy);
    	end

    end
else
    area = polyarea(x, y);
end

