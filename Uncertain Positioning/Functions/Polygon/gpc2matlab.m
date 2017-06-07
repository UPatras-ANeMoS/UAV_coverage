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

function polys = gpc2matlab(filename, read_single_polygon, read_hole_flags)

% if no option is provided, assume hole flags should be read
if nargin == 1
    read_single_polygon = 0;
    read_hole_flags = 0;
elseif nargin == 2
    read_hole_flags = 0;
end

a = importdata(filename);

% the data index points to the first element of the data that hasn't been
% accessed
data_index = 1;

% if read_multiple_polygons is true, the first number in the file is the
% number of polygons
if read_single_polygon == 0
    num_polys = a(1);
    data_index = data_index + 1;
else
    num_polys = 1;
end

polys = cell( [1 num_polys] );

% loop over each polygon
for i=1:num_polys
    
    num_contours_i = a(data_index);
    data_index = data_index + 1;
    
    % loop over each contour of polygon i
    for j=1:num_contours_i
        
        % number of vertices of contour j
        num_vertices_j = a(data_index);
        data_index = data_index + 1;
        
        % is hole of contour j (can be turned off)
        if read_hole_flags == 1
            is_hole_j = a(data_index);
            data_index = data_index + 1;
        else
            % if no hole flag is present, assume it's not a hole
            is_hole_j = 0;
        end
        
        contour_j = zeros(2,num_vertices_j);
        
        % loop over every point coordinate of contour j
        for k=1:2*num_vertices_j
            contour_j(k) = a(data_index);
            data_index = data_index + 1;
        end
        
        % make CW or CCW depending on is_hole_j
        if is_hole_j == 0
            [ contour_j(1,:), contour_j(2,:) ] = poly2cw(...
                contour_j(1,:), contour_j(2,:) );
        else
            [ contour_j(1,:), contour_j(2,:) ] = poly2ccw(...
                contour_j(1,:), contour_j(2,:) );
        end
        
        % add the contour to the contour list (the cell of polygon i)
        polys{i} = [ polys{i} contour_j ];
        
        % if there is another contour to be added, add NaNs to the cell
        if j ~= num_contours_i
            polys{i} = [ polys{i} [NaN ; NaN] ];
        end
        
        
    end
    
end

% if only one polygon was read, return it as an array
% if read_single_polygon == 1
%     polys = polys{1};
% end
        
