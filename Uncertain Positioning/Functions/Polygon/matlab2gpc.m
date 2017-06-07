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

function matlab2gpc(polys, filename, write_single_polygon, write_hole_flags)

file1 = fopen(filename,'w');

% if no option is provided, assume hole flags should not be written
if nargin == 2
    write_single_polygon = 0;
    write_hole_flags = 0;
elseif nargin == 3
    write_hole_flags = 0;
end

% if read_multiple_polygons is true, the first number in the file is the
% number of polygons
if write_single_polygon == 0
    num_polys = length(polys);
    fprintf(file1,'%d\n',num_polys);
else
    num_polys = 1;
end

% loop over each polygon
for i=1:num_polys
    
    % number of contours, make it use NaNs - NOT COMPLETE
    % Only works for single contour polygons
    num_contours_i = 1; % find the number of contours from the number of NaNs in polys
    fprintf(file1,'%d\n',num_contours_i);
    
    % loop over each contour of polygon i
    for j=1:num_contours_i
        
        % number of vertices of contour j
        num_vertices_j = length( polys{i}(1,:) );
        fprintf(file1,'%d\n',num_vertices_j);
        
        % is hole of contour j (can be turned off)
        if write_hole_flags == 1
            % add correct hole detection - NOT COMPLETE
            is_hole_j = 0;
            fprintf(file1,'%d\n',is_hole_j);
        end
        
        % If the last vertex is the same as the first, remove it
        if isequal(polys{i}(:,1), polys{i}(:,end))
            polys{i} = polys{i}(:,1:end-1);
        end
        
        % loop over every point coordinate of contour j
        for k=1:num_vertices_j
            fprintf(file1,'%f %f\n', polys{i}(1,k), polys{i}(2,k));
        end
        
    end
    
end
        
