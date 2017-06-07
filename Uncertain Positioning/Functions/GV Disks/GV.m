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

% Constructs the Guaranteed Voronoi diagram of the input disks
function [GVcells, rcells, covered_area, overlap] = ...
GV(region_vert, max_line, x, uncert, guarsrad)
%%%%%%%%%%%%%%%%%%%%%% TO DO %%%%%%%%%%%%%%%%%%%%%%
% It doesnt work on three or more colinear nodes
%%%%%%%%%%%%%%%%%%%%%% TO DO %%%%%%%%%%%%%%%%%%%%%%

% Number of nodes
N = length( x(1,:) );

% Voronoi Cells -----------------------------------------------------------
p = 201; % points per edge
% increase max_line
max_line = 1.01 * max_line;

[GVcells, overlap] = guar_voronoi_cells(region_vert , x , uncert, p, max_line);


% Radius constrained cells ------------------------------------------------
covered_area = 0;
rcells = cell( [N 1] ); % Initialization
% parfor - slower
for i=1:N
    if ~isempty( GVcells{i} )
        rcells{i} = rad_cell( x(:,i) , GVcells{i} , guarsrad(i) );
        
        % Add covered area
        if ~isempty( rcells{i} )
            covered_area = covered_area + polyarea( rcells{i}(1,:), rcells{i}(2,:) );
        end
        
    end
end




