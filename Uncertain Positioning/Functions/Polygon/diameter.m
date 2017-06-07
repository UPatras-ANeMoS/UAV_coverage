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

% find the largest line segment that can fit inside polygon polv
% poly_vert = [x ; y] the vertices of the polygon
function maxl = diameter(poly_vert)

verts = length( poly_vert(1,:) ); % number of vertices
maxl = 0; % maximum line segment length
for i=1:verts
    for j=1:verts
        d = eucl_dist( poly_vert(:,i) , poly_vert(:,j) );
        if d > maxl
            maxl = d;
        end
    end
end

