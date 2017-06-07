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

function x_new = no_collisions(region, x, i, uncertx, DEBUG)
% Keep node i inside its no collision cell by projecting the move vector on
% the boundary
if nargin < 5
    DEBUG = 0;
end
region_diameter = diameter(region);
N = length(x);

% Find the Voronoi cell of node i
Vcelli = region;
for j=1:N
    if j ~= i
        midpt = (x(:,i) + x(:,j)) / 2;
        theta = atan2( x(2,j)-x(2,i), x(1,j)-x(1,i) );
        % Create halfplane
        hp = [ 0 -region_diameter -region_diameter 0 ; 
            region_diameter region_diameter -region_diameter -region_diameter];
        % Rotate halfplane
        hp = rot(hp, theta);
        hp = bsxfun(@plus,hp,midpt);

        %%%%%%%%% DEBUG %%%%%%%%%
%         if DEBUG
%             plot_poly( hp, 'b');
%         end
        %%%%%%%%% DEBUG %%%%%%%%%

        % Intersect with the current voronoi cell
        [ipolyx , ipolyy] =...
                    polybooland(Vcelli(1,:), Vcelli(2,:),...
                                hp(1,:), hp(2,:) );
        Vcelli = [ipolyx ; ipolyy];
        if isempty(Vcelli)
            x_new = x(:,i);
            return;
        end
    end
end

%%%%%%%%% DEBUG %%%%%%%%%
    if DEBUG
        plot_poly( Vcelli, 'g');
    end
%%%%%%%%% DEBUG %%%%%%%%%

% Use keep in region to keep the disk inside its voronoi cell
% x_new = keep_in_region_old(Vcelli, x(:,i), uncertx, move_vector);
x_new = keep_in_region(Vcelli, x(:,i), uncertx);
