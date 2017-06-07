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

% makes sure the move_vector results in movement inside region omega
function x_new = keep_in_region_old(region_vert, x, uncertx, move_vector)

% offset the region inwards so that the disk of node x will be inside
region_offset = offset(region_vert, uncertx, 'in');
if isempty(region_offset)
    x_new = x;
    return;
end
%%%%%%%%% DEBUG %%%%%%%%%
% plot_poly( region_offset, 'm');
%%%%%%%%% DEBUG %%%%%%%%%

% use uncertx to guarantee the disk will not move out of omega
% Maybe use polyxpoly instead of inpolygon (might be slower)
x_new = x + move_vector;

if inpolygon(x_new(1), x_new(2), region_offset(1,:), region_offset(2,:) )
    % movement inside omega
    return;
    
else
    % movement outside omega
    
    % move the node to the intersection point
    % ONLY ONE INTERSECTION POINT, ok for CONVEX region
    [xi, yi, ii] = polyxpoly([x_new(1) x(1)], [x_new(2) x(2)],...
        region_offset(1,:), region_offset(2,:) );
    
    % if x is on the boundary of omega, polyxpoly sometimes returns no
    % intersection
    if isempty(xi)
        x_new = x;
        % should find which edge x is on but polyxpoly doesnt work
        return;
        
    else
        x_new = [xi(1) ; yi(1)];
        
        % subtract the above move from the move vector
        move_vector = move_vector - (x_new - x);

        % find the edge where the intersection happened
        % ONLY ONE INTERSECTION POINT
        if isequal( x_new, region_offset(:,ii(2)) )
            % the intersection happens at a vertex of omega, no other move
            % is possible
            return;
        else
            edge = region_offset(:, ii(2)+1) - region_offset(:, ii(2));
            move_vector = projv(move_vector, edge);

            % check if the new move vector is still inside the polygon
            x_temp = x_new + move_vector;
            % this check isnt always correct
            if inpolygon(x_temp(1), x_temp(2), region_offset(1,:), region_offset(2,:) )
                x_new = x_temp;
                return;
            else
                % find the vertex of omega the node must move to
                [xi, yi] = polyxpoly([x_temp(1) x_new(1)], [x_temp(2) x_new(2)],...
                    region_offset(1,:), region_offset(2,:) );
                % one of the two vertices returned will be x_new and the
                % other one the desired
                if isempty(xi)
                    % x_temp is on the boundary of region_offset and
                    % inpolygon isnt accurate enough
                    x_new = x_temp;
                    return;
                elseif isequal( [xi(1) ; yi(1)] , x_new) || length(xi) == 1
                    x_new = [xi(1) ; yi(1)];
                else
                    x_new = [xi(2) ; yi(2)];
                end
            end
        end
        
        
    end
    
    
end

