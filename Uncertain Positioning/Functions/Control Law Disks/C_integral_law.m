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

function move_vector = C_integral_law( xi, GVcelli, sradiusi, DEBUG )
% Calculate the integral over Ci and return the move vector
% POLAR COORDINATES ARE (TH,R)
% Handle optional debug argument
if nargin == 3
    DEBUG = 0;
end

% If the cell is empty, dont move
if isempty(GVcelli)
    move_vector = [0 ; 0];
    return;
end
% Points Per Circle
ppc = 120;

% Translate the cell vertices so that point x is [0 0]'
celli_tr = bsxfun(@minus, GVcelli, xi);

% Create the sensing circle
C = create_circle([0 ; 0], sradiusi, ppc);

% Find intersection points of the cell and sensing circle
[xint, yint] = polyxpoly(celli_tr(1,:), celli_tr(2,:), C(1,:), C(2,:));
% Sort intersection points CCW, helps in finding valid pairs
[xint, yint] = poly2ccw(xint, yint);

% If there are no intersections or one intersection,
% the node either covers the whole cell, or covers as much as possible.
% Either way it must not move.
if length(xint) <= 1
    move_vector = [0 ; 0];
    return;
end
%%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
if DEBUG
    plot(xint+xi(1), yint+xi(2), 'mo');
end
%%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%

% Organize the points in pairs
% Each pair corresponds to an arc of the sensing circle
% If the arc is inside the voronoi cell it must be integrated
% The order of the points DOES matter
% The first is the start of the arc

% convert intersection points to polar coordinates
[ intp_pol(1,:), intp_pol(2,:) ] = cart2pol(xint, yint);
% Convert all angles from [-pi,pi] to [0,2pi]
neg_ang = intp_pol(1,:) < 0;
neg_ang = 2*pi*neg_ang;
intp_pol(1,:) = intp_pol(1,:) + neg_ang;
% add the first point as last as well
intp_pol = [intp_pol intp_pol(:,1)];


% Loop over each arc and find its midpoint
% If the midpoint is inside the cell, integrate the arc
move_vector = zeros(2,1);
for l=1:length(intp_pol)-1
    A = intp_pol(:,l);
    B = intp_pol(:,l+1);
    % If the endpoint has a smaller angle than the starting point increase
    % it
    if A(1) > B(1)
        B(1) = B(1) + 2*pi;
    end
    
    % Find arc midpoint
    arc_midpoint = [ (A(1) + B(1)) / 2 ; sradiusi ];
    % Return theta to [-pi,pi]
    if arc_midpoint(1) > 2*pi
        arc_midpoint(1) = arc_midpoint(1) - 2*pi;
    end
    if arc_midpoint(1) > pi
        arc_midpoint(1) = arc_midpoint(1) - 2*pi;
    end
    
    % Go back to cartesian to chek if it's inside the cell
    [amcx, amcy] = pol2cart( arc_midpoint(1), arc_midpoint(2) );
    if inpolygon( amcx, amcy, celli_tr(1,:), celli_tr(2,:) )
        % X part of control law
        % Make into a unit vector and multiply with the arc length
        arc_midpoint(2) = abs( sin(B(1)) - sin(A(1)) );
        % Add to the move vector
        [amcx, amcy] = pol2cart( arc_midpoint(1), arc_midpoint(2) );
        move_vector = move_vector + [amcx ; amcy];
        % Y part of control law
        % Make into a unit vector and multiply with the arc length
        arc_midpoint(2) = abs( cos(B(1)) - cos(A(1)) );
        % Add to the move vector
        [amcx, amcy] = pol2cart( arc_midpoint(1), arc_midpoint(2) );
        move_vector = move_vector + [amcx ; amcy];
        
        %%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
        if DEBUG
            m = norm([amcx ; amcy])/sradiusi;
            plot(amcx/m+xi(1), amcy/m+xi(2), 'mo');
        end
        %%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
    end
end
