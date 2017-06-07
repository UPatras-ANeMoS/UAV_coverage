% Copyright 2017 Sotiris Papatheodorou
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%    http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% Return the Additively Weighted Guaranteed Voronoi diagram of node i
% region is the region within which the cell is constrained
% q is the vector of nodes
% r is the vector of uncertainty radii
% R is the vector of weights
function W = AWGV_cell(region, q, r, R, i, diam)
	if nargin < 6
		diam = diameter(region);
	end
	
	% Number of nodes
	N = length(q(1,:));
	
	% Initialize cell to region
	W = region;
	% Loop over all other nodes
	for j=1:N
		if i ~= j
			% Calculate hyperbola parameters
			a = (r(i) + r(j) + R(j) - R(i))/2;
			c = norm(q(:,i)-q(:,j))/2;
			b = sqrt(c^2-a^2);
%             fprintf('a = %f  b = %f  c = %f\n', a, b, c);
			
            % The cell of i is empty (or a ray for ==)
            if a >= c
                W = [];
                return
            % The cell of i is the whole region (or minus a ray for ==)
            elseif a <= -c
                % Node j doesn't affect node i's cell, continue with the
                % next node
                continue
            end
            
            if a == 0
                H = [0 -diam -diam 0 ; diam diam -diam -diam];
            else
                % Create the cell of i wrt j
                t = linspace(asinh(diam/abs(a)),-asinh(diam/abs(a)), 100);
                H = [a*cosh(t) ; b*sinh(t)];
                % Extend the cell towards qi
                H = [H(:,1)+[0;diam] H(:,1)+[2*diam;diam] ...
                    H(:,end)+[2*diam;-diam] H(:,end)+[0;-diam] H];
            end
            % Rotate and translate the cell
			H = rot(H, atan2(q(2,i)-q(2,j),q(1,i)-q(1,j)));
			H = bsxfun(@plus, H, (q(:,i)+q(:,j))/2);
			
			
			% Intersect with current cell
			if ~isempty(W)
				[pbx, pby] = polybool('and', ...
                    W(1,:), W(2,:), H(1,:), H(2,:));
				W = [pbx ; pby];
			else
				return
			end
		end
	end
end

function B = rot(A,theta)
	R = [    cos(theta)  -sin(theta);
				sin(theta)  cos(theta)];
	B = zeros( size(A) );
	N = length(A(1,:));
	for i=1:N
		B(:,i) = R * A(:,i);
	end
end

function maxl = diameter(poly)
    verts = length( poly(1,:) ); % number of vertices
    maxl = 0; % maximum line segment length
    for i=1:verts
        for j=1:verts
            d = norm( poly(:,i) - poly(:,j) );
            if d > maxl
                maxl = d;
            end
        end
    end
end