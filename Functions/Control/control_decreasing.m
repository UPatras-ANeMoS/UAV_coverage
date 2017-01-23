% MIT License
% 
% Copyright (c) 2017 Sotiris Papatheodorou
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function [uX, uY, uZ] = control_decreasing(region, zmin, zmax, a, b, ...
    W, C, f_u, i, Xi, Yi, Zi, gridsize)

N = length(W);
uXY = zeros(2,1);
uZ = 0;
Ri = Zi * tan(a);

% Keep only CW (external) contours
% If Wi has holes, remove the corresponding contours since they
% dont contribute to the control law
% Find NaN indices
nanindex = find( isnan( W{i}(1,:) ) );
if ~isempty( nanindex )
    Wi = []; % It will contain the external contour of Wi
    indices = [ 0 nanindex length( W{i}(1,:) )+1 ];
    for k=1:length(nanindex)+1
        % Keep a part of Wi
        tempx = W{i}(1, indices(k)+1 : indices(k+1)-1 );
        tempy = W{i}(2, indices(k)+1 : indices(k+1)-1 );
        if ispolycw(tempx, tempy)
            Wi = [Wi [tempx ; tempy]];
        end

    end
else
    % Wi has no holes
    Wi = W{i};
end


% Integrate over the arcs
% Loop over all line segments of Wi
% Wi is used instead of W{i} to include cases with NaNs

% Remove duplicate last vertex from Wi
if isequal( Wi(:,1), Wi(:,end) )
    Wi = Wi(:,1:end-1);
end
for k=1:length(Wi(1,:))
    % endpoints of the current line segment
    pt1 = Wi(:,k);
    pt2 = Wi(:,mod(k,length(Wi(1,:))) + 1);

    % Initialize arc
    arc = 0;

    % Check if they are on the boundary. If both are on it dont
    % integrate
    [~, onB1] = inpolygon( pt1(1), pt1(2), region(1,:), region(2,:) );
    [~, onB2] = inpolygon( pt2(1), pt2(2), region(1,:), region(2,:) );

    if ~(onB1 && onB2)
        % Check if they are on Ci
        [~, onCi1] = inpolygon( pt1(1), pt1(2), C{i}(1,:), C{i}(2,:) );
        [~, onCi2] = inpolygon( pt2(1), pt2(2), C{i}(1,:), C{i}(2,:) );
        if onCi1 && onCi2

            % This is either a free arc or a dominant arc
            % Loop over all overlaping nodes
            for j=1:N
                if i ~= j
                    [inCj1] = inpolygon( pt1(1), pt1(2), C{j}(1,:), C{j}(2,:) );
                    [inCj2] = inpolygon( pt2(1), pt2(2), C{j}(1,:), C{j}(2,:) );

                    if inCj1 && inCj2
                        arc = 2;
                        % Break keeping the value of j for later use
                        break
                    end
                end % overlap check
            end % Loop over all nodes in overlap

            if arc ~= 2
                % This is a free arc
                arc = 1;
            end
        end % On Ci check
    end % Boundary of omega check



    % Find the control law depending on the arc
    switch arc
        case 0
            % Boundary or arc of Cj
        case 1
            % Free arc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Normal vector on Ci
            n1 = (pt1-[Xi ; Yi]) / Ri;
            n2 = (pt2-[Xi ; Yi]) / Ri;
            nvector = (n1 + n2) / 2;
            dl = norm( [pt1(1)-pt2(1) , pt1(2)-pt2(2)] );
            % Do not use fp because of numerical accuracy
            fi = b * f_u(i);

            % X-Y control law
            uXY = uXY + fi * dl * nvector;

            % Z control law
            uZ = uZ + fi*tan(a)*dl;

        case 2
            % Dominant arc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This is a dominant arc, normal vector with
            % magnitude fi-fj
            n1 = (pt1-[Xi ; Yi]) / Ri;
            n2 = (pt2-[Xi ; Yi]) / Ri;
            nvector = (n1 + n2) / 2;
            dl = norm( [pt1(1)-pt2(1) , pt1(2)-pt2(2)] );
            % Do not use fp because of numerical accuracy
            fi = b * f_u(i);
            fj = b * f_u(j);
            % The value of j has been kept from the break statement

            % X-Y control law
            uXY = uXY + (fi-fj) * dl * nvector;

            % Z control law
            uZ = uZ + (fi-fj)*tan(a)*dl;
            
    end % End of arc selection switch
end % Loop over each edge of Wi

% Area integrals
% Create grid
minxWi = min( W{i}(1,:) );
maxxWi = max( W{i}(1,:) );
minyWi = min( W{i}(2,:) );
maxyWi = max( W{i}(2,:) );
mingrid = min( [minxWi minyWi] );
maxgrid = max( [maxxWi maxyWi] );
lx = linspace(mingrid, maxgrid, gridsize);
ly = linspace(mingrid, maxgrid, gridsize);
[gx, gy] = meshgrid( lx, ly );
dx = lx(2)-lx(1);

% Start integration
Ix = 0;
Iy = 0;
Iz = 0;
for k=1:gridsize^2
    if inpolygon( gx(k), gy(k), W{i}(1,:), W{i}(2,:) )
        Ix = Ix + dx^2 * dfp_dxi(gx(k), gy(k), Xi, Yi, Zi, zmin, zmax, a, b);
        % Use the same function as for xi but with swapped inputs
        Iy = Iy + dx^2 * dfp_dxi(gy(k), gx(k), Yi, Xi, Zi, zmin, zmax, a, b);
        Iz = Iz + dx^2 * dfp_dzi(gx(k), gy(k), Xi, Yi, Zi, zmin, zmax, a, b);
    end
end

uX = uXY(1) + Ix;
uY = uXY(2) + Iy;
uZ = uZ + Iz;
