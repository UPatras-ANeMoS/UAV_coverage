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

function [uX, uY, uZ] = control_uniform_phi(region, phi, gridsize, zmin, zmax, a, ...
    W, C, f, i, Xi, Yi, Zi)

N = length(W);
uXY = zeros(2,1);
uZ = 0;
Ri = Zi * tan(a);

if ~isempty(W{i})
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
    % Wi is a closed list of points
    for k=1:length(Wi(1,:))-1
        % endpoints of the current line segment
        pt1 = Wi(:,k);
        pt2 = Wi(:,k+1);

        % Check if they are on the boundary. If both are on it dont
        % integrate
        [~, onB1] = inpolygon( pt1(1), pt1(2), region(1,:), region(2,:) );
        [~, onB2] = inpolygon( pt2(1), pt2(2), region(1,:), region(2,:) );

        if ~(onB1 && onB2)
            % Check if any of them is on Ci, if not, dont integrate
            [~, onCi1] = inpolygon( pt1(1), pt1(2), C{i}(1,:), C{i}(2,:) );
            [~, onCi2] = inpolygon( pt2(1), pt2(2), C{i}(1,:), C{i}(2,:) );

            if (onCi1 || onCi2)
                % Check if they are both inside Cj for all j in overlap
                % If they are, then this is a dominant arc
                free_arc = 1; % Free arc flag
                % Loop over all overlapping nodes
                for j=1:N
                    if i~=j
                        [inCj1] = inpolygon( pt1(1), pt1(2), C{j}(1,:), C{j}(2,:) );
                        [inCj2] = inpolygon( pt2(1), pt2(2), C{j}(1,:), C{j}(2,:) );

                        if inCj1 && inCj2
                            % This is a dominant arc, normal vector with
                            % magnitude fi-fj

                            d = norm( [pt1(1)-pt2(1) , pt1(2)-pt2(2)] );
                            n1 = (pt1-[Xi ; Yi]) / Ri;
                            n2 = (pt2-[Xi ; Yi]) / Ri;
                            nvector = (n1 + n2) / 2;
                            midpt = Ri * nvector + [Xi ; Yi];
                            phi_mpt = phi(midpt(1), midpt(2));

                            % X-Y control law
                            uXY = uXY + (f(i)-f(j)) * d * nvector * phi_mpt;

                            % Z control law
                            uZ = uZ + (f(i)-f(j)) * tan(a) * d * phi_mpt;
                        end

                        % If any of the points is inside a Cj, this is
                        % not a free arc
                        if inCj1 || inCj2
                            free_arc = 0;
                        end
                    end
                end % All other node for

                if free_arc
                    % This is a free arc, normal vector
                    d = norm( [pt1(1)-pt2(1) , pt1(2)-pt2(2)] );
                    n1 = (pt1-[Xi ; Yi]) / Ri;
                    n2 = (pt2-[Xi ; Yi]) / Ri;
                    nvector = (n1 + n2) / 2;
                    midpt = Ri * nvector + [Xi ; Yi];
                    phi_mpt = phi(midpt(1), midpt(2));

                    % X-Y control law
                    uXY = uXY + f(i) * d * nvector * phi_mpt;

                    % Z control law
                    uZ = uZ + f(i) * tan(a) * d * phi_mpt;
                end

            end
        end
    end % line segment for

    % Area integral for Z control law
    % Numerically integrate phi on Wi
    I = 0;
    [xm, ym] = meshgrid(linspace(Xi-Ri,Xi+Ri,gridsize), ...
        linspace(Yi-Ri,Yi+Ri,gridsize));
    dx = abs(xm(1,1)-xm(1,2));
    dy = abs(ym(1,1)-ym(2,1));
    for l=1:gridsize^2
        if inpolygon(xm(l), ym(l), W{i}(1,:), W{i}(2,:))
            ds = dx*dy;

            I = I + ds * phi(xm(l), ym(l));
        end
    end
    uZ = uZ + dfu(Zi, zmin, zmax) * I;
end

uX = uXY(1);
uY = uXY(2);
