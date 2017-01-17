% MIT License
% 
% Copyright (c) 2016-2017 Sotiris Papatheodorou
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD CHECKs FOR EMPTY CELLS Wi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear variables
close all

PLOT_STATE_3D = 0;
PLOT_STATE_2D = 1;
PLOT_STATE_QUALITY = 0;
SAVE_PLOTS = 0;

% Add function path
addpath( genpath('Functions') );

% Bullo region
Xb=[ 0, 2.125, 2.9325, 2.975, 2.9325, 2.295, 0.85, 0.17 ];
Yb=[ 0, 0, 1.5, 1.6, 1.7, 2.1, 2.3, 1.2 ];
[Xb, Yb] = poly2cw(Xb, Yb);
region = [Xb ; Yb];
region_area = polyarea( Xb, Yb );
omega_diameter = diameter( region );

% ----------------- Network parameters -----------------
% Altitude constraints
zmin = 0.3;
zmax = 2.3;

% Quality-coverage tradeoff
Q = 1;

% Sensing cone angle (half the angle of the cone)
a = 20*pi/180;

% Sensing quality for zmin at the boundary of Ci
b = 0.5;

% Optimal altitude
% zopt = z_optimal_decreasing(zmin, zmax, a, b);
% zopt = 1.563;
zopt = z_optimal_decreasing(zmin, zmax);


% Initial positions - 3 nodes - 12 seconds
% X = [0.40, 0.60, 0.55];
% Y = [0.50, 0.60, 0.50];
% Z = [0.45, 0.55, 0.50];

% Initial positions - 9 nodes - 24 seconds
% X = [0.40, 0.60, 0.55, 0.60, 0.50, 0.70, 0.60, 0.90, 0.80];
% Y = [0.50, 0.60, 0.50, 0.40, 0.60, 0.50, 0.75, 0.85, 0.95];
% Z = [0.45, 0.55, 0.50, 0.60, 0.40, 0.52, 0.57, 0.63, 0.65];

% Initial positions - 9 nodes - smax = 50+
% X = [1, 1.2, 1.4, 1.6, 1.8, 0.3, 2, 2.1, 0.6];
% Y = [1, 1.5, 1.2, 0.8, 1.1, 0.5, 1.5, 1.3, 0.4];
% Z = [0.4, 0.6, 0.7, 0.8, 1, 0.38, 0.35, 0.5, 0.55];
% Z = [0.4, 0.6, 0.7, 0.8, 1, 0.55, 0.35, 0.5, 0.55]; % fi = fj

% Initial positions 
% X = [0.4, 0.7, 1.7, 1.8, 1.2, 1.7, 1.8];
% Y = [0.5, 0.4, 0.5, 0.4, 1.5, 1.9, 1.4];
% Z = [0.7, 0.7, 1.3, 0.6, 1.3, 1.0, 0.7];

% Initial positions - Sensing in sensing - movement
% X = [1, 1.1];
% Y = [0.8, 0.85];
% Z = [0.38, 0.75];
% zmax = 1.3;

% Initial positions - Sensing in sensing - no movement
% X = [1, 1.05];
% Y = [0.8, 0.8];
% Z = [0.38, 0.75];
% zmax = 1.3;

% Initial positions - Same zi
% X = [1, 1.05];
% Y = [0.8, 0.8];
% Z = [0.5, 0.5];

% Initial positions - Close zi
% X = [1, 1.05];
% Y = [0.8, 0.8];
% Z = [0.5, 0.6];

% Initial position at zmax
X = [1.6];
Y = [1.1];
Z = [zmax-0.01];
Z = [0.9];

% All cases
% X = [0.4, 0.7, 1.7, 1.8, 1.2, 1.7, 1.8];
% Y = [0.5, 0.4, 0.5, 0.4, 1.5, 1.9, 1.4];
% Z = [0.7, 0.7, 1.3, 0.6, 1.3, 1.0, 0.7];

% Number of nodes
N = length(X);

% Sensing radii
R = tan(a) * Z;


% ----------------- Simulation parameters -----------------
% Simulation steps
smax = 50;
% Simulation time step
Tstep = 0.1;
% Points Per Circle
PPC = 60;
% Grid size for integrals (size of side)
gridsize = 100;
% Radius for disks on plots
disk_rad = 0.02;
% vector for circle parametrization
t = linspace(0, 2*pi, PPC+1);
t = t(1:end-1); % remove duplicate last element
t = fliplr(t); % flip to create CW ordered circles
% Simulation data storage
Xs = zeros(smax, N);
Ys = zeros(smax, N);
Zs = zeros(smax, N);
Rs = zeros(smax, N);
cov_area = zeros(smax,1);
H = zeros(smax,1);
% Initialize (cell) arrays
% f is just the zi dependent component of the coverage quality
f_u = zeros(1, N);
% Sensing disks
C = cell([1 N]);
% Sensed space partitioning
W = cell([1 N]);
% Overlap is used to store which sensing disks overlap with eachother
overlap = zeros(N, N);





% ----------------- Simulation -----------------
if PLOT_STATE_3D || PLOT_STATE_2D || PLOT_STATE_QUALITY
	figure
end
tic;
for s=1:smax
    fprintf('%.2f%% complete\n',100*s/smax);
    
    % Sensing radii
    R = tan(a) * Z;
    
    % Coverage quality
    f_u = fu(Z, zmin, zmax);
    % Sensing disks
    for i=1:N
        C{i} = [X(i) + R(i) * cos(t) ; Y(i) + R(i) * sin(t)];
        % Initialize cells to sensing disks
        W{i} = C{i};
    end
    
    % Store simulation data
    Xs(s,:) = X;
    Ys(s,:) = Y;
    Zs(s,:) = Z;
    Rs(s,:) = R;
    
    % Sensed space partitioning
    for i=1:N
        % Loop over all other nodes
        for j=1:N
            if j ~= i
                % Degenerate case, split the common region at the
                % intersection points
                if Z(i) == Z(j)
                    % Find intersection points
                    [xx, yy, ii] = polyxpoly( C{i}(1,:), C{i}(2,:),...
                                C{j}(1,:), C{j}(2,:) );
                    % only if the sensing circles intersect
                    if ~isempty(xx)
                        % set the corresponding elements of the overlap matrix
                        overlap(i,j) = 1;
                        overlap(j,i) = 1;
                        
                        % Add the intersection points
                        Ci = C{i};
                        for k=1:length(xx)
                            % assuming ii is sorted
                            Ci = [Ci(:,1:ii(k)+(k-1)) [xx(k) ; yy(k)] ...
                                Ci(:,ii(k)+k:end)];
                        end
                        % Remove the points of Ci that are inside Cj
                        [inCj, onCj] = inpolygon( Ci(1,:), Ci(2,:),...
                            C{j}(1,:), C{j}(2,:) );
                        Ci = Ci(:, ~inCj | onCj);
                        % Polybool the new Ci with Wi
                        [pbx, pby] = polybool( 'and', W{i}(1,:), W{i}(2,:),...
                            Ci(1,:), Ci(2,:) );
                        W{i} = [pbx ; pby];
                    end
                % Remove a portion of the sensing disk if zi ~= zj
                else
                    % Overlap check
                    [pbx, ~] = polybool( 'and', C{i}(1,:), C{i}(2,:),...
                    C{j}(1,:), C{j}(2,:) );
                    % Only change Wi on overlap
                    if ~isempty(pbx)
                        % set the corresponding elements of the overlap matrix
                        overlap(i,j) = 1;
                        overlap(j,i) = 1;
                        
                        % Create the Intersection Circle
                        Ki = (1-b) / R(i)^2;
                        Kj = (1-b) / R(j)^2;
                        ICx = (Ki*f_u(i)*X(i) - Kj*f_u(j)*X(j)) / (Ki*f_u(i) - Kj*f_u(j));
                        ICy = (Ki*f_u(i)*Y(i) - Kj*f_u(j)*Y(j)) / (Ki*f_u(i) - Kj*f_u(j));
                        ICr = sqrt( ICx^2 + ICy^2 + (-Ki*f_u(i)*(X(i)^2 + Y(i)^2)...
                            + Kj*f_u(j)*(X(j)^2 + Y(j)^2) + f_u(i)-f_u(j)) / (Ki*f_u(i) - Kj*f_u(j)));
                        % If the radius of IC is large, use more points for
                        % the creation of IC
                        ep = floor(2*ICr/(R(i)+R(j)));
                        tt = linspace(0, 2*pi, (ep+1)*PPC+1);
                        tt = tt(1:end-1); % remove duplicate last element
                        tt = fliplr(tt); % flip to create CW ordered circles
                        IC = [ICx + ICr * cos(tt) ; ICy + ICr * sin(tt)];

%                         plot_poly(IC, 'g');
%                         hold on

                        % If i is the dominant node, intersect
                        if Z(i) < Z(j)
                                [pbx, pby] = polybool( 'and', W{i}(1,:), W{i}(2,:),...
                                    IC(1,:), IC(2,:) );
                        % else remove the intersection of Cj with IC
                        else
                                [pbx, pby] = polybool( 'and', C{j}(1,:), C{j}(2,:),...
                                    IC(1,:), IC(2,:) );
                                [pbx, pby] = polybool( 'minus', W{i}(1,:), W{i}(2,:),...
                                    pbx, pby );
                        end

                        % Change the currecnt cell
                        W{i} = [pbx ; pby];

                    end % overlap check
                    
                end % Z(i),Z(j) check
            end % j~=i
        end % for j
        
        % AND with the region omega
        [pbx, pby] = polybool( 'and', W{i}(1,:), W{i}(2,:), Xb, Yb );
        W{i} = [pbx ; pby];
    end
    
    
    % Find covered area and H objective
    for i=1:N
        cov_area(s) = cov_area(s) + polyarea_nan(W{i}(1,:), W{i}(2,:));
        
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
        
        for k=1:gridsize^2
            if inpolygon( gx(k), gy(k), W{i}(1,:), W{i}(2,:) )
                 H(s) = H(s) + dx^2 * fp(gx(k), gy(k), X(i), Y(i), Z(i), zmin, zmax, a, b);
            end
        end
    end
    cov_area(s) = cov_area(s)/region_area;
    
    
    
    % ----------------- Control law -----------------
    move_vectors = zeros(2, N); % [uX ; uY]
    uZ = zeros(1,N);
    for i=1:N
        
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
            [~, onB1] = inpolygon( pt1(1), pt1(2), Xb, Yb );
            [~, onB2] = inpolygon( pt2(1), pt2(2), Xb, Yb );

            if ~(onB1 && onB2)
                % Check if they are on Ci
                [~, onCi1] = inpolygon( pt1(1), pt1(2), C{i}(1,:), C{i}(2,:) );
                [~, onCi2] = inpolygon( pt2(1), pt2(2), C{i}(1,:), C{i}(2,:) );
                if onCi1 && onCi2
                    
                    % This is either a free arc or a dominant arc
                    % Loop over all overlaping nodes
                    for j=1:N
                        if overlap(i,j)
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
                    n1 = (pt1-[X(i) ; Y(i)]) / R(i);
                    n2 = (pt2-[X(i) ; Y(i)]) / R(i);
                    nvector = (n1 + n2) / 2;
                    dl = norm( [pt1(1)-pt2(1) , pt1(2)-pt2(2)] );
					% Do not use fp because of numerical accuracy
					fi = b * f_u(i);

                    % X-Y control law
                    move_vectors(:,i) = move_vectors(:,i) + fi * dl * nvector;

                    % Z control law
                    uZ(i) = uZ(i) + fi*tan(a)*dl;

                    % DEBUG PLOTS
%                     plot( (pt1(1)+pt2(1))/2, (pt1(2)+pt2(2))/2, 'r.');
%                     hold on
                case 2
                    % Dominant arc %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % This is a dominant arc, normal vector with
                    % magnitude fi-fj
                    n1 = (pt1-[X(i) ; Y(i)]) / R(i);
                    n2 = (pt2-[X(i) ; Y(i)]) / R(i);
                    nvector = (n1 + n2) / 2;
                    dl = norm( [pt1(1)-pt2(1) , pt1(2)-pt2(2)] );
					% Do not use fp because of numerical accuracy
					fi = b * f_u(i);
					fj = b * f_u(j);
                    % The value of j has been kept from the break statement

                    % X-Y control law
                    move_vectors(:,i) = move_vectors(:,i) + (fi-fj) * dl * nvector;

                    % Z control law
                    uZ(i) = uZ(i) + (fi-fj)*tan(a)*dl;

                    % DEBUG PLOTS
%                     plot( (pt1(1)+pt2(1))/2, (pt1(2)+pt2(2))/2, 'b.');
%                     hold on
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
                Ix = Ix + dx^2 * dfp_dxi(gx(k), gy(k), X(i), Y(i), Z(i), zmin, zmax, a, b);
                % Use the same function as for xi but with swapped inputs
                Iy = Iy + dx^2 * dfp_dxi(gy(k), gx(k), Y(i), X(i), Z(i), zmin, zmax, a, b);
                Iz = Iz + dx^2 * dfp_dzi(gx(k), gy(k), X(i), Y(i), Z(i), zmin, zmax, a, b);
            end
        end
        %%%%%%%%%%%% DEBUG %%%%%%%%%%%%
% 		fprintf('r%d = %.8f\n', i, R(i))
%         fprintf('z%d = %.8f\n', i, Z(i))
%         fprintf('Ib = %.8f  Ic = %.8f\n', uZ(i), Iz)
        %%%%%%%%%%%% DEBUG %%%%%%%%%%%%
        
        move_vectors(1,i) = move_vectors(1,i) + Ix;
        move_vectors(2,i) = move_vectors(2,i) + Iy;
        uZ(i) = uZ(i) + Iz;
        
    end % node for
    
    % Control inputs
    uX = move_vectors(1,:);
    uY = move_vectors(2,:);
    
    
    
    % ----------------- Simulate with ode -----------------
    
    Tspan = [s*Tstep (s+1)*Tstep];
    IC = [X Y Z]';
    u = [uX uY uZ]';
    [T, XYZ] = ode45(@(t,y) DYNAMICS_simple(t, y, u), Tspan, IC);
    
    % Check if the movement kept the nodes inside omega
    for i=1:N
        inOmega = inpolygon( XYZ(end, i), XYZ(end, N+i), Xb, Yb);
        if inOmega
            % We want the last row of XYZ
            X = XYZ(end, 1:N );
            Y = XYZ(end, N+1:2*N );
            Z = XYZ(end, 2*N+1:3*N );
        end
        % Else keep the previous position
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Move in omega can be used here
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	end
    
        
    
    
	if PLOT_STATE_3D || PLOT_STATE_2D || PLOT_STATE_QUALITY
		% ----------------- Plot network 2D -----------------
		if PLOT_STATE_2D
			clf
			hold on
			% Region
			plot_poly( region, 'k');
			% Sensing disks and cells
			for i=1:N
				plot_poly( C{i}, 'r--');
				plot_poly( W{i}, 'k');
			end
			% Node positions
            for i=1:N
%                 tmpc = [Xs(s,i) + disk_rad * cos(t) ; Ys(s,i) + disk_rad * sin(t)];
%                 fill( tmpc(1,:), tmpc(2,:), 'k', 'EdgeColor', 'none' );
                plot( Xs(s,i), Ys(s,i), 'k.' )
                hold on
            end
            plot_AABB([-0.5 3 -0.5 3], 'w.');
			
			set( gca, 'Units', 'normalized', 'Position', [0 0 1 1] );
			axis([-0.5 3 -0.5 3])
			axis equal
			axis off
			
			if SAVE_PLOTS
				fname = strcat( '~/Frames/', sprintf('2D_frame_%d.png', s) );
% 				print(fname, '-dpng');
                saveas(gcf, fname);
			else
				pause(0.01);
			end
		end

		% ----------------- Plot network 3D -----------------
		if PLOT_STATE_3D
			clf
			hold on
			% Sensing disks and cells
			for i=1:N
				plot3_poly( [C{i} ; zeros(size(C{i}(1,:)))], 'r--');
				plot3_poly( [W{i} ; zeros(size(W{i}(1,:)))], 'k');
			end
			% Node positions and cones
			for i=1:N
				plot3( Xs(s,i), Ys(s,i), Zs(s,i), 'ko' )
				plot3( [Xs(s,i) Xs(s,i)], [Ys(s,i) Ys(s,i)], [Zs(s,i) 0], 'k--' )
				for j=1:24:PPC
					plot3([C{i}(1,j) Xs(s,i)], [C{i}(2,j) Ys(s,i)], [0 Zs(s,i)], 'r--');
				end
			end
			% Plot region
			plot3_poly( [region ; zeros(size(region(1,:)))], 'k' );
%             plot3(3, 3, zmax, 'w');
            plot3_AABB([-0.5 3 -0.5 3 0 zmax], 'w.');
			
			set( gca, 'Units', 'normalized', 'Position', [0 0 1 1] );
			view(-16, 34);
			axis([-0.5 3 -0.5 3 0 zmax])
			axis equal
			axis off
			
			if SAVE_PLOTS
				fname = strcat( '~/Frames/', sprintf('3D_frame_%d.png', s) );
% 				print(fname, '-dpng');
                saveas(gcf, fname);
			else
				pause(0.01);
			end
		end
	end % End of plot if
end
elapsed_time = toc;
average_iteration = elapsed_time / smax;
fprintf('\nSimulation time: %.4f s\n', elapsed_time)
fprintf('Average iteration time: %.4f s\n', average_iteration)




% ----------------- Final plots -----------------
% Plot covered area
figure;
plot( Tstep*linspace(1,smax,smax), 100*cov_area, 'b');
hold on
area_opt = 100 * N * pi * (zopt * tan(a))^2 / region_area;
plot( Tstep*[1 smax], [area_opt area_opt], 'k--');
axis([0 Tstep*smax 0 100]);
% axis([0 Tstep*smax 0 140]);
h = xlabel('$Time ~(s)$');
set(h,'Interpreter','latex')
h = ylabel('$A_{cov}~(\%)$');
title('Covered Area')
set(h,'Interpreter','latex')



% Plot objective
% Create grid
minxWi = min( C{i}(1,:) );
maxxWi = max( C{i}(1,:) );
minyWi = min( C{i}(2,:) );
maxyWi = max( C{i}(2,:) );
mingrid = min( [minxWi minyWi] );
maxgrid = max( [maxxWi maxyWi] );
lx = linspace(mingrid, maxgrid, gridsize);
ly = linspace(mingrid, maxgrid, gridsize);
[gx, gy] = meshgrid( lx, ly );
dx = lx(2)-lx(1);

% Start integration
H_opt = 0;
for k=1:gridsize^2
    if inpolygon( gx(k), gy(k), C{i}(1,:), C{i}(2,:) )
        H_opt = H_opt + dx^2 * fp(gx(k), gy(k), X(i), Y(i), Z(i), zmin, zmax, a, b);
    end
end
H_opt = N*H_opt;
figure;
plot( Tstep*linspace(1,smax,smax), 100*H / H_opt, 'b');
hold on
plot( Tstep*[1 smax], [100 100], 'k--');
% axis([0 Tstep*smax 0 100]);
h = xlabel('$Time ~(s)$');
set(h,'Interpreter','latex')
h = ylabel('$\frac{\mathcal{H}}{\mathcal{H}_{opt}} ~(\%)$');
set(h,'Interpreter','latex')


% Plot trajectories
traj = zeros(3,smax,N);
traj(1,:,:) = Xs;
traj(2,:,:) = Ys;
traj(3,:,:) = Zs;
figure;
plot3( [Xb Xb(1)], [Yb Yb(1)], zeros(length(Xb)+1), 'k');
hold on
for i=1:N
    % Initial
    plot3( traj(1,1,i), traj(2,1,i), traj(3,1,i), 'bo');
    hold on
    % Trajectory
    plot3( traj(1,:,i), traj(2,:,i), traj(3,:,i), 'b');
    hold on
    % Final
    plot3( traj(1,smax,i), traj(2,smax,i), traj(3,smax,i), 'bs');
    hold on
    
    % Projections
    % Initial
    plot3( traj(1,1,i), traj(2,1,i), 0, 'ko');
    hold on
    % Trajectory projection
    plot3( traj(1,:,i), traj(2,:,i), zeros(size(traj(3,:,i))), 'k--');
    hold on
    % Final
    plot3( traj(1,smax,i), traj(2,smax,i), 0, 'ks');
    hold on
end
axis equal
% axis([0 3 0 2.5 0 zmax]);
xlabel('x');
ylabel('y');
zlabel('z');

% ------------------- Save Results -------------------------
filename = ...
	strcat( 'results_decreasing_' , datestr(clock,'yyyymmdd_HHMM') , '.mat' );
save(filename);
