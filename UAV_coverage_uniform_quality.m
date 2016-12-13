% MIT License
% 
% Copyright (c) 2016 Sotiris Papatheodorou
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT LOCATION ON CANVAS CHANGES ON EACH ITERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables
close all

PLOT_STATE_3D = 0;
PLOT_STATE_2D = 1;
PLOT_STATE_QUALITY = 0;
SAVE_PLOTS = 0;

% Bullo region
Xb=[ 0, 2.125, 2.9325, 2.975, 2.9325, 2.295, 0.85, 0.17 ];
Yb=[ 0, 0, 1.5, 1.6, 1.7, 2.1, 2.3, 1.2 ];
[Xb, Yb] = poly2cw(Xb, Yb);
region = [Xb ; Yb];
region_area = polyarea( Xb, Yb );

% ----------------- Network parameters -----------------
% Altitude constraints
zmin = 0.5;
zmax = 2.5;

% Optimal altitude
zopt = (2*zmin)/3 + (3*zmax^2 - 6*zmax*zmin + 4*zmin^2)^(1/2)/3;

% Quality-coverage tradeoff
Q = 1;

% Sensing cone angle (half the angle of the cone)
a = 20*pi/180;



% Initial positions - 3 nodes - 15 seconds - zopt = 1.359
% zoffset = 0.15;
% X = [0.40, 0.60, 0.55];
% Y = [0.50, 0.60, 0.50];
% Z = [0.45, 0.55, 0.50] + zoffset;

% Initial positions - 9 nodes - 20 seconds
% zoffset = 0.15;
% X = [0.40, 0.60, 0.55, 0.60, 0.50, 0.70, 0.60, 0.90, 0.80];
% Y = [0.50, 0.60, 0.50, 0.40, 0.60, 0.50, 0.75, 0.85, 0.95];
% Z = [0.45, 0.55, 0.50, 0.60, 0.40, 0.52, 0.57, 0.63, 0.65] + zoffset;

% Initial positions - 9 nodes - 5+ seconds
% X = [1, 1.2, 1.4, 1.6, 1.8, 0.3, 2, 2.1, 0.6];
% Y = [1, 1.5, 1.2, 0.8, 1.1, 0.5, 1.5, 1.3, 0.4];
% Z = [0.4, 0.6, 0.7, 0.8, 1, 0.38, 0.35, 0.5, 0.55];
% Z = [0.4, 0.6, 0.7, 0.8, 1, 0.55, 0.35, 0.5, 0.55]; % fi = fj

% Initial positions - Sensing in sensing - movement
% X = [1, 1.1];
% Y = [0.8, 0.85];
% Z = [0.38, 0.75];
% zmax = 1.3;

% Initial positions - Sensing in sensing - no movement
% X = [1, 1.05];
% Y = [0.8, 0.8];
% Z = [0.38, 1.1];
% zmax = 1.3;

% Initial positions - Same zi
X = [1, 1.05];
Y = [0.8, 0.8];
Z = [0.6, 0.6];

% Initial position at zmax
% X = [0.6];
% Y = [0.1];
% Z = [zmax-0.01];

% All cases
% X = [0.4, 0.7, 1.7, 1.8, 1.2, 1.7, 1.8];
% Y = [0.5, 0.4, 0.5, 0.4, 1.5, 1.9, 1.4];
% Z = [0.7, 0.7, 1.3, 0.6, 1.3, 1.0, 0.7];

% Almost coincident nodes
% X = [1, 1] + 0.2;
% Y = [1, 1];
% Z = [zmax-0.001, zmax-0.002];

% Empty cell
% X = [2.2, 2, 2.1, 2.4, 2.3] - 1;
% Y = [0.5, 0.35, 0.6, 0.6, 0.3] + 0.5;
% Z = [1.3, 0.9, 1.2, 1, 0.95];

% Number of nodes
N = length(X);

% Sensing radii
R = tan(a) * Z;


% ----------------- Simulation parameters -----------------
% Simulation steps
smax = 200;
% Simulation time step
Tstep = 0.1;
% Points Per Circle
PPC = 60;
% Radius for disks on plots
disk_rad = 0.02;
% vector for circle parameterization
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
% Coverage quality
f = zeros(1, N);
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
    
    % Coverage quality and sensing disks
    for i=1:N
        % Coverage quality, saturate to 0,1 for z outside (zmin,zmax)
        if Z(i) <= zmin
            f(i) = 1;
        elseif Z(i) >= zmax
            f(i) = 0;
        else
            f(i) = ( (Z(i)-zmin)^2 - (zmax-zmin)^2 )^2 / (zmax-zmin)^4;
        end
        
        % Sensing disks
        C{i} = [X(i) + R(i) * cos(t) ; Y(i) + R(i) * sin(t)];
    end
    
    % Store simulation data
    Xs(s,:) = X;
    Ys(s,:) = Y;
    Zs(s,:) = Z;
    Rs(s,:) = R;
    
    % Sensed space partitioning
    [W, overlap] = sensed_partitioning_uniform(Xb, Yb, C, f);
    
    % Find covered area and H objective
    for i=1:N
        if ~isempty(W{i})
            cov_area(s) = cov_area(s) + polyarea_nan(W{i}(1,:), W{i}(2,:));
            H(s) = H(s) + f(i) * polyarea_nan(W{i}(1,:), W{i}(2,:));
        end
    end
    cov_area(s) = cov_area(s)/region_area;
    
    
    
    % ----------------- Control law -----------------
    move_vectors = zeros(2, N); % [uX ; uY]
    uZ = zeros(1,N);
    for i=1:N
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
                [~, onB1] = inpolygon( pt1(1), pt1(2), Xb, Yb );
                [~, onB2] = inpolygon( pt2(1), pt2(2), Xb, Yb );

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
                            if overlap(i,j)
                                [inCj1] = inpolygon( pt1(1), pt1(2), C{j}(1,:), C{j}(2,:) );
                                [inCj2] = inpolygon( pt2(1), pt2(2), C{j}(1,:), C{j}(2,:) );

                                if inCj1 && inCj2
                                    % This is a dominant arc, normal vector with
                                    % magnitude fi-fj

                                    d = norm( [pt1(1)-pt2(1) , pt1(2)-pt2(2)] );
                                    n1 = (pt1-[X(i) ; Y(i)]) / R(i);
                                    n2 = (pt2-[X(i) ; Y(i)]) / R(i);
                                    nvector = (n1 + n2) / 2;

                                    % X-Y control law
                                    move_vectors(:,i) = move_vectors(:,i) + (f(i)-f(j)) * d * nvector;

                                    % Z control law
                                    uZ(i) = uZ(i) + (f(i)-f(j))*tan(a)*d;

                                    % DEBUG PLOTS
    %                                     plot( (pt1(1)+pt2(1))/2, (pt1(2)+pt2(2))/2, 'b.');
    %                                     hold on
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
                            n1 = (pt1-[X(i) ; Y(i)]) / R(i);
                            n2 = (pt2-[X(i) ; Y(i)]) / R(i);
                            nvector = (n1 + n2) / 2;

                            % X-Y control law
                            move_vectors(:,i) = move_vectors(:,i) + f(i) * d * nvector;

                            % Z control law
                            uZ(i) = uZ(i) + f(i)*tan(a)*d;

                            % DEBUG PLOTS
    %                             plot( (pt1(1)+pt2(1))/2, (pt1(2)+pt2(2))/2, 'g.');
    %                             hold on
                        end

                    end
                end
            end % line segment for


            % Area integral for Z control law
            % Find derivative of fi
            if Z(i)<=zmin || Z(i)>=zmax
                dfi=0;
            else
                dfi=4*(Z(i)-zmin)*((Z(i)-zmin)^2-(zmax-zmin)^2)/(zmax-zmin)^4;
            end

            uZ(i) = uZ(i) + Q * dfi * polyarea_nan( W{i}(1,:), W{i}(2,:) );
        end
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
		
		% ----------------- Plot network quality -----------------
		if PLOT_STATE_QUALITY
			clf
			hold on
			% Plot cylinders
			for i=1:N
				plot3_cell_quality(W{i}, f(i), 'r');
			end
			% Plot region
			plot3_poly( [region ; zeros(size(region(1,:)))], 'k' );
%             plot3(3, 3, zmax, 'w');
            plot3_AABB([-0.5 3 -0.5 3 0 zmax], 'w.');
			
			set( gca, 'Units', 'normalized', 'Position', [0 0 1 1] );
			view(-16, 34);
			axis([-0.5 3 -0.5 3 0 1])
			axis equal
			axis off
			
			if SAVE_PLOTS
				fname = strcat( '~/Frames/', sprintf('Q_frame_%d.png', s) );
% 				print(fname, '-dpng');
                saveas(gcf, fname);
			else
				pause(0.01);
			end
		end
	end

end
elapsed_time = toc
average_iteration = elapsed_time / smax




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
set(h,'Interpreter','latex')

% Plot objective
H_opt = N * pi * (zopt * tan(a))^2 * fi_u(zopt, zmin, zmax);
figure;
plot( Tstep*linspace(1,smax,smax), 100*H / H_opt, 'b');
hold on
plot( Tstep*[1 smax], [100 100], 'k--');
axis([0 Tstep*smax 0 100]);
h = xlabel('$Time ~(s)$');
set(h,'Interpreter','latex')
h = ylabel('$\frac{\mathcal{H}}{\mathcal{H}_{opt}} ~(\%)$');
set(h,'Interpreter','latex')

% Save trajectories
traj = zeros(3,smax,N);
traj(1,:,:) = Xs;
traj(2,:,:) = Ys;
traj(3,:,:) = Zs;

% ------------------- Save Results -------------------------
save('uniform_results.mat');
