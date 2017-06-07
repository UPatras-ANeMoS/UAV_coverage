% Copyright 2016-2017 Sotiris Papatheodorou
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

clear variables
close all
addpath( genpath('Functions') );
rng(1);

%%%%%%%%%%%%%%%%%%%%%%%% Set simulation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Control law options %%%%%%
% Uncomment to select control law
% CTRL_LAW = 'CELL_CENTROID';
% CTRL_LAW = 'RCELL_CENTROID';
% CTRL_LAW = 'FREE_ARCS';
CTRL_LAW = 'GV_COMPLETE';
% CTRL_LAW = 'GV_COMPROMISE';
% CTRL_LAW = 'AWGV_COMPLETE';

% Prevent robots from colliding when coming close
STOP_COLLISIONS = 1;
% Keep the nodes inside the region at all times
KEEP_IN_REGION = 1;

% Control law gain
a = 1;


%%%%%% Simulation options %%%%%%
% Show the network state in each iteration
PLOT_STATE = 1;
% Save the network state in each iteration
SAVE_FRAMES = 0;
% Save the results to file
SAVE_RESULTS = 0;

% Simulation duration in seconds
Tfinal = 3;
% Time step in seconds
Tstep = 0.01;


%%%%%% Region %%%%%%
% Load region
[region, rdiameter, rarea] = read_region( 'Input Files/region.txt' );
axis_scale = [ -0.5 3.5 -0.5 3.5 ];

% [region, rdiameter, rarea] = read_region( 'Input Files/region_sq.txt' );
% [region, rdiameter, rarea] = read_region( 'Input Files/region_pi.txt' );
% sc = 1;
% axis_scale = [-sc sc -sc sc];

%%%%%% Nodes %%%%%%
% Load nodes
% [x, sradii, uradii, cradii] = read_nodes( 'Input Files/2_nodes.txt' );
% [x, sradii, uradii, cradii] = read_nodes( 'Input Files/3_nodes.txt' );
% [x, sradii, uradii, cradii] = read_nodes( 'Input Files/4_nodes.txt' );
[x, sradii, uradii, cradii] = read_nodes( 'Input Files/10_nodes.txt' );
% [x, sradii, uradii, cradii] = read_nodes( 'Input Files/6_nodes_inv_tri.txt' );
% [x, sradii, uradii, cradii] = read_nodes( 'Input Files/8_nodes_heterogeneous.txt' );
N = length(sradii);

% Inverted triangle
% N = 6;
% x = zeros(2,N);
% bi = 0.07;
% bo = 0.14;
% t = [30 150 270 90 210 330];
% x(1,:) = cosd(t);
% x(2,:) = sind(t);
% x(:,1:3) = bi * x(:,1:3);
% x(:,4:6) = bo * x(:,4:6);


% Sensing radius - overwrite radius from file
% sensing_rad = 0.5;
% sradii = sensing_rad * ones(1,N);
% sradii = sensing_rad * rand([1,N]) + sensing_rad/10;

% Uncertainty radius - overwrite radius from file
% uncert_rad = 0.05;
% uradii = uncert_rad * ones(1,N);
% uradii = uncert_rad * rand([1,N]) + uncert_rad/10;

% Communication radius - overwrite radius from file
comm_rad = 1;
cradii = comm_rad * ones(1,N);

% Use a finite communication range
FINITE_COMM_RANGE = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%% No need to change anything below %%%%%%%%%%%%%%%%%%%%%%
% Set communication range to region diameter if not relevant
if ~FINITE_COMM_RANGE
    cradii = rdiameter * ones(1,N);
end
% Use guaranteed sensing radii if there is uncertainty
sradii = sradii - uradii;
% Whether the robots have uncertain positions
if sum(uradii) == 0
    UNCERT = 0;
else
    UNCERT = 1;
end
% Increase uncertainty radius for collision avoidance
ur_e = 0.00001;


% Initializations
smax = floor(Tfinal/Tstep);
s = 1;
move_vector = zeros(2,N);
x_new = zeros(2,N);
xstorage = zeros(2,N,smax);
covered_area = zeros(1,smax);
H = zeros(1,smax);
in_range = cell([1 N]);
cells = cell([1 N]);
rcells = cell([1 N]);

% Create sim struct that contains all information, used for plots
sim = struct;
sim.plot_cells = 1;
sim.plot_rcells = 0;
sim.plot_comm = FINITE_COMM_RANGE;
sim.plot_vel = 0;
sim.region = region;
sim.x = x;
sim.uradii = uradii;
sim.sradii = sradii;
sim.cradii = cradii;
sim.cells = cells;
sim.rcells = rcells;
sim.velocity = zeros(2,N);
sim.in_range = in_range;
sim.axis = axis_scale;


if PLOT_STATE
    figure;
    hold on
end

tic;
while s <= smax
    eta = (smax-s) * toc/s;
    fprintf('%d/%d  %.2f%%  %.0f seconds remaining\n', s, smax, 100*s/smax, eta)
    
    
	% Find the neighbors and cell of each node
	for i=1:N
		% Find the nodes in communication range of each node i
		in_range{i} = in_comms_range( x, i, cradii(i) );
		% Put node i first in the list of neighbors
		in_range{i} = [i in_range{i}];
		tmpx = x(:,in_range{i});

		% Find the cell of each node i based on its neighbors
		cells{i} = AWGV_cell(region, ...
            tmpx, uradii(in_range{i}), sradii(in_range{i}), 1, rdiameter);
% 		cells{i} = GV_cell( region, tmpx, uradii, 1 );
		rcells{i} = rad_cell( x(:,i) , cells{i} , sradii(i));
    end
    
    % Update sim struct
    sim.cells = cells;
    sim.rcells = rcells;
    sim.in_range = in_range;
	
	% Store values and covered area - objective function H
    xstorage(:,:,s) = x;
    covered_area(s) = total_covered_area(region, x, sradii);
	% If there is uncertainty, the objective H is not the covered area
    if UNCERT
		for i=1:N
            if ~isempty(rcells{i})
                H(s) = H(s) + polyarea_nan(rcells{i}(1,:), rcells{i}(2,:));
            end
		end
    else
        H(s) = covered_area(s);
	end
	
    
    
    % Control law
    move_vector = zeros(2,N);
    for i=1:N % parfor slower
		% Select control law
        switch CTRL_LAW
            case 'CELL_CENTROID'
                move_vector(:,i) = ...
                    a * centroid_law( x(:,i), cells{i} );
            case 'RCELL_CENTROID'
                move_vector(:,i) = ...
                    a * centroid_law( x(:,i), rcells{i} );
            case 'FREE_ARCS'
                move_vector(:,i) = ...
                    a * C_integral_law_num( x(:,i), cells{i}, sradii(i) );
            case 'GV_COMPLETE'
                move_vector(:,i) = ...
                    a * C_integral_law_num( x(:,i), cells{i}, sradii(i) );
                %%%%%%% ADD DELAUNAY CHECK HERE %%%%%%%
				for j=1:N
					if i ~= j
						% Integral over Hij and Hji
						move_vector(:,i) = move_vector(:,i) + ...
						a * H_integral_law...
						( x(:,i), uradii(i), cells{i}, sradii(i),...
						  x(:,j), uradii(j), cells{j}, sradii(j), true );
					end
				end
            case 'GV_COMPROMISE'
                move_vector(:,i) = ...
                    a * C_integral_law_num( x(:,i), cells{i}, sradii(i) );
                %%%%%%% ADD DELAUNAY CHECK HERE %%%%%%%
				for j=1:N
					if i ~= j
						% Integral over Hij
						move_vector(:,i) = move_vector(:,i) + ...
						a * H_integral_law...
						( x(:,i), uradii(i), cells{i}, sradii(i),...
						  x(:,j), uradii(j), cells{j}, sradii(j), false );
					end
                end
            case 'AWGV_COMPLETE'
                move_vector(:,i) = ...
                    a * C_integral_law_num( x(:,i), cells{i}, sradii(i) );
                for j=1:N
					if i ~= j
						% Integral over Hij and Hji
						move_vector(:,i) = move_vector(:,i) + ...
						a * AWGV_cell_integrals...
						( x(:,i), uradii(i), cells{i}, sradii(i),...
						  x(:,j), uradii(j), cells{j}, sradii(j), true );
					end
                end
        end
    end
    
    
    % Simulate with ode45
    Tspan = [s*Tstep (s+1)*Tstep];
    IC = [x(1,:) x(2,:)]';
    u = [move_vector(1,:) move_vector(2,:)]';
    [T, Xsim] = ode45(@(t,y) DYNAMICS_integrator(t, y, u), Tspan, IC);
	% Keep only the final state of the ODE simulation
    x(1,:) = Xsim(end , 1:N);
    x(2,:) = Xsim(end , N+1:end);

    
    % Make sure no nodes collide
    if STOP_COLLISIONS
        x_temp = x;
        for i=1:N
			%%%%%%%%%%%% uradii has been increased here %%%%%%%%%%%%
            x_temp(:,i) = no_collisions(region, x, i, uradii(i)+ur_e);
        end
        x = x_temp;
    end

    
    % Make sure the new position is inside the region
    if KEEP_IN_REGION
        for i=1:N
            x(:,i) = keep_in_region(region, x(:,i), uradii(i));
        end
	end
    
	% Update the sim velocity
	sim.velocity = move_vector;
%     sim.velocity = x - sim.x;
	
	
	
	% Plot network state
    if PLOT_STATE
        clf
        plot_sim( sim );
    end

    % Pause for plot
    if PLOT_STATE
		if SAVE_FRAMES
			plot_AABB(axis_scale, 'w.');
			set( gca, 'Units', 'normalized', 'Position', [0 0 1 1] );
			fname = sprintf('~/Frames/frame_%05d.png',s);
			saveas(gcf, fname);
		else
			pause(0.01)
		end
	end
	
	% Update the sim struct with the new positions
	sim.x = x;
    
    
    s = s + 1;
end

elapsed_time = toc;
average_iteration = elapsed_time/s;
fprintf('Simulation time %.3f s\n', elapsed_time)
fprintf('Average iteration %.3f s\n', average_iteration)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Final plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show final state if it was not shown during the simulation
if ~PLOT_STATE
    figure
    plot_sim( sim )
end

% Create time vector
t = Tstep*linspace(0,smax-1,smax);
% Total objective H
figure('name','Objective')
hold on
if UNCERT
    plot(t, 100 * H / rarea, 'b');
    plot(t, 100 * covered_area / rarea, 'r');
    legend('H', 'Covered area', 'Location','southeast')
else
    plot(t, 100 * H / rarea, 'b');
    legend('H - Covered area', 'Location','southeast')
end
xlabel('Time (s)')
ylabel('H - Covered area (% of total area)')



% Save results
if SAVE_RESULTS
	% Restore the original sensing radius
	sradii = sradii + uradii;
	
    % Generate filemane
    filename = strcat( 'sim_' , CTRL_LAW );
    if FINITE_COMM_RANGE
        filename = strcat( filename , '_fincomm' );
    end
    if KEEP_IN_REGION
        filename = strcat( filename , '_inreg' );
    end
    if STOP_COLLISIONS
        filename = strcat( filename , '_nocol' );
    end
    filename = ...
        strcat( filename , '_' , datestr(clock,'yyyymmdd_HHMM') , '.mat' );


    save(filename,'elapsed_time','smax','Tstep','region','N', 'xstorage', ...
		'uradii','sradii','cradii','H','covered_area','axis_scale', ...
        'CTRL_LAW', 'FINITE_COMM_RANGE', 'KEEP_IN_REGION', 'STOP_COLLISIONS');
end
