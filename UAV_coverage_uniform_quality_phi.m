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

%%%%%%%%%%%%%%%%%%% Set Simulation Options %%%%%%%%%%%%%%%%%%%
% Network options
% Altitude constraints
zmin = 0.3;
zmax = 2.3;
% Sensing cone angle (half the angle of the cone)
a = 20*pi/180;

% A priori importance function
% phi = @PHI_uniform;
phi = @PHI_gaussian1;
% phi = @PHI_gaussian2;

% Simulation options
% Simulation duration in seconds
Tfinal = 15;
% Time step in seconds
Tstep = 0.01;

% Control law options
% Planar control law gain
axy = 1;
% Altitude control law gain
az = 1;

% Network plots to show during simulation
PLOT_STATE_2D = 0;
PLOT_STATE_3D = 1;
PLOT_STATE_PHI = 1;
PLOT_STATE_QUALITY = 1;
SAVE_PLOTS = 1;

% Save simulation results to file
SAVE_RESULTS = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







% Add function path
addpath( genpath('Functions') );

% ---------------- Region ----------------
% Bullo region
Xb=[ 0, 2.125, 2.9325, 2.975, 2.9325, 2.295, 0.85, 0.17 ];
Yb=[ 0, 0, 1.5, 1.6, 1.7, 2.1, 2.3, 1.2 ];
[Xb, Yb] = poly2cw(Xb, Yb);
region = [Xb ; Yb];
region_area = polyarea( Xb, Yb );
axis_scale = [-0.5 3 -0.5 3];

% ---------------- Initial State ----------------
% Initial positions - 3 nodes - 15 seconds
X = [0.40, 0.60, 0.55];
Y = [0.50, 0.40, 0.50];
Z = [0.45, 0.55, 0.50];

% Initial positions - 3 nodes spread out - 5 seconds
% X = [0.40, 2.00, 0.85];
% Y = [0.50, 0.40, 2.00];
% Z = [0.45, 0.55, 0.50];


% Initial positions - 9 nodes - 20 seconds
% X = [0.40, 0.60, 0.55, 0.60, 0.50, 0.70, 0.60, 0.90, 0.80];
% Y = [0.50, 0.60, 0.50, 0.40, 0.60, 0.50, 0.75, 0.85, 0.95];
% Z = [0.45, 0.55, 0.50, 0.60, 0.40, 0.52, 0.57, 0.63, 0.65];

% Initial positions - 7 nodes (deleted two from 9 nodes) - 8 seconds
% X = [1.6910359269619974 1.0671555014553136 ...
% 	1.3060376187415315 0.50870489556525222 2.0464590369292974 ...
% 	1.0609262551606342 1.8494381205125359];
% Y = [1.2370376839041961 0.98786408275888737 ...
% 	0.33136168814863487 1.2225014356751873 0.57709034562636752 ...
% 	1.8083501035247116 1.8935528313332963];
% Z = [1.0966326999239728 1.0448629252962573 ...
% 	1.1043691751220159 1.04059568718053 1.1572027919849177 ...
% 	1.2338162957229133 1.022310013652235];

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
% X = [1, 1.05];
% Y = [0.8, 0.8];
% Z = [0.6, 0.6];

% Initial position at zmax
% X = [0.6];
% Y = [0.1];
% Z = [zmax-0.01];

% Almost coincident nodes
% X = [1, 1] + 0.2;
% Y = [1, 1];
% Z = [zmax-0.001, zmax-0.002];

% Empty cell
% X = [2.2, 2, 2.1, 2.4, 2.3] - 1;
% Y = [0.5, 0.35, 0.6, 0.6, 0.3] + 0.5;
% Z = [1.3, 0.9, 1.2, 1, 0.95];


% ---------------- Simulation initializations ----------------
% Caclulate optimal altitude
zopt = z_optimal_uniform(zmin, zmax);
% Number of nodes
N = length(X);
% Simulation steps
smax = floor(Tfinal/Tstep);
% Grid size for double integrals
gridsize = 50;
% Points Per Circle
PPC = 60;
% Radius for points on plots
disk_rad = 0.02;
% Vector for circle parameterization
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
% Control inputs
uX = zeros(1,N);
uY = zeros(1,N);
uZ = zeros(1,N);
% Communication range for each node
r_comm = zeros(1,N);
% Adjacency matrix for communication graph
A = zeros(N,N);

% Create sim struct that contains all information, used for plots
sim = struct;
sim.region = region;
sim.axis = axis_scale;
sim.phi = phi;
sim.PPC = PPC;
sim.zmin = zmin;
sim.zmax = zmax;
sim.X = X;
sim.Y = Y;
sim.Z = Z;
sim.N = N;
sim.C = C;
sim.W = W;
sim.f = f;
sim.A = A;
sim.PLOT_COMMS = 0;
sim.PLOT_STATE_3D = PLOT_STATE_3D;
sim.PLOT_STATE_2D = PLOT_STATE_2D;
sim.PLOT_STATE_PHI = PLOT_STATE_PHI;
sim.PLOT_STATE_QUALITY = PLOT_STATE_QUALITY;
sim.SAVE_PLOTS = SAVE_PLOTS;






%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%
if PLOT_STATE_3D || PLOT_STATE_2D || PLOT_STATE_QUALITY || PLOT_STATE_PHI
	figure
end
tic;
for s=1:smax
	fprintf('%.2f%% complete\n',100*s/smax);

    % ----------------- Partitioning -----------------
    % Sensing radii
    R = tan(a) * Z;
    % Coverage quality
    f = fu(Z, zmin, zmax);
    % Sensing disks
    for i=1:N
        C{i} = [X(i) + R(i) * cos(t) ; Y(i) + R(i) * sin(t)];
    end
    % Communication range
    r_comm = communication_range(Z, zmin, zmax, a);
    
    % Store simulation data
    Xs(s,:) = X;
    Ys(s,:) = Y;
    Zs(s,:) = Z;
    Rs(s,:) = R;
    
    % Sensed space partitioning
    for i=1:N
        % Find the nodes in communication range of each node i
		A(i,:) = in_comms_range3( X, Y, Z, i, r_comm(i) );
        
        % The index of i in the reduced state vector is
        ind = sum(A(i,1:i));
        
		% Find the cell of each node i based on its neighbors
		W{i} = sensed_partitioning_uniform_cell(region, ...
            C( logical(A(i,:)) ), f( logical(A(i,:)) ), ind);
    end
    
    
    % ----------------- Plots -----------------
    sim.X = X;
    sim.Y = Y;
    sim.Z = Z;
    sim.C = C;
    sim.W = W;
    sim.f = f;
    sim.A = A;
    sim.s = s;
    clf
    plot_sim_UAV( sim );
    
    
    % ----------------- Objective -----------------
    % Find covered area and H objective
    for i=1:N
        if ~isempty(W{i})
            % Numerically integrate phi on Wi
			I = 0;
			[xm, ym] = meshgrid(linspace(X(i)-R(i),X(i)+R(i),gridsize), ...
				linspace(Y(i)-R(i),Y(i)+R(i),gridsize));
			dx = abs(xm(1,1)-xm(1,2));
			dy = abs(ym(1,1)-ym(2,1));
            for l=1:gridsize^2
                if inpolygon(xm(l), ym(l), W{i}(1,:), W{i}(2,:))
                    ds = dx*dy;
                    I = I + ds * phi(xm(l), ym(l));
                end
            end
            H(s) = H(s) + f(i) * I;
            cov_area(s) = cov_area(s) + polyarea_nan(W{i}(1,:), W{i}(2,:));
        end
    end
    
    
    % ----------------- Control law -----------------
    for i=1:N
        % The index of i in the reduced state vector is
        ind = sum(A(i,1:i));
        
        % Give correct info based on adjacency matrix A
        [uX(i), uY(i), uZ(i)] = ...
            control_uniform_phi(region, phi, gridsize, zmin, zmax, a, ...
            W(logical(A(i,:))), C(logical(A(i,:))), f(logical(A(i,:))), ...
            i, X(i), Y(i), Z(i));
    end
    
    % Control inputs with gain
    uX = axy * uX;
    uY = axy * uY;
    uZ = az * uZ;
    
    
    % ----------------- Simulate with ode -----------------
    Tspan = [s*Tstep (s+1)*Tstep];
    IC = [X Y Z]';
    u = [uX uY uZ]';
    [T, XYZ] = ode45(@(t,y) DYNAMICS_simple(t, y, u), Tspan, IC);
    
    % Keep the last row of XYZ
    X = XYZ(end, 1:N );
    Y = XYZ(end, N+1:2*N );
    Z = XYZ(end, 2*N+1:3*N );
    
end
elapsed_time = toc;
average_iteration = elapsed_time / smax;
fprintf('\nSimulation time: %.4f s\n', elapsed_time)
fprintf('Average iteration time: %.4f s\n', average_iteration)




%%%%%%%%%%%%%%%%%%% Final plots %%%%%%%%%%%%%%%%%%%
% Plot covered area
figure;
plot( Tstep*linspace(1,smax,smax), 100*cov_area/region_area, 'b');
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
figure;
plot( Tstep*linspace(1,smax,smax), H, 'b');
axis([0 Tstep*smax 0 ceil(max(H))]);
h = xlabel('$Time ~(s)$');
set(h,'Interpreter','latex')
h = ylabel('$\mathcal{H}$');
set(h,'Interpreter','latex')

% Save trajectories
traj = zeros(3,smax,N);
traj(1,:,:) = Xs;
traj(2,:,:) = Ys;
traj(3,:,:) = Zs;

%%%%%%%%%%%%%%%%%%% Save Results %%%%%%%%%%%%%%%%%%%
if SAVE_RESULTS
    filename = ...
        strcat( 'results_uniform_phi_' , datestr(clock,'yyyymmdd_HHMM') , '.mat' );
    save(filename);
end
