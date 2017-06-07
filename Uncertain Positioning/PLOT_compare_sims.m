% Copyright 2016-2017 Sotiris Papatheodorou
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

% Compare simulations for area coverage of a convex region by circular sensors
% Created 8/10/2016

clear variables
close all


%%%%%%%%%%%%%%%%%%%%%%%%%% Set plot settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulations to load
sim1 = 'sim_FREE_ARCS_inreg_20161209_0006.mat';
sim2 = 'sim_GV_COMPLETE_inreg_20161209_0004.mat';
% Simulation labels for plots
s1_label = 'Free arcs';
s2_label = 'Complete law';
% Simulation colors
s1_color = 'b';
s2_color = 'r';
% Plots to show
TRAJECTORIES = 1;
OBJECTIVE = 1;
AREA = 0;
OBJECTIVE_AREA = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Add function path
addpath( genpath('Common Functions') );

% Load simulations
s1 = load(sim1);
s2 = load(sim2);
s1.axis_scale = [-0.5 2.5 0 2.5];
s2.axis_scale = [-0.5 2.5 0 2.5];

%%%%%%%%%% Plot trajectories %%%%%%%%%%
if TRAJECTORIES
	figure('name','Trajectories')
	hold on
	
	% Simulation 1
	h1 = plot( squeeze( s1.xstorage(1,:,:) )', ...
		squeeze( s1.xstorage(2,:,:) )', s1_color);
	hi = plot( squeeze( s1.xstorage(1,:,1) )', ...
		squeeze( s1.xstorage(2,:,1) )', 'ko');
	plot( squeeze( s1.xstorage(1,:,end) )', ...
		squeeze( s1.xstorage(2,:,end) )', strcat(s1_color,'s'));
	plot_poly( s1.region, 'k' );
	
	% Simulation 2
	h2 = plot( squeeze( s2.xstorage(1,:,:) )', ...
		squeeze( s2.xstorage(2,:,:) )', s2_color);
	plot( squeeze( s2.xstorage(1,:,1) )', ...
		squeeze( s2.xstorage(2,:,1) )', 'ko');
	plot( squeeze( s2.xstorage(1,:,end) )', ...
		squeeze( s2.xstorage(2,:,end) )', strcat(s2_color,'s'));
	plot_poly( s2.region, 'k' );
	
	axis( max(s1.axis_scale, s2.axis_scale) );
    axis off
	axis equal
% 	legend([hi h1(1) h2(1)], ...
% 		'Initial positions', s1_label, s2_label, 'Location','southeast')
    
    % Final state distance
    mean_dist = 0;
    for i=1:s1.N
        mean_dist = mean_dist + norm(s1.xstorage(:,i,end) - s2.xstorage(:,i,end));
    end
    mean_dist = mean_dist / s1.N;
    mean_dist = 100 * mean_dist / diameter( s1.region );
    fprintf('Final state mean distance %.2f %% of region diameter\n',mean_dist);
end


%%%%%%%%%% Plot objective %%%%%%%%%%
if OBJECTIVE
	figure('name','H')
	hold on
	
	% Simulation 1
	t1 = linspace(0, s1.Tstep*s1.smax, s1.smax);
	Hmax1 = 0;
    for i=1:s1.N
        Hmax1 = Hmax1 + pi * (s1.sradii(i)-s1.uradii(i))^2;
    end
	plot(t1, 100*s1.H/Hmax1, s1_color);
	
	% Simulation 2
	t2 = linspace(0, s2.Tstep*s2.smax, s2.smax);
	Hmax2 = 0;
    for i=1:s1.N
        Hmax2 = Hmax2 + pi * (s1.sradii(i)-s1.uradii(i))^2;
    end
	plot(t2, 100*s2.H/Hmax2, s2_color);

	xlabel('Time (s)')
	ylabel('H (% of maximum)')
% 	legend(s1_label, s2_label, 'Location','southeast')
	axis([0 max(t1(end), t2(end)) 0 100])
end


%%%%%%%%%% Plot area %%%%%%%%%%
if AREA
	figure
	hold on
	
	% Simulation 1
	t1 = linspace(0, s1.Tstep*s1.smax, s1.smax);
	Amax1 = polyarea_nan(s1.region(1,:), s1.region(2,:));
	plot(t1, 100*s1.covered_area/Amax1, s1_color);
	
	% Simulation 2
	t2 = linspace(0, s2.Tstep*s2.smax, s2.smax);
	Amax2 = polyarea_nan(s2.region(1,:), s2.region(2,:));
	plot(t2, 100*s2.covered_area/Amax2, s2_color);

	xlabel('Time (s)')
	ylabel('Area (% of maximum)')
	title('Covered area')
	legend(s1_label, s2_label, 'Location','southeast')
	axis([0 max(t1(end), t2(end)) 0 100])
end


%%%%%%%%%% Plot objective and area %%%%%%%%%%
if OBJECTIVE_AREA
	figure
	hold on
	
	% Simulation 1
	t1 = linspace(0, s1.Tstep*s1.smax, s1.smax);
	Amax1 = polyarea_nan(s1.region(1,:), s1.region(2,:));
    Hmax1 = 0;
    for i=1:s1.N
        Hmax1 = Hmax1 + pi * (s1.sradii(i)-s1.uradii(i))^2;
    end
	plot(t1, 100*s1.H/Hmax1, s1_color);
	plot(t1, 100*s1.covered_area/Amax1, [s1_color '--']);
	
	% Simulation 2
	t2 = linspace(0, s2.Tstep*s2.smax, s2.smax);
	Amax2 = polyarea_nan(s2.region(1,:), s2.region(2,:));
    Hmax2 = 0;
    for i=1:s1.N
        Hmax2 = Hmax2 + pi * (s1.sradii(i)-s1.uradii(i))^2;
    end
	plot(t2, 100*s2.H/Hmax2, s2_color);
	plot(t2, 100*s2.covered_area/Amax2, [s2_color '--']);

	xlabel('Time (s)')
	ylabel('% of maximum')
	title('Covered area and coverage objective')
	legend([s1_label ' objective'], [s1_label ' area'],...
		[s2_label ' objective'], [s2_label ' area'],...
		'Location','southeast')
	axis([0 max(t1(end), t2(end)) 0 100])
end
