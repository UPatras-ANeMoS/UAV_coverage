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

% Plots the state of the network in the input struct
function plot_sim( sim )

N = length(sim.x(1,:));

hold on

% Plot cells
if sim.plot_cells
    for i=1:N
        plot_poly( sim.cells{i}, 'b');
    end
end

% Plot r limited cells
if sim.plot_rcells
    for i=1:N
        plot_poly( sim.rcells{i}, 'g');
    end
end

% Plot the region
plot_poly( sim.region, 'k');

% Plot nodes and uncertainty
plot(sim.x(1,:), sim.x(2,:), 'k.');
for i=1:N
    plot_circle( sim.x(1:2,i) , sim.uradii(i) , 'k');
end

% Plot sensing
for i=1:N
    plot_circle( sim.x(1:2,i) , sim.sradii(i) , 'r');
end

% Plot communication graph and radii
if sim.plot_comm
    for i=1:N
        plot_circle( sim.x(:,i) , sim.cradii(i) , 'r--');
    end
    % Plot communicating pairs
    for i=1:N
        for j=1:length(sim.in_range{i})
            if j ~= i
                % Only plot neighbors with higher index to avoid repetitions
                plot([sim.x(1,i) sim.x(1,sim.in_range{i}(j))], ...
                    [sim.x(2,i) sim.x(2,sim.in_range{i}(j))], 'g');
            end
        end
    end
end

% Plot velocity vectors
if sim.plot_vel
	for i=1:N
		plot_poly([sim.x(:,i) sim.x(:,i) + sim.velocity(:,i)], 'm');
	end
end

axis(sim.axis)
axis equal
axis off
