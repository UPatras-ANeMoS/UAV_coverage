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

% Plots the GV with guaranteed cells
function plot_GV(region_vert, x, GVcells, uncert, guarsrad, rcells )
%%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
hold on
%%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%

N = length( x(1,:) );

% Plot guaranteed cells
for i=1:N
    if ~isempty(GVcells{i})
        plot_poly( GVcells{i}, 'b');
%         fill( GVcells{i}(1,:), GVcells{i}(2,:) , 'b', 'EdgeColor','None');
        %fill_nan( GVcells{i}(1,:), GVcells{i}(2,:) , 'b');
        hold on;
    end
end


if nargin > 5
    % Plot covered cells
    for i=1:N
        if ~isempty(rcells{i})
%             fill( rcells{i}(1,:), rcells{i}(2,:) , 'g', 'EdgeColor','None');
            %fill_nan( rcells{i}(1,:), rcells{i}(2,:) , 'g');
        end
    end
end

% Plot region - Black outline
plot_poly( region_vert, 'k');

% plot nodes
scatter( x(1,:) , x(2,:) , '.k');
hold on;

% plot node uncertainty
for i=1:N
    plot_circle( x(1:2,i) , uncert(i) , 'k');
end

if nargin > 4
    % Plot guaranteed node radius
    for i=1:N
        plot_circle( x(1:2,i) , guarsrad(i), 'r' );
    end

%     % plot centroids
%     for i=1:N
%         if ~isempty( GVcells{i} )
%             centroids = centroid( GVcells{i} );
%             scatter( centroids(1,:), centroids(2,:), 'm');
%         end
%     end
end



axis square;
% axis([ -0.5 3.5 -0.5 3.5 ]);
% hold off;

