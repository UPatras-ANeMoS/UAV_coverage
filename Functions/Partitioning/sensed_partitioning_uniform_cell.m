% Copyright 2017 Sotiris Papatheodorou
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

function W = sensed_partitioning_uniform_cell(region, C, f, i)

N = length(f);


% Initialize cell to sensing disk
W = C{i};


% Loop over all other nodes
for j=1:N
    if j ~= i
        % Remove a portion of the sensing region if fi < fj
        if f(i) < f(j)
            % Remove Cj
            [pbx, pby] = polybool( 'minus', W(1,:), W(2,:),...
            C{j}(1,:), C{j}(2,:) );

            % Save to the current cell
            W = [pbx ; pby];

        % Degenerate case, split the common region at the
        % intersection points
        elseif f(i) == f(j)
            % Find intersection points
            [xx, yy, ii] = polyxpoly( W(1,:), W(2,:),...
                        C{j}(1,:), C{j}(2,:) );
            % only if the sensing regions intersect
            if ~isempty(xx)
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
                [pbx, pby] = polybool( 'and', W(1,:), W(2,:),...
                    Ci(1,:), Ci(2,:) );
                W = [pbx ; pby];
            end
        end
    end
end

% AND with the region omega
if ~isempty(W)
    [pbx, pby] = polybool( 'and', W(1,:), W(2,:), region(1,:), region(2,:) );
    W = [pbx ; pby];
end