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

% Return the GV cell of node i with respect to all nodes in x
function C = GV_cell( region, x, uradii, i )

% Initialize the cell to the region
C = region;

% Loop over all points in x
for j=1:length(x(1,:))
    if j ~= i
        % Calculate the hyperbola branch
        [E, ~, ~] = GV_edge(x(:,i), x(:,j), uradii(i), uradii(j), 101, diameter(region));
        if ~isempty(E) && ~isempty(C)
            E = E(1:2,:);
            % Intersect with the current cell
            [ipolyx , ipolyy] =...
                    polybooland(C(1,:), C(2,:),...
                    E(1,:), E(2,:) );

            C = [ ipolyx; ipolyy];
        else
            C = [];
            return
        end
    end
end
