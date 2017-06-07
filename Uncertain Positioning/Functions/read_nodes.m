% Copyright 2017 Sotiris Papatheodorou
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

function [q, sradii, uradii, cradii] = read_nodes( filename )

% Each line of the file contains a single node
% x y sensing_radius uncertainty_radius communication_radius

d = importdata( filename );
% Get number of nodes and number of fields per node
[N, Nf] = size(d);

% Get nodes
q = d(:,1:2)';
% Get sensing radii
if Nf >= 3
    sradii = d(:,3)';
else
    sradii = zeros(1,N);
end
% Get uncertainty radii
if Nf >= 4
    uradii = d(:,4)';
else
    uradii = zeros(1,N);
end
% Get communication radii
if Nf == 5
    cradii = d(:,5)';
else
    cradii = zeros(1,N);
end
