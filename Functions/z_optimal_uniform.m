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

% Returns the optimal altitude for the given network parameters
% and uniform coverage quality (third argument is a placeholder for a)
function zopt = z_optimal_uniform(zmin, zmax)

zopt = (2*zmin)/3 + (3*zmax^2 - 6*zmax*zmin + 4*zmin^2)^(1/2)/3;
