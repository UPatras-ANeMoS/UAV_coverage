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

function rc = communication_range(Z, zmin, zmax, a)

% The three worst case scenarios
r1 = 2*Z*tan(a);
r2 = sqrt( (Z+zmin).^2 * tan(a)^2 + (Z-zmin).^2 );
r3 = sqrt( (Z+zmax).^2 * tan(a)^2 + (Z-zmax).^2 );

% Use tha largest distance
rc = max(r1, r2);
rc = max(rc, r3);
