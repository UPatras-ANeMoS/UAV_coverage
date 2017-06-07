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

% Uniform coverage quality function derivative with respect to z
% Return a vector with the coverage quality for each altitude in z
function f = dfu(z, zmin, zmax)
    
    f = zeros(size(z));
    
    for i=1:length(z)
        if z(i) < zmin || z(i) > zmax
            f(i) = 0;
        else
            f(i) = (4*(z(i) - zmin)*((z(i) - zmin)^2 - (zmax - zmin)^2))/(zmax - zmin)^4;
        end
    end
