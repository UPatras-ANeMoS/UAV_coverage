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

% returns a vector with two differently linearly spaced portions
% fracl [0,1] defines how large the first portion will be compared to
% the total vector length
% fracp [0,1] defines what part of the total points the first portion will
% contain
function x = linspace2(a, b, N, fracl, fracp)

% point defining the two portions [a,c] [c,b]
c = a + fracl*(b-a);

% points in the first portion
N1 = ceil( fracp*N );
N2 = N-N1;

x1 = linspace(a, c, N1);
x2 = linspace(c, b, N2+1);
% N2+1 because point c on x2 will be removed afterwards

x = [x1 x2(2:end)];

