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

% Plots a circle with center O and radius r
function plot_handle = plot_circle( O , r , linespec)

if nargin < 3
    linespec = 'b';
end

t = linspace(0, 2*pi, 360);

x = r*cos(t) + O(1);
y = r*sin(t) + O(2);

plot_handle = plot(x,y,linespec);

