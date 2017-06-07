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

function [region, rdiameter, rarea] = read_region( filename )
% Read region from file
region = importdata( filename );
% File contains one vertex on each line, need to get the transpose
region = region';
% Make region clockwise
[xr , yr] = poly2cw(region(1,:), region(2,:));
region = [xr ; yr];
% Calculate diameter and area
rdiameter = diameter(region);
rarea = polyarea_nan(region(1,:), region(2,:));
