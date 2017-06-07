% Copyright 2016-2017 Sotiris Papatheodorou
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

% Assume the origin is the center of rotation of the sensing pattern C
% Find the greatest distance from the origin to C
function [dmin, dmax] = sensing_pattern_radii(C)
dmin = Inf;
dmax = 0;
for i=1:length(C(1,:))
	if norm(C(:,i)) > dmax
		dmax = norm(C(:,i));
	end
	if norm(C(:,i)) < dmin
		dmin = norm(C(:,i));
	end
end
