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

function contour = keep_single_contour( P, x )
% Keep only the contour of P that x is in
contour = [];
% The number of contours is the number of NaNs + 1
Nc = sum(isnan(P(1,:)))+1;

nan_indices = find( isnan(P(1,:)) );
nan_indices = [0 nan_indices length(P)+1];

for i=1:Nc
   % Create contour
   C = P(:,nan_indices(i)+1:nan_indices(i+1)-1);
   % If x is inside the current contour, return it
   if inpolygon( x(1),x(2), C(1,:), C(2,:) )
       contour = C;
       return
   end
end
