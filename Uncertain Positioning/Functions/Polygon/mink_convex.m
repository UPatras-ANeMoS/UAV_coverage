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

function M = mink_convex( A, B )


NvA = length(A(1,:));
NvB = length(B(1,:));

M = zeros(2,NvA*NvB);
k=1;

% Find the sums of all vertices of A with all vertices of B
for i=1:NvA
   for j=1:NvB
       M(:,k) = A(:,i) + B(:,j);
       k = k + 1;
   end
end

Mind = convhull( M(1,:), M(2,:) );

M = [M(1,Mind) ; M(2,Mind)];
