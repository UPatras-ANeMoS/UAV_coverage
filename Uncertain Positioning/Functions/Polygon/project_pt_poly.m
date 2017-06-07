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

% If x is outside polygon P, it projects it on its boundary and returns the
% projection. 
function pr = project_pt_poly( x, P )

% If x is in P do nothing
if inpolygon( x(1), x(2), P(1,:), P(2,:) )
    pr = x;
    return
end

% Add the first vertex to the end if needed
if ~isequal(P(:,1), P(:,end))
    P = [P P(:,1)];
end
N = length(P(1,:));

% Find the vertex of poly that is closest to x
minvd = inf;
for i=1:N-1
    if norm(x-P(:,i)) < minvd
        minvd = norm(x-P(:,i));
        minv = P(:,i);
    end
end

% Find the point on an edge of P that is closest to x
mine = [0;0];
mined = inf;
for i=1:N-1
    line = P(:,i:i+1);
    [pr, on] = project_pt_line( x, line );
    if on && norm(x-pr) < mined
       mined = norm(x-pr);
       mine = pr;
    end
end

% Return the closest point
if minvd < mined
    pr = minv;
else
    pr = mine;
end
