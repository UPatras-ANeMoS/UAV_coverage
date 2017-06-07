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

function Poff = offset_out(P, r)
% homothetic transformation of polygon P by r wrt to it's centroid
% B = O + k(A-O)

centroidP = centroid(P);
Poff = r.*bsxfun(@minus, P, centroidP);
Poff = bsxfun(@plus, Poff, centroidP);
