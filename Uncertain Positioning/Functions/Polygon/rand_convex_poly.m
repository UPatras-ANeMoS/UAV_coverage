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

% returns a random convex polygon with M vertices in xmin xmax ymin ymax
function [poly, mbc_center, mbc_radius] = rand_convex_poly(region, M, scale)

% Maximum and minimum values of region multiplied with scaling factor
% used to generate a polygon inside the region
xminsc = scale*min( region(1,:) );
xmaxsc = scale*max( region(1,:) );
yminsc = scale*min( region(2,:) );
ymaxsc = scale*max( region(2,:) );

% generate three random points
poly = [ xminsc+xmaxsc*rand([1 3]) ;
        yminsc+ymaxsc*rand([1 3]) ];

% add points until the polygon has M vertices
while length( poly(1,:) ) < M
    % add a random point
    new_pt = [ xminsc+xmaxsc*rand ; yminsc+ymaxsc*rand ];
    poly_temp = [poly  new_pt];
    
    % find the convex hull
    k = convhull(poly_temp(1,:), poly_temp(2,:));
    k = k(1:end-1);
    
    if length(k) == length(poly_temp(1,:))
        poly = [poly_temp(1,k) ; poly_temp(2,k)];
    end
    
end

% make the vertices clockwise
[poly(1,:), poly(2,:)] = poly2cw(poly(1,:), poly(2,:));


% find the minimum bound circle
[mbc_center, mbc_radius] = minboundcircle(poly);
% center is returned as row vector
mbc_center = mbc_center';

% offset the region by mbc_radius
region_offset = offset(region, mbc_radius, 'in');

% randomly select a new point in region_offest and move mbc_center (and poly) to it
while true

    % Maximum and minimum values of region
    % used to generate points in the region
    xmin = min( region_offset(1,:) );
    xmax = max( region_offset(1,:) );
    ymin = min( region_offset(2,:) );
    ymax = max( region_offset(2,:) );

    xnew = max([ xmax-xmin ymax-ymin ]) * rand(2,1) + min([ xmin ymin ]);
    if inpolygon( xnew(1), xnew(2), region_offset(1,:), region_offset(2,:) )
        break;
    end
end

% move the polygon
poly = bsxfun(@plus, poly, xnew-mbc_center);

% move mbc_center
mbc_center = xnew;

