% MIT License
% 
% Copyright (c) 2016 Sotiris Papatheodorou
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% Finds the area of a polygon with holes
% polygon contours are Nan delimited
function area = polyarea_nan(x, y)
% The outer contour is the last NaN delimited segment

nan_indices = find(isnan( x ));
indices = [ 0 nan_indices length(x)+1 ];
area = 0;

if ~isempty(nan_indices)
    for i=1:length(nan_indices)+1

        tempx = x( indices(i)+1 : indices(i+1)-1 );
        tempy = y( indices(i)+1 : indices(i+1)-1 );
        if ispolycw(tempx, tempy)
            % external contour
            area = area + polyarea(tempx, tempy);
        else
            % internal contour
            area = area - polyarea(tempx, tempy);
        end

    end
else
    area = polyarea(x, y);
end

