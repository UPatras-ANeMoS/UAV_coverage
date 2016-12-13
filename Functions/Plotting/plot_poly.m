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

% Plot NaN delimited polygons
% The first line of P contains the x coordinates of the polygon's vertices
% and the second line their y cordinate
% 2016/3/4
function H = plot_poly( P , colorspec )

hh = ishold;
hold on

if nargin == 1
    colorspec = 'b';
end

if ~isempty(P)
    x = P(1,:);
    y = P(2,:);

    nan_indices = find( isnan(P(1,:)) ); % The indices are the same for x and y
    indices = [ 0 nan_indices length(x)+1 ];

    % Collect the handles of all the plots
    H = zeros(1, length( nan_indices )+1);
    for i=1:length( nan_indices )+1

        tmpP = [x( indices(i) + 1 : indices(i+1) - 1 ) ; ...
                y( indices(i) + 1 : indices(i+1) - 1 )];
        % if the last vertex is not the same as the last, add the first vertex to
        % the end
        if ~isequal( tmpP(:,1), tmpP(:,end) )
            tmpP = [tmpP tmpP(:,1)];
        end

        h = plot( tmpP(1,:), tmpP(2,:), colorspec );

        H(i) = h;

    end

    if ~hh
		hold off;
	end
end
