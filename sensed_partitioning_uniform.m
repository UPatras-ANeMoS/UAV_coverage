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

function [W, overlap] = sensed_partitioning_uniform(Xregion, Yregion, C, f)

N = length(f);

W = cell([1 N]);
overlap = zeros(N, N);

for i=1:N
	% Initialize cells to sensing disks
	W{i} = C{i};
end

for i=1:N
	% Loop over all other nodes
	for j=1:N
		if j ~= i
			% Remove a portion of the sensing disk if fi < fj
			if f(i) < f(j)
				% Overlap check
				[pbx, ~] = polybool( 'and', C{i}(1,:), C{i}(2,:),...
				C{j}(1,:), C{j}(2,:) );
				% Only change Wi on overlap
				if ~isempty(pbx)
					% set the corresponding elements of the overlap matrix
					overlap(i,j) = 1;
					overlap(j,i) = 1;

					% Remove Cj only if there is overlap
					[pbx, pby] = polybool( 'minus', W{i}(1,:), W{i}(2,:),...
					C{j}(1,:), C{j}(2,:) );

					% Save to the current cell
					W{i} = [pbx ; pby];
				end

			% Degenerate case, split the common region at the
			% intersection points
			elseif f(i) == f(j)
				% Find intersection points
				[xx, yy, ii] = polyxpoly( C{i}(1,:), C{i}(2,:),...
							C{j}(1,:), C{j}(2,:) );
				% only if the sensing circles intersect
				if ~isempty(xx)
					% set the corresponding elements of the overlap matrix
					overlap(i,j) = 1;
					overlap(j,i) = 1;
					% Add the intersection points
					Ci = C{i};
					for k=1:length(xx)
						% assuming ii is sorted
						Ci = [Ci(:,1:ii(k)+(k-1)) [xx(k) ; yy(k)] ...
							Ci(:,ii(k)+k:end)];
					end
					% Remove the points of Ci that are inside Cj
					[inCj, onCj] = inpolygon( Ci(1,:), Ci(2,:),...
						C{j}(1,:), C{j}(2,:) );
					Ci = Ci(:, ~inCj | onCj);
					% PolYregionool the new Ci with Wi
					[pbx, pby] = polybool( 'and', W{i}(1,:), W{i}(2,:),...
						Ci(1,:), Ci(2,:) );
					W{i} = [pbx ; pby];
				end
			end
		end
	end

	% AND with the region omega
	if ~isempty(W{i})
		[pbx, pby] = polybool( 'and', W{i}(1,:), W{i}(2,:), Xregion, Yregion );
		W{i} = [pbx ; pby];
	end
end