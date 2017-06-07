function [x3, y3] = polybooland(x1, y1, x2, y2)
% polybool AND

% Apply polygon set operations in the case where x1, y1, x2, and y2 are
% numerical vectors, possibly containing NaN-separators to distinguish
% between different parts within multipart polygons. The inputs may be
% either row vectors or column vectors. x3 and y3 will be column vectors
% unless x1 is a row vector.


% checkxy(x1, y1, mfilename, 'X1', 'Y1', 2, 3)
% checkxy(x2, y2, mfilename, 'X2', 'Y2', 2, 3)

[x3, y3, emptyInput] = handleEmptyInputs(x1, y1, x2, y2);
if ~emptyInput
    p1 = vectorsToGPC(x1, y1, mfilename, 'X1', 'Y1');
    p2 = vectorsToGPC(x2, y2, mfilename, 'X2', 'Y2');
    
    p3 = gpcmex('int', p1, p2);
    
    [x3, y3] = vectorsFromGPC(p3);
    
    rowVectorInput = (size(x1,2) > 1);
    if ~isempty(x3) && rowVectorInput
        x3 = x3';
        y3 = y3';
    end
end

%-----------------------------------------------------------------






function [x3, y3, emptyInput] = handleEmptyInputs(x1, y1, x2, y2)
% Assuming that x1 and y1, and x2 and y2, are pairs of inputs having
% consistent sizes, return the appropriate values for x3 and y3 in the
% event that x1 and/or x2 are empty (or contain only NaN), and set
% emptyInput to true. Otherwise, set x3 and y3 to empty and set
% emptyInput to false. Operation has been validated and equals one of
% the following strings: 'int','union','xor','diff'.

% NaN-only arrays should behave the same way as empty arrays, so filter
% them up-front. Be careful, because all(isnan([])) evaluates to true.
% Also, be careful to preserve shape: return 1-by-0 given a row
% vector of NaN and a 0-by-1 given a column vector of NaN.
if  all(isnan(x1)) && ~isempty(x1)
    x1(1:end) = [];
    y1(1:end) = [];
end

if all(isnan(x2)) && ~isempty(x2)
    x2(1:end) = [];
    y2(1:end) = [];
end

if isempty(x2)
    % Intersection is empty, but preserve shape
    % by using x2 and y2 rather than [].
    x3 = x2;
    y3 = y2;
    emptyInput = true;
elseif isempty(x1)
    % Intersection or difference is empty, but preserve
    % shape by using x1 and y1 rather than [].
    x3 = x1;
    y3 = y1;        
    emptyInput = true;
else
    x3 = [];
    y3 = [];
    emptyInput = false;
end

%-----------------------------------------------------------------

function p = vectorsToGPC(x, y, varargin)
% Convert polygon representation from a pair of NaN-separated coordinate
% vectors to a structure array that can be passed to gpcmex.
%
%   p = vectorsToGPC(x, y)
%   p = vectorsToGPC(x, y, func_name, var_name_1, var_name_2)

% Copyright 2009-2011 The MathWorks, Inc.

% Original code
% Locate the vertex coordinates within the NaN-separated arrays.
[first,last] = internal.map.findFirstLastNonNan(x);

% Pre-allocate structure array.
p(1, numel(first)) = struct('x',[],'y',[],'ishole',[]);

% Loop over parts and assign fields.
% There are no holes in the polygons I use
% ishole = ~ispolycw(x,y);
for k = 1:numel(first);
    p(k).x = x(first(k):last(k));
    p(k).y = y(first(k):last(k));
%     p(k).ishole = ishole(k);
    p(k).ishole = false;
end

% If all three optional inputs have been provided, check to see if
% there's at least one clockwise ring and warn if there is not.
% if numel(varargin) == 3
%     warnNoClockwiseRing(p, varargin{:})
% end

% % Faster version, no NaN detection, does not work on disjoint polygons
% p(1, 1) = struct('x',[],'y',[],'ishole',[]);
% 
% p(1).x = x;
% p(1).y = y;
% p(1).ishole = false;

%-----------------------------------------------------------------------

function [x, y] = vectorsFromGPC(p)
% Convert a structure array of the sort returned from gpcmex to a pair of
% NaN-separated coordinate vectors representing a multipart polygon.

% Copyright 2009 The MathWorks, Inc.

if isempty(p)
    x = [];
    y = [];
else
    % Old code
    % Determine the indices at which the parts will start and end in the
    % output vectors.
    numberOfVerticesPerPart = arrayfun(@(s) numel(s.x), p);
    
    first = cumsum(1 + [0; numberOfVerticesPerPart(1:end-1)]);
    last  = cumsum( ...
        [numberOfVerticesPerPart(1); 1 + numberOfVerticesPerPart(2:end)]);
    
    % Initialize the outputs as column vectors fill with NaN. After the
    % vertex coordinates are assigned, the remaining NaNs will serve to
    % separate the parts. Terminating NaNs are not provided.
    x = NaN(last(end),1);
    y = NaN(last(end),1);
    
    % Copy the vertex coordinates into the output vectors for all parts.
    for k = 1:numel(p)
        x(first(k):last(k)) = p(k).x;
        y(first(k):last(k)) = p(k).y;
    end
    
%     % Set the vertex direction around holes to be counter-clockwise. (This
%     % could be performed in the preceding loop, but that would entail
%     % multiple calls to ispolycw. If we had a lightweight "isringcw"
%     % function, however, then we might best call that within the loop.)
%     reverseVertices = find(~xor(ispolycw(x,y), [p.ishole]'));
%     for t = 1:numel(reverseVertices)
%         k = reverseVertices(t);
%         x(first(k):last(k)) = p(k).x(end:-1:1);
%         y(first(k):last(k)) = p(k).y(end:-1:1);
%     end
    
%     % Faster version, no NaN detection, does not work on disjoint polygons
%     x = p(numel(p)).x;
%     y = p(numel(p)).y;
end
