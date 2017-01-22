% Centroid of polygon P as a column vector
function C = centroid(P)
% The first row of P are the x values and the second the y values
% of the polygon vertices. The vertices must be CW or CCW.

% The first vertex must be the same as the last.
if P(:,1) ~= P(:,end)
    P = [P P(:,1)];
end
    
% http://es.mathworks.com/matlabcentral/answers/55056-i-m-trying-to-find-the-area-and-the-center-of-a-polygon
x = P(1,:);
y = P(2,:);
A = x(1:end-1).*y(2:end)-x(2:end).*y(1:end-1);
x_bar = (sum((x(2:end)+x(1:end-1)).*A)*1/6)/(sum(A)/2);
y_bar = (sum((y(2:end)+y(1:end-1)).*A)*1/6)/(sum(A)/2);

C = [x_bar ; y_bar];

