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

function move_vector = H_integral_law2( Xi, uradiusi, GVcelli, sradiusi, ...
    Xj, uradiusj, GVcellj, sradiusj, DEBUG )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IT DOESN'T WORK FOR MORE THAN TWO NODES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the integral over Hij, Hji and return the move vector
% Handle optional debug argument
if nargin == 8
    DEBUG = 0;
end
% Number of points
p = 501;

% Translate and rotate Xi and Xj so that the hyperbola is East-West
tvector = (Xi+Xj)/2;
XXi = Xi - tvector;
XXj = Xj - tvector;
theta = -atan2( Xj(2)-Xi(2), Xj(1)-Xi(1) );
XXi = rot(XXi, theta);
XXj = rot(XXj, theta);

% Also translate and rotate both cells
GVi = bsxfun(@minus, GVcelli, (Xi+Xj)/2);
GVj = bsxfun(@minus, GVcellj, (Xi+Xj)/2);
for k=1:length( GVi(1,:) )
    GVi(:,k) = rot( GVi(:,k), theta );
end
for k=1:length( GVj(1,:) )
    GVj(:,k) = rot( GVj(:,k), theta );
end

% Hyperbola parameters
a = (uradiusi + uradiusj) / 2;
c = norm( Xi-Xj ) / 2;
b = sqrt(c^2 - a^2);

move_vector = zeros(2,1);





% Integral over Hji
% Find the intersection points of the sensing circle and the hyperbola
% Create a linspace and integrate over them
% Create sensing circle
C = create_circle(XXi, sradiusi, 360);
[~, yint] = polyxpoly( C(1,:), C(2,:), GVi(1,:), GVi(2,:) );
if length(yint) == 2
    % Get t values of the intersection points
    t1 = asinh( yint(1)/b );
    t2 = asinh( yint(2)/b );
    tmin = min([t1 t2]);
    tmax = max([t1 t2]);

    % Create a vector with p t values between tmin and tmax
    tt = linspace(tmin,tmax,p);
    dt = abs(tt(1)-tt(2));
    for k=1:length(tt)
        % Substitute the symbolic expressions
        t = tt(k);
        xi = XXi(1);
        yi = XXi(2);
        xj = XXj(1);
        yj = XXj(2);
        uni = FJni(a,t,xi,xj,yi,yj);

        % Rotate back to the correct orientation
        uni = rot( uni, -theta );

        % Find the length of the arc (half the distance from each adjacent vertex)
        % Find the current, previous and next points on the hyperbola
        Hc = [-a*cosh(tt(k)) ; b*sinh(tt(k))];
        Hp = [-a*cosh(tt(k)-dt) ; b*sinh(tt(k)-dt)];
        Hn = [-a*cosh(tt(k)+dt) ; b*sinh(tt(k)+dt)];
        d = norm(Hc-Hp) / 2 + ...
            norm(Hc-Hn) / 2;
        % Add to the move vector
        move_vector = move_vector + d*uni;
        %%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
        if DEBUG
            Hc = rot(Hc, -theta)+tvector;
            plot(Hc(1), Hc(2), 'g.');
            plot([Hc(1) Hc(1)+d*uni(1)], ...
                 [Hc(2) Hc(2)+d*uni(2)], 'g');
        end
        %%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
    end
    
end




% Integral over Hji
% Find the intersection points of the sensing circle and the hyperbola
% Create a linspace and integrate over them
% Create sensing circle
% Create sensing circle
C = create_circle(XXj, sradiusj, 360);
[~, yint] = polyxpoly( C(1,:), C(2,:), GVj(1,:), GVj(2,:) );
if length(yint) == 2
    % Get t values of the intersection points
    t1 = asinh( yint(1)/b );
    t2 = asinh( yint(2)/b );
    tmin = min([t1 t2]);
    tmax = max([t1 t2]);
    
    % Create a vector with p t values between tmin and tmax
    tt = linspace(tmin,tmax,p);
    dt = abs(tt(1)-tt(2));
    for k=1:length(tt)
        % Substitute the symbolic expressions
        t = tt(k);
        xi = XXi(1);
        yi = XXi(2);
        xj = XXj(1);
        yj = XXj(2);
        unj = FJnj(a,t,xi,xj,yi,yj);

        % Rotate back to the correct orientation
        unj = rot( unj, -theta );

        % Find the length of the arc (half the distance from each adjacent vertex)
        % Find the current, previous and next points on the hyperbola
        Hc = [a*cosh(tt(k)) ; b*sinh(tt(k))];
        Hp = [a*cosh(tt(k)-dt) ; b*sinh(tt(k)-dt)];
        Hn = [a*cosh(tt(k)+dt) ; b*sinh(tt(k)+dt)];
        d = norm(Hc-Hp) / 2 + ...
            norm(Hc-Hn) / 2;
        % Add to the move vector
        move_vector = move_vector + d*unj;
        %%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
        if DEBUG
            Hc = rot(Hc, -theta)+tvector;
            plot(Hc(1), Hc(2), 'g.');
            plot([Hc(1) Hc(1)+d*unj(1)], ...
                 [Hc(2) Hc(2)+d*unj(2)], 'g');
        end
        %%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
    end
    
end
