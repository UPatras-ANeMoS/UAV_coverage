% Copyright 2016-2017 Sotiris Papatheodorou
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%    http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

clear variables
close all

addpath( genpath('../Functions') );

PLOT_NORMALS = 0;
PLOT_MOVE_VECTORS = 0;
COMPARE_GV = 0;

% Region
region = 8*[1 1 -1 -1 ; 1 -1 -1 1];

% Uncertainty centers
q = [-2 2 0;
    0 0 2];
q = [-2 1.65;
    -6 -6];
% Uncertainty radii
r = [0.2 0.4 0.5];
% r = [0.4 0.4 0.5];
% Guaranteed sensing radii
R_original = [3 1 2.4];
% R_original = [3.5 3.5 2.4];
R = R_original - r;
N = length(q(1,:));

% Select nodes
ii = 1;
jj = 2;


% AWGV cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AWGV = cell([1 N]);
for i=1:N
    AWGV{i} = AWGV_cell(region, q, r, R, i);
end


% Normal vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PLOT_NORMALS
    ni = zeros(size(AWGV{ii}));
    % Own cell normal vectors
    if ~isempty(AWGV{ii})
        for k=1:length(AWGV{ii}(1,:))
            a = (r(ii) + r(jj) + R(jj) - R(ii))/2;
            t = hyperbola_point2parameter(AWGV{ii}(:,k), q(:,ii), q(:,jj), a);
            ni(:,k) = Fni_AWGV(R(ii),R(jj),r(ii),r(jj),t,q(1,ii),q(1,jj),q(2,ii),q(2,jj));
        end


        % Point on hyperobola test
        p1 = AWGV{ii}(:,1);
        a = (r(ii) + r(jj) + R(jj) - R(ii))/2;
        c = norm(q(:,ii)-q(:,jj))/2;
        b = sqrt(c^2-a^2);
        theta = atan2(q(2,ii)-q(2,jj),q(1,ii)-q(1,jj));
        center = (q(:,ii)+q(:,jj))/2;
        t = hyperbola_point2parameter(p1, q(:,ii), q(:,jj), a);
        pp1 = rot([a*cosh(t) ; b*sinh(t)], theta) + center;
    end
end











% Free arc move vectors
mv_FA = C_integral_law_num( q(:,ii), AWGV{ii}, R(ii) );



% AWGV Move vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jni_AWGV = zeros(size(AWGV{ii}));
% Own cell move vectors
if ~isempty(AWGV{ii})
    for k=1:length(AWGV{ii}(1,:))
        a = (r(ii) + r(jj) + R(jj) - R(ii))/2;
        t = hyperbola_point2parameter(AWGV{ii}(:,k), q(:,ii), q(:,jj), a);
        Jni_AWGV(:,k) = FJni_AWGV(R(ii),R(jj),r(ii),r(jj),t,q(1,ii),q(1,jj),q(2,ii),q(2,jj));
    end
end

Jnj_AWGV = zeros(size(AWGV{jj}));
% Other cell move vectors
if ~isempty(AWGV{jj})
    for k=1:length(AWGV{jj}(1,:))
        a = (r(jj) + r(ii) + R(ii) - R(jj))/2;
        t = hyperbola_point2parameter(AWGV{jj}(:,k), q(:,jj), q(:,ii), a);
        Jnj_AWGV(:,k) = FJnj_AWGV(R(ii),R(jj),r(ii),r(jj),t,q(1,ii),q(1,jj),q(2,ii),q(2,jj));
    end
end


% Calculate total move vector
mv_AWGV = zeros(2,1);
mv_AWGV = mv_AWGV + mv_FA;
% Use Jn elements that are inside the sensing
if ~isempty(AWGV{ii})
    for k=1:length(AWGV{ii}(1,:))
        if norm(q(:,ii)-AWGV{ii}(:,k)) < R(ii) && ~sum(isnan(Jni_AWGV(:,k)))
            % Number of vertices
            Nv = length(AWGV{ii}(1,:));
            % Calculate length of current arc
            dl = ( norm(AWGV{ii}(:,k)-AWGV{ii}(:,mod(k-2,Nv)+1)) + ...
                    norm(AWGV{ii}(:,k)-AWGV{ii}(:,mod(k,Nv)+1)) ) / 2;
            % Add Jacobian * normal weighted by arc length
            mv_AWGV = mv_AWGV + Jni_AWGV(:,k) * dl;
        end
    end
end
if ~isempty(AWGV{jj})
    for k=1:length(AWGV{jj}(1,:))
        if norm(q(:,jj)-AWGV{jj}(:,k)) < R(jj) && ~sum(isnan(Jnj_AWGV(:,k)))
            % Number of vertices
            Nv = length(AWGV{jj}(1,:));
            % Calculate length of current arc
            dl = ( norm(AWGV{jj}(:,k)-AWGV{jj}(:,mod(k-2,Nv)+1)) + ...
                    norm(AWGV{jj}(:,k)-AWGV{jj}(:,mod(k,Nv)+1)) ) / 2;
            % Add Jacobian * normal weighted by arc length
            mv_AWGV = mv_AWGV + Jnj_AWGV(:,k) * dl;
        end
    end
end










% GV Move vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if COMPARE_GV
    Jni_GV = zeros(size(AWGV{ii}));
    % Own cell move vectors
    if ~isempty(AWGV{ii})
        for k=1:length(AWGV{ii}(1,:))
            a = (r(ii) + r(jj) + R(jj) - R(ii))/2;
            t = hyperbola_point2parameter(AWGV{ii}(:,k), q(:,ii), q(:,jj), a);
    %         Jni_GV(:,k) = FJni(a,t,q(1,ii),q(1,jj),q(2,ii),q(2,jj));
            Jni_GV(:,k) = FJni_GV(r(ii),r(jj),t,q(1,ii),q(1,jj),q(2,ii),q(2,jj));
        end
    end

    Jnj_GV = zeros(size(AWGV{jj}));
    % Other cell move vectors
    if ~isempty(AWGV{jj})
        for k=1:length(AWGV{jj}(1,:))
            a = (r(jj) + r(ii) + R(ii) - R(jj))/2;
            t = hyperbola_point2parameter(AWGV{jj}(:,k), q(:,jj), q(:,ii), a);
    %         Jnj_GV(:,k) = FJnj(a,t,q(1,ii),q(1,jj),q(2,ii),q(2,jj));
            Jnj_GV(:,k) = FJnj_GV(r(ii),r(jj),t,q(1,ii),q(1,jj),q(2,ii),q(2,jj));
        end
    end

    % Calculate total move vector
    mv_GV = zeros(2,1);
    mv_GV = mv_GV + mv_FA;
    % Use Jn elements that are inside the sensing
    if ~isempty(AWGV{ii})
        for k=1:length(AWGV{ii}(1,:))
            if norm(q(:,ii)-AWGV{ii}(:,k)) < R(ii) && ~sum(isnan(Jni_GV(:,k)))
                % Number of vertices
                Nv = length(AWGV{ii}(1,:));
                % Calculate length of current arc
                dl = ( norm(AWGV{ii}(:,k)-AWGV{ii}(:,mod(k-2,Nv)+1)) + ...
                        norm(AWGV{ii}(:,k)-AWGV{ii}(:,mod(k,Nv)+1)) ) / 2;
                % Add Jacobian * normal weighted by arc length
                mv_GV = mv_GV + Jni_GV(:,k) * dl;
            end
        end
    end
    if ~isempty(AWGV{jj})
        for k=1:length(AWGV{jj}(1,:))
            if norm(q(:,jj)-AWGV{jj}(:,k)) < R(jj) && ~sum(isnan(Jnj_GV(:,k)))
                % Number of vertices
                Nv = length(AWGV{jj}(1,:));
                % Calculate length of current arc
                dl = ( norm(AWGV{jj}(:,k)-AWGV{jj}(:,mod(k-2,Nv)+1)) + ...
                        norm(AWGV{jj}(:,k)-AWGV{jj}(:,mod(k,Nv)+1)) ) / 2;
                % Add Jacobian * normal weighted by arc length
                mv_GV = mv_GV + Jnj_GV(:,k) * dl;
            end
        end
    end
end








% PLOT PLOT PLOT PLOT PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
figure
hold on
plot(q(1,:), q(2,:), 'k.');
for i=1:N
    plot_circle(q(:,i), r(i), 'k');
    plot_circle(q(:,i), R(i), 'r--');
end
for i=1:N
    if ~isempty(AWGV{i})
        plot(AWGV{i}(1,:), AWGV{i}(2,:), 'b');
    end
end

% Normals
if PLOT_NORMALS
    for k=1:length(n(1,:))
        plot([AWGV{ii}(1,k) AWGV{ii}(1,k)+ni(1,k)], [AWGV{ii}(2,k) AWGV{ii}(2,k)+ni(2,k)], 'k');
    end
end

% Move vectors
if PLOT_MOVE_VECTORS
    for k=1:length(Jni_AWGV(1,:))
        plot([AWGV{ii}(1,k) AWGV{ii}(1,k)+Jni_AWGV(1,k)], [AWGV{ii}(2,k) AWGV{ii}(2,k)+Jni_AWGV(2,k)], 'r');
    end
    for k=1:length(Jnj_AWGV(1,:))
        plot([AWGV{jj}(1,k) AWGV{jj}(1,k)+Jnj_AWGV(1,k)], [AWGV{jj}(2,k) AWGV{jj}(2,k)+Jnj_AWGV(2,k)], 'r');
    end
    plot([q(1,ii) q(1,ii)+mv_AWGV(1)], [q(2,ii) q(2,ii)+mv_AWGV(2)], 'r');
end

% GV Move vectors
if PLOT_MOVE_VECTORS && COMPARE_GV
    for k=1:length(Jni_AWGV(1,:))
        plot([AWGV{ii}(1,k) AWGV{ii}(1,k)+Jni_GV(1,k)], [AWGV{ii}(2,k) AWGV{ii}(2,k)+Jni_GV(2,k)], 'k--');
    end
    for k=1:length(Jnj_AWGV(1,:))
        plot([AWGV{jj}(1,k) AWGV{jj}(1,k)+Jnj_GV(1,k)], [AWGV{jj}(2,k) AWGV{jj}(2,k)+Jnj_GV(2,k)], 'k--');
    end
    plot([q(1,ii) q(1,ii)+mv_GV(1)], [q(2,ii) q(2,ii)+mv_GV(2)], 'k--');
end

axis equal
