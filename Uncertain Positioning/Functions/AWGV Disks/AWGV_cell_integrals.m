% Copyright 2017 Sotiris Papatheodorou
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

function move_vector = AWGV_cell_integrals( Xi, uradiusi, celli, ...
    sradiusi, Xj, uradiusj, cellj, sradiusj, FULL, DEBUG )
% Calculate the integral over Hij, Hji and return the move vector
% Handle optional debug argument
if nargin < 10
    DEBUG = 0;
end
if nargin < 9
    FULL = 1;
end



% Hyperbola parameters
aij = (uradiusi + uradiusj + sradiusj - sradiusi) / 2;
aji = (uradiusi + uradiusj + sradiusi - sradiusj) / 2;


% c = norm( Xi-Xj ) / 2;
% bij = sqrt(c^2 - aij^2);
% bji = sqrt(c^2 - aji^2);
% if ~isreal(bij)
%     disp('bij is not real');
%     disp(aij);
%     disp(bij);
%     disp(c);
% end
% if ~isreal(bji)
%     disp('bji is not real');
%     disp(aji);
%     disp(bji);
%     disp(c);
% end

move_vector = zeros(2,1);




% Integral over Hij
% Loop over all vertices of GVi (the rotated and translated cell, as if
% produced by an East-West hyperbola). For each vertex check if it
% satisfies the East-West hyperbola equation. If it does find t and
% integrate
if ~isempty(celli)
    % Loop over all vertices
    for k=1:length( celli(1,:) )
        % Check if it's inside the sensing circle
        if norm( celli(:,k)-Xi ) <= sradiusi
            % Find the t parameter
            t = hyperbola_point2parameter(celli(:,k), Xi, Xj, aij);
            
            % Check whether it's a point of the hyperbola
            % t is NaN if it isn't, NaN ~= NaN check
            if t == t
                % Substitute the symbolic expressions
                xi = Xi(1);
                yi = Xi(2);
                xj = Xj(1);
                yj = Xj(2);

                uni = FJni_AWGV(sradiusi,sradiusj,uradiusi,uradiusj,...
                    t,xi,xj,yi,yj);

                % Find the length of the arc (half the distance from each adjacent vertex)
                % Cases for first and last vertices
                if k == 1
                    d = norm(celli(:,k)-celli(:,end)) / 2 + ...
                        norm(celli(:,k)-celli(:,k+1)) / 2;
                elseif k == length( celli(1,:) )
                    d = norm(celli(:,k)-celli(:,k-1)) / 2 + ...
                        norm(celli(:,k)-celli(:,1)) / 2;
                else
                    d = norm(celli(:,k)-celli(:,k-1)) / 2 + ...
                        norm(celli(:,k)-celli(:,k+1)) / 2;
                end
                % Add to the move vector
                move_vector = move_vector + d*uni;
                %%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
                if DEBUG
                    plot(celli(1,k), celli(2,k), 'm.');
                    plot([celli(1,k) celli(1,k)+d*uni(1)], ...
                         [celli(2,k) celli(2,k)+d*uni(2)], 'm');
                end
                %%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
            end
        end
    end
end




if FULL
    % Integral over Hji
    % Loop over all vertices of GVj (the rotated and translated cell, as if
    % produced by an East-West hyperbola). For each vertex check if it
    % satisfies the East-West hyperbola equation. If it does find t and
    % integrate
    if ~isempty(cellj)
        % Loop over all vertices
        for k=1:length( cellj(1,:) )
            % Check if it's inside the sensing circle
            if norm( cellj(:,k)-Xj ) <= sradiusj
                % Find the t parameter
                t = hyperbola_point2parameter(cellj(:,k), Xj, Xi, aji);

                % Check whether it's a point of the hyperbola
                % t is NaN if it isn't, NaN ~= NaN check
                if t == t
                    % Substitute the symbolic expressions
                    xi = Xi(1);
                    yi = Xi(2);
                    xj = Xj(1);
                    yj = Xj(2);

                    unj = FJnj_AWGV(sradiusi,sradiusj,uradiusi,uradiusj,...
                        t,xi,xj,yi,yj);

                    % Find the length of the arc (half the distance from each adjacent vertex)
                    % Cases for first and last vertices
                    if k == 1
                        d = norm(cellj(:,k)-cellj(:,end)) / 2 + ...
                            norm(cellj(:,k)-cellj(:,k+1)) / 2;
                    elseif k == length( cellj(1,:) )
                        d = norm(cellj(:,k)-cellj(:,k-1)) / 2 + ...
                            norm(cellj(:,k)-cellj(:,1)) / 2;
                    else
                        d = norm(cellj(:,k)-cellj(:,k-1)) / 2 + ...
                            norm(cellj(:,k)-cellj(:,k+1)) / 2;
                    end
                    % Add to the move vector
                    move_vector = move_vector + d*unj;
                    %%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
                    if DEBUG
                        plot(cellj(1,k), cellj(2,k), 'm.');
                        plot([cellj(1,k) cellj(1,k)+d*unj(1)], ...
                             [celli(2,k) celli(2,k)+d*unj(2)], 'm');
                    end
                    %%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%
                end
            end
        end
    end
end
