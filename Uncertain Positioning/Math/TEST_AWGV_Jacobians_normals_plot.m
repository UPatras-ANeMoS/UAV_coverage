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

clear variables
close all

PLOT_JACOBIAN = 1;
PLOT_NORMAL = 1;
PLOT_PRODUCT = 1;
cpp_result_dir = '../../NRobot/bin/';

% Setup 2 agents
c = 1;
xi = -c;
yi = 0;
ri = 0.05;
Ri = 0.9;
xj = c;
yj = 0;
rj = 0.1;
Rj = 0.8;
ai = ri + rj + Rj - Ri;
bi = sqrt(c^2-ai^2);
aj = ri + rj + Ri - Rj;
bj = sqrt(c^2-aj^2);

tmax = ceil( acosh(15/abs(ai)) );
t = linspace(-tmax, tmax, 100);

% Compute hyperbola branches
Hij = [ai*cosh(t) ; bi*sinh(t)];
Hji = [aj*cosh(t) ; bj*sinh(t)];

% Compute Jacobian matrices
Ji = zeros(2,2,length(t));
Jj = zeros(2,2,length(t));
for k=1:length(t)
	Ji(:,:,k) = FJi_AWGV(Ri,Rj,ri,rj,t(k),xi,xj,yi,yj);
	Jj(:,:,k) = FJj_AWGV(Ri,Rj,ri,rj,t(k),xi,xj,yi,yj);
end

% Compute normal vectors
ni = Fni_AWGV(Ri,Rj,ri,rj,t,xi,xj,yi,yj);
nj = Fnj_AWGV(Ri,Rj,ri,rj,t,xi,xj,yi,yj);

% Compute Jacobian-normal products
Jni = FJni_AWGV(Ri,Rj,ri,rj,t,xi,xj,yi,yj);
Jnj = FJnj_AWGV(Ri,Rj,ri,rj,t,xi,xj,yi,yj);

% Load C++ results
d = importdata([cpp_result_dir 'Jn_values.txt']);
t_cpp = d(:,1)';
Jni_x_cpp = d(:,2)';
Jni_y_cpp = d(:,3)';
Jnj_x_cpp = d(:,4)';
Jnj_y_cpp = d(:,5)';
d = importdata([cpp_result_dir 'n_values.txt']);
t_cpp = d(:,1)';
ni_x_cpp = d(:,2)';
ni_y_cpp = d(:,3)';
nj_x_cpp = d(:,4)';
nj_y_cpp = d(:,5)';
d = importdata([cpp_result_dir 'J_values.txt']);
t_cpp = d(:,1)';
Ji_x_xi_cpp = d(:,2)';
Ji_x_yi_cpp = d(:,3)';
Ji_y_xi_cpp = d(:,4)';
Ji_y_yi_cpp = d(:,5)';
Jj_x_xi_cpp = d(:,6)';
Jj_x_yi_cpp = d(:,7)';
Jj_y_xi_cpp = d(:,8)';
Jj_y_yi_cpp = d(:,9)';

% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Jacobian
if PLOT_JACOBIAN
	figure('Name','Jacobian')
	hold on
	
	% Ji
	subplot(2,4,1)
	hold on
	plot(t, squeeze(Ji(1,1,:)), 'b')
	plot(t_cpp, Ji_x_xi_cpp, 'r--')
	plot(t, zeros(size(t)), 'k--')
	ylabel('J_i x x_i')
	xlabel('parameter t')

	subplot(2,4,2)
	hold on
	plot(t, squeeze(Ji(1,2,:)), 'b')
	plot(t_cpp, Ji_x_yi_cpp, 'r--')
	plot(t, zeros(size(t)), 'k--')
	ylabel('J_i x y_i')
	xlabel('parameter t')

	subplot(2,4,5)
	hold on
	plot(t, squeeze(Ji(2,1,:)), 'b')
	plot(t_cpp, Ji_y_xi_cpp, 'r--')
	plot(t, zeros(size(t)), 'k--')
	ylabel('J_i y x_i')
	xlabel('parameter t')

	subplot(2,4,6)
	hold on
	plot(t, squeeze(Ji(2,2,:)), 'b')
	plot(t_cpp, Ji_y_yi_cpp, 'r--')
	plot(t, zeros(size(t)), 'k--')
	ylabel('J_i y y_i')
	xlabel('parameter t')
	
	% Jj
	subplot(2,4,3)
	hold on
	plot(t, squeeze(Jj(1,1,:)), 'b')
	plot(t_cpp, Jj_x_xi_cpp, 'r--')
	plot(t, zeros(size(t)), 'k--')
	ylabel('J_j x x_i')
	xlabel('parameter t')

	subplot(2,4,4)
	hold on
	plot(t, squeeze(Jj(1,2,:)), 'b')
	plot(t_cpp, Jj_x_yi_cpp, 'r--')
	plot(t, zeros(size(t)), 'k--')
	ylabel('J_j x y_i')
	xlabel('parameter t')

	subplot(2,4,7)
	hold on
	plot(t, squeeze(Jj(2,1,:)), 'b')
	plot(t_cpp, Jj_y_xi_cpp, 'r--')
	plot(t, zeros(size(t)), 'k--')
	ylabel('J_j y x_i')
	xlabel('parameter t')

	subplot(2,4,8)
	hold on
	plot(t, squeeze(Jj(2,2,:)), 'b')
	plot(t_cpp, Jj_y_yi_cpp, 'r--')
	plot(t, zeros(size(t)), 'k--')
	ylabel('J_j y y_i')
	xlabel('parameter t')
end

% Plot normals
if PLOT_NORMAL
	figure('Name','Normal')
	hold on
	
	plot(xi, yi, 'ro');
	plot(xj, yj, 'ko');
	plot(Hij(1,:), Hij(2,:), 'r');
	plot(Hji(1,:), Hji(2,:), 'k');
	
	for k=1:length(t_cpp)
		% Matalb
		plot([Hij(1,k) Hij(1,k)+ni(1,k)], ...
			 [Hij(2,k) Hij(2,k)+ni(2,k)], 'r');
		plot([Hji(1,k) Hji(1,k)+nj(1,k)], ...
			 [Hji(2,k) Hji(2,k)+nj(2,k)], 'k');
		% C++
		plot([Hij(1,k) Hij(1,k)+ni_x_cpp(k)], ...
			 [Hij(2,k) Hij(2,k)+ni_y_cpp(k)], 'm--');
		plot([Hji(1,k) Hji(1,k)+nj_x_cpp(k)], ...
			 [Hji(2,k) Hji(2,k)+nj_y_cpp(k)], 'b--');
	end
	
	axis square
	axis([-5 5 -5 5])
end

% Plot Jacobian-normal product
if PLOT_PRODUCT
	figure('Name','Product')

	subplot(2,2,1)
	hold on
	plot(t, Jni(1,:), 'b')
	plot(t_cpp, Jni_x_cpp, 'r--')
	plot(t, zeros(size(t)), 'k--')
	ylabel('Jn_i x')
	xlabel('parameter t')

	subplot(2,2,3)
	hold on
	plot(t, Jni(2,:), 'b')
	plot(t_cpp, Jni_y_cpp, 'r--')
	plot(t, zeros(size(t)), 'k--')
	ylabel('Jn_i y')
	xlabel('parameter t')

	subplot(2,2,2)
	hold on
	plot(t, Jnj(1,:), 'b')
	plot(t_cpp, Jnj_x_cpp, 'r--')
	plot(t, zeros(size(t)), 'k--')
	ylabel('Jn_j x')
	xlabel('parameter t')

	subplot(2,2,4)
	hold on
	plot(t, Jnj(2,:), 'b')
	plot(t_cpp, Jnj_y_cpp, 'r--')
	plot(t, zeros(size(t)), 'k--')
	ylabel('Jn_j y')
	xlabel('parameter t')

end
