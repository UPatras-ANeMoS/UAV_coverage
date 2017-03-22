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

% Add function path
addpath( genpath('Functions') );

% Load results
S1 = load('results_uniform_anisotropic_approx_20170310_1731');
S2 = load('results_uniform_anisotropic_norot_20170310_1759');
S3 = load('results_uniform_anisotropic_20170310_1837');
S = [S1 S2 S3];
colors = {'b', 'r', 'g'};


%%%%%%%%%%%%%%% Plot covered area %%%%%%%%%%%%%%%
figure
hold on
Tstep = 0;
smax = 0;
for i=1:length(S)
	plot( S(i).Tstep*linspace(1,S(i).smax,S(i).smax), ...
		100*S(i).cov_area/S(i).region_area, colors{i});
	if S(i).smax > smax
		smax = S(i).smax;
		Tstep = S(i).Tstep;
	end
end
axis([0 Tstep*smax 0 100]);
h = xlabel('$Time ~(s)$');
set(h,'Interpreter','latex')
h = ylabel('$A_{cov}~(\%)$');
set(h,'Interpreter','latex')

%%%%%%%%%%%%%%% Plot objective %%%%%%%%%%%%%%%
figure
hold on
Tstep = 0;
smax = 0;
H = 0;
for i=1:length(S)
	plot( S(i).Tstep*linspace(1,S(i).smax,S(i).smax), S(i).H, colors{i});
	if S(i).smax > smax
		smax = S(i).smax;
		Tstep = S(i).Tstep;
	end
	if max(S(i).H) > H
		H = max(S(i).H);
	end
end
axis([0 Tstep*smax 0 ceil(H)]);
h = xlabel('$Time ~(s)$');
set(h,'Interpreter','latex')
h = ylabel('$\mathcal{H}$');
set(h,'Interpreter','latex')
