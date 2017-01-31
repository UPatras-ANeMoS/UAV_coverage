% MIT License
% 
% Copyright (c) 2017 Sotiris Papatheodorou
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

function plot_sim_UAV(sim)

if sim.PLOT_STATE_3D || sim.PLOT_STATE_2D || sim.PLOT_STATE_QUALITY || sim.PLOT_STATE_PHI
    
    
    
    % ----------------- Plot network 2D -----------------
    if sim.PLOT_STATE_2D
        clf
        hold on
        % Region
        plot_poly( sim.region, 'k');
        % Sensing disks and cells
        for i=1:sim.N
            plot_poly( sim.C{i}, 'r--');
            plot_poly( sim.W{i}, 'k');
        end
        % Node positions
        for i=1:sim.N
%                 tmpc = [sim.X(i) + disk_rad * cos(t) ; sim.Y(i) + disk_rad * sin(t)];
%                 fill( tmpc(1,:), tmpc(2,:), 'k', 'EdgeColor', 'none' );
            plot( sim.X(i), sim.Y(i), 'k.' )
            hold on
        end
        % Connectivity
        if sim.PLOT_COMMS
            for i=1:sim.N
                for j=1:sim.N
                    if sim.A(i,j)
                        plot([sim.X(i) sim.X(j)], [sim.Y(i) sim.Y(j)], 'g')
                    end
                end
            end
        end
        plot_AABB(sim.axis, 'w.');

        set( gca, 'Units', 'normalized', 'Position', [0 0 1 1] );
        axis(sim.axis)
        axis equal
        axis off

        if sim.SAVE_PLOTS
            fname = strcat( '~/Frames/', sprintf('2D_frame_%d.png', sim.s) );
            print(fname,'-dpng','-r600');
%             saveas(gcf, fname);
        else
            pause(0.01);
        end
    end

    
    
    % ----------------- Plot network 3D -----------------
    if sim.PLOT_STATE_3D
        clf
        hold on
        % Sensing disks and cells
        for i=1:sim.N
            plot3_poly( [sim.C{i} ; zeros(size(sim.C{i}(1,:)))], 'r--');
            plot3_poly( [sim.W{i} ; zeros(size(sim.W{i}(1,:)))], 'k');
        end
        % Node positions and cones
        for i=1:sim.N
            plot3( sim.X(i), sim.Y(i), sim.Z(i), 'ko' )
            plot3( [sim.X(i) sim.X(i)], [sim.Y(i) sim.Y(i)], [sim.Z(i) 0], 'k--.' )
            for j=1:24:sim.PPC
                plot3([sim.C{i}(1,j) sim.X(i)], [sim.C{i}(2,j) sim.Y(i)], [0 sim.Z(i)], 'r--');
            end
        end
        % Connectivity
        if sim.PLOT_COMMS
            for i=1:sim.N
                for j=1:sim.N
                    if sim.A(i,j)
                        plot3([sim.X(i) sim.X(j)], [sim.Y(i) sim.Y(j)], [sim.Z(i) sim.Z(j)], 'g')
                    end
                end
            end
        end
        % Plot region
        plot3_poly( [sim.region ; zeros(size(sim.region(1,:)))], 'k' );
        plot3_AABB([sim.axis 0 sim.zmax], 'w.');

        set( gca, 'Units', 'normalized', 'Position', [0 0 1 1] );
        view(-16, 34);
        axis([sim.axis 0 sim.zmax])
        axis equal
        axis off

        if sim.SAVE_PLOTS
            fname = strcat( '~/Frames/', sprintf('3D_frame_%d.png', sim.s) );
            print(fname,'-dpng','-r600');
%             saveas(gcf, fname);
        else
            pause(0.01);
        end
    end

    
    
    % ----------------- Plot network phi -----------------
    if sim.PLOT_STATE_PHI
        clf
        hold on
        plot_phi( sim.phi , sim.region );
        % Region
        plot_poly( sim.region, 'k');
        % Sensing disks and cells
        for i=1:sim.N
            plot_poly( sim.C{i}, 'r--');
            plot_poly( sim.W{i}, 'k');
        end
        % Node positions
        for i=1:sim.N
%                 tmpc = [sim.X(i) + disk_rad * cos(t) ; sim.Y(i) + disk_rad * sin(t)];
%                 fill( tmpc(1,:), tmpc(2,:), 'k', 'EdgeColor', 'none' );
            plot( sim.X(i), sim.Y(i), 'k.' )
            hold on
        end
        % Connectivity
        if sim.PLOT_COMMS
            for i=1:sim.N
                for j=1:sim.N
                    if sim.A(i,j)
                        plot([sim.X(i) sim.X(j)], [sim.Y(i) sim.Y(j)], 'g')
                    end
                end
            end
        end
        plot_AABB(sim.axis, 'w.');

        set( gca, 'Units', 'normalized', 'Position', [0 0 1 1] );
        axis(sim.axis)
        axis equal
        axis off

        if sim.SAVE_PLOTS
            fname = strcat( '~/Frames/', sprintf('PHI_frame_%d.png', sim.s) );
            print(fname,'-dpng','-r600');
%             saveas(gcf, fname);
        else
            pause(0.01);
        end
    end
    
    
    
    % ----------------- Plot network quality -----------------
    if sim.PLOT_STATE_QUALITY
        clf
        hold on
        % Plot cylinders
        for i=1:sim.N
            plot3_cell_quality(sim.W{i}, sim.f(i), 'r');
        end
        % Plot region
        plot3_poly( [sim.region ; zeros(size(sim.region(1,:)))], 'k' );
        plot3_AABB([sim.axis 0 sim.zmax], 'w.');

        set( gca, 'Units', 'normalized', 'Position', [0 0 1 1] );
        view(-16, 34);
        axis([sim.axis 0 1])
        axis equal
        axis off

        if sim.SAVE_PLOTS
            fname = strcat( '~/Frames/', sprintf('Q_frame_%d.png', sim.s) );
            print(fname,'-dpng','-r600');
%             saveas(gcf, fname);
        else
            pause(0.01);
        end
    end
end
    