clc; clear; close all;
tic

%% 1. Adjustable parameters
% This section sets the computational domain, spatial grid sizes, and time parameters.
% x0,x1 define the finite computational interval [x0,x1] in the x direction;
% y0,y1 define the finite computational interval [y0,y1] in the y direction.
% Nx,Ny are the Fourier grid sizes in x and y, respectively, and need not be equal.
x0 = -15;  x1 = 15;  Nx = 64*2;
y0 = -15;  y1 = 15;  Ny = 64*2;

% t0 is the initial time, and t_end is the final time.
% The angular velocity is 1/4, so one circular orbit takes 8*pi units of time.
% To end the calculation exactly at the end of one period, first specify
% the number of steps per period, then calculate the RK4 time step dt.
% save_step is the number of RK4 time steps between saved frames;
% here save_step = 100, so save_dt = 100*dt.
t0 = 0;
t_end = 8*pi;
period_steps = 10000;
dt = (t_end-t0)/period_steps;
save_step = 100;
save_dt = save_step*dt;

% The Fourier wavenumbers are ordered as [0,1,...,N/2-1,-N/2,...,-1],
% so Nx and Ny must be even.
if mod(Nx,2) ~= 0 || mod(Ny,2) ~= 0
    error('Nx and Ny must be even for the Fourier wave-number ordering used here.');
end

%% 2. Grid and Fourier wavenumbers
% dx,dy are the spatial step sizes in the two directions. A periodic Fourier
% pseudospectral method is used, so the grids end at x1-dx and y1-dy without duplicating endpoints.
Lx = x1 - x0;  dx = Lx/Nx;
Ly = y1 - y0;  dy = Ly/Ny;

% X,Y are one-dimensional column vectors; meshgrid generates the two-dimensional grids x,y.
% In the MATLAB arrays, the row index corresponds to y and the column index to x,
% so the two-dimensional arrays x,y,A all have size Ny-by-Nx.
X = (x0 + (0:Nx-1)*dx).';
Y = (y0 + (0:Ny-1)*dy).';
[x,y] = meshgrid(X,Y);

% Kx,Ky are the wavenumbers in Fourier space.
% For a periodic function A,
%   F[A_xx + A_yy] = -(kx^2+ky^2) F[A].
Kx = (2*pi/Lx)*[0:Nx/2-1 -Nx/2:-1].';
Ky = (2*pi/Ly)*[0:Ny/2-1 -Ny/2:-1].';
[kx,ky] = meshgrid(Kx,Ky);
k2 = kx.^2 + ky.^2;

fprintf('space: x in [%g,%g], Nx = %d; y in [%g,%g], Ny = %d\n', ...
    x0,x1,Nx,y0,y1,Ny);
fprintf('time : t0 = %g, t_end = %g, dt = %.6g, save_dt = %.6g\n', ...
    t0,t_end,dt,save_dt);

%% 3. Parameters of the circularly moving dromion
% The peak of the exact solution moves along a circle of radius R.
% Omega is the angular velocity, R*Omega is the tangential speed, and density_peak is the constant peak intensity.
R = 8;
Omega = 1/4;
c = (2 + sqrt(2))/2;
density_peak = 1;

% Fix the z-axis and color limits of the three-dimensional plots.
% The numerical and exact intensity plots use the same limits; the difference plot is symmetric about zero.
% If the difference plot is clipped, increase diff_plot_max.
density_plot_range = [0, 1.05*density_peak];
diff_plot_max = 0.005*density_peak;
diff_plot_range = [-diff_plot_max, diff_plot_max];

% By default, a full-period run saves only the error plot and the manuscript's
% one-row, two-column validation plot. Set this flag to true to save three-panel frames as well.
save_frame_images = strcmpi(getenv('CIRCULAR_SAVE_FRAMES'),'true');

%% 4. Output times and output directory
% Output times must lie on the RK4 time grid.
% Since save_dt = save_step*dt, generate the output times from integer step indices.
Nt = round((t_end-t0)/dt);
if abs(t0 + Nt*dt - t_end) > 100*eps(max(1,abs(t_end)))
    error('t_end - t0 must be an integer multiple of dt.');
end
if mod(Nt,save_step) ~= 0
    error('t_end - t0 must be an integer multiple of save_dt.');
end

save_steps = 0:save_step:Nt;
t_save = t0 + save_steps*dt;
nSave = numel(t_save);

% The four snapshots correspond to angles 0, pi/2, pi, and 3*pi/2.
snapshot_steps = [0,period_steps/4,period_steps/2,3*period_steps/4];
if any(mod(snapshot_steps,save_step) ~= 0)
    error('The four quarter-period snapshots must fall on saved time steps.');
end
snapshot_times = t0 + snapshot_steps*dt;
snapshot_density = zeros(Ny,Nx,numel(snapshot_steps));
snapshot_found = false(size(snapshot_steps));

outDir = fullfile(fileparts(mfilename('fullpath')),'circular_dromion_results');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

%% 5. Exact initial condition
% Evaluate the exact circular-motion solution on the grid at t=t0 to initialize the main field.
% The two exponential variables describe localization in x and y, and den is their common denominator.
xc = R*cos(Omega*t0);
yc = R*sin(Omega*t0);
ex = exp(2*c*(x-xc) - log(2)/2);
ey = exp(2*c*(y-yc) - log(2)/2);
den = 1 + ex + ey + 2*ex.*ey;
phase = -x*sin(Omega*t0) + y*cos(Omega*t0) ...
    + (2 + 2*sqrt(2))*t0;
A0 = 2*(1 + sqrt(2))*exp(c*(x+y-xc-yc) - log(2)/2) ...
    .*exp(1i*phase)./den;

% fut contains the two-dimensional Fourier coefficients of A. All subsequent time stepping is performed in Fourier space.
fut = fft2(A0);

% I1 is the first conserved quantity, and E_I1 is its relative error.
% E_inf,E_2 are the infinity and L2 norms of the intensity difference |A_num|^2-|A_ex|^2, respectively.
I1 = zeros(nSave,1);
E_I1 = zeros(nSave,1);
E_inf = zeros(nSave,1);
E_2 = zeros(nSave,1);

% Simpson's rule requires the midpoint value in each subinterval.
% These midpoints are not on the original Fourier grid; their values are obtained by cubic spline interpolation below.
y_mid = Y(1:end-1) + dy/2;
x_mid = X(1:end-1) + dx/2;

%% 6. Time stepping and output
% Create the three-panel frame window only when save_frame_images=true.
if save_frame_images
fig = figure('Color','w','Visible','on','Position',[100 100 1500 450]);
tiledlayout(fig,1,3,'Padding','compact','TileSpacing','compact');

ax_num = nexttile;
p_num = surf(x,y,zeros(Ny,Nx),'EdgeColor','none');
title(ax_num,'|A_{num}|^2'); xlabel(ax_num,'x'); ylabel(ax_num,'y'); zlabel(ax_num,'|A|^2');
view(ax_num,42,28); axis(ax_num,'tight'); colorbar(ax_num);
xlim(ax_num,[x0 x1]); ylim(ax_num,[y0 y1]);
zlim(ax_num,density_plot_range);
clim(ax_num,density_plot_range);

ax_exact = nexttile;
p_exact = surf(x,y,zeros(Ny,Nx),'EdgeColor','none');
title(ax_exact,'|A_{ex}|^2'); xlabel(ax_exact,'x'); ylabel(ax_exact,'y'); zlabel(ax_exact,'|A|^2');
view(ax_exact,42,28); axis(ax_exact,'tight'); colorbar(ax_exact);
xlim(ax_exact,[x0 x1]); ylim(ax_exact,[y0 y1]);
zlim(ax_exact,density_plot_range);
clim(ax_exact,density_plot_range);

ax_diff = nexttile;
p_diff = surf(x,y,zeros(Ny,Nx),'EdgeColor','none');
title(ax_diff,'|A_{num}|^2-|A_{ex}|^2'); xlabel(ax_diff,'x'); ylabel(ax_diff,'y'); zlabel(ax_diff,'difference');
view(ax_diff,42,28); axis(ax_diff,'tight'); colorbar(ax_diff);
xlim(ax_diff,[x0 x1]); ylim(ax_diff,[y0 y1]);
zlim(ax_diff,diff_plot_range);
clim(ax_diff,diff_plot_range);
end

% The outer loop visits all output times.
% These times are integer multiples of dt, so the inner loop advances a fixed number of RK4 steps.
current_step = 0;
for si = 1:nSave
    target_step = save_steps(si);
    target_t = t_save(si);

    for step = current_step+1:target_step
        % Every step uses the fixed time step dt.
        % step is the index of the step to be completed, starting at t_now.
        t_now = t0 + (step-1)*dt;
        fut0 = fut;

        % Classical fourth-order Runge--Kutta uses four stages K1,K2,K3,K4.
        % At each stage, recover A, update the boundaries at the stage time,
        % reconstruct U,V, and evaluate the right-hand side in Fourier space.
        for rk = 1:4
            if rk == 1
                t_stage = t_now;
                fut_stage = fut0;
            elseif rk == 2
                t_stage = t_now + dt/2;
                fut_stage = fut0 + dt*K1/2;
            elseif rk == 3
                t_stage = t_now + dt/2;
                fut_stage = fut0 + dt*K2/2;
            else
                t_stage = t_now + dt;
                fut_stage = fut0 + dt*K3;
            end

            A_stage = ifft2(fut_stage);
            absA2 = abs(A_stage).^2;

            % Fourier pseudospectral differentiation:
            %   partial_x |A|^2 = F^{-1}[ i*kx*F(|A|^2) ],
            %   partial_y |A|^2 = F^{-1}[ i*ky*F(|A|^2) ].
            % real(...) removes tiny imaginary parts caused by roundoff.
            absA2_hat = fft2(absA2);
            absA2_x = real(ifft2(1i*kx.*absA2_hat));
            absA2_y = real(ifft2(1i*ky.*absA2_hat));

            % Mean-flow boundary values used in the numerical calculation.
            % U_bottom is a 1-by-Nx row on the lower boundary y=y0;
            % V_left is a Ny-by-1 column on the left boundary x=x0.
            U_bottom = (3 + 2*sqrt(2))*sech(c*(X.' ...
                - R*cos(Omega*t_stage)) - log(2)/4).^2 ...
                - (X.'/4)*cos(Omega*t_stage);
            V_left = (3 + 2*sqrt(2))*sech(c*(Y ...
                - R*sin(Omega*t_stage)) - log(2)/4).^2 ...
                - (Y/4)*sin(Omega*t_stage);

            % Reconstruct U:
            % The potential equation is U_y = partial_x |A|^2.
            % At each fixed x_m, integrate in y from the lower boundary.
            % Apply Simpson's rule on each subinterval, using cubic spline values at the midpoints.
            absA2_x_mid = interp1(Y,absA2_x,y_mid,'spline');
            U_step = dy*(absA2_x(1:end-1,:) + 4*absA2_x_mid ...
                + absA2_x(2:end,:))/6;
            U = [U_bottom; U_bottom + cumsum(U_step,1)];

            % Reconstruct V:
            % The potential equation is V_x = partial_y |A|^2.
            % At each fixed y_n, integrate in x from the left boundary.
            % Again, use Simpson's rule with cubic spline values at the midpoints.
            absA2_y_mid = interp1(X,absA2_y.',x_mid,'spline').';
            V_step = dx*(absA2_y(:,1:end-1) + 4*absA2_y_mid ...
                + absA2_y(:,2:end))/6;
            V = [V_left, V_left + cumsum(V_step,2)];

            % The DS-I main-field equation is written as
            %   A_t = i*(A_xx + A_yy + (U+V)A).
            % In Fourier space, A_xx+A_yy corresponds to -k2.*fut_stage;
            % evaluate the nonlinear potential term in physical space, then apply fft2.
            K_stage = 1i*(-k2.*fut_stage + fft2((U+V).*A_stage));

            if rk == 1
                K1 = K_stage;
            elseif rk == 2
                K2 = K_stage;
            elseif rk == 3
                K3 = K_stage;
            else
                K4 = K_stage;
            end
        end

        % RK4 weighted average: the weights of K1,K2,K3,K4 are 1:2:2:1.
        fut = fut0 + dt*(K1 + 2*K2 + 2*K3 + K4)/6;
    end
    current_step = target_step;

    %% 7. Exact solution, errors, and conserved quantity at the current output time
    A_num = ifft2(fut);

    % Reevaluate the exact circular-motion solution at the same output time target_t
    % for intensity plots and error comparisons with the numerical solution.
    xc = R*cos(Omega*target_t);
    yc = R*sin(Omega*target_t);
    ex = exp(2*c*(x-xc) - log(2)/2);
    ey = exp(2*c*(y-yc) - log(2)/2);
    den = 1 + ex + ey + 2*ex.*ey;
    phase = -x*sin(Omega*target_t) + y*cos(Omega*target_t) ...
        + (2 + 2*sqrt(2))*target_t;
    A_exact = 2*(1 + sqrt(2))*exp(c*(x+y-xc-yc) - log(2)/2) ...
        .*exp(1i*phase)./den;

    density_num = abs(A_num).^2;
    density_exact = abs(A_exact).^2;
    density_diff = density_num - density_exact;

    % Discrete approximation of the conserved quantity I1 = integral |A|^2 dxdy.
    % E_I1 measures its relative drift with respect to the initial value of I1.
    I1(si) = sum(density_num,'all')*dx*dy;
    if si == 1
        I1_0 = I1(si);
    end
    E_I1(si) = abs(I1(si)-I1_0)/abs(I1_0);
    E_inf(si) = max(abs(density_diff(:)));
    E_2(si) = sqrt(sum(density_diff(:).^2)*dx*dy);

    snapshot_index = find(target_step == snapshot_steps,1);
    if ~isempty(snapshot_index)
        snapshot_density(:,:,snapshot_index) = density_num;
        snapshot_found(snapshot_index) = true;
    end

    %% 8. Display and save three-dimensional plots at the current output time
    % At each output time, update the three-panel display and save it as a PNG.
    % Frames are numbered frame_001.png, frame_002.png, ...;
    % rerunning the script overwrites images with the same names.
    if save_frame_images
    set(p_num,'ZData',density_num,'CData',density_num);
    set(p_exact,'ZData',density_exact,'CData',density_exact);
    set(p_diff,'ZData',density_diff,'CData',density_diff);

    title(ax_num,sprintf('|A_{num}|^2, t = %.3f',target_t));
    title(ax_exact,sprintf('|A_{ex}|^2, t = %.3f',target_t));
    title(ax_diff,sprintf('|A_{num}|^2-|A_{ex}|^2, t = %.3f',target_t));

    % Use the same z-axis and color limits at every output time.
    % This allows direct comparison across times without automatic rescaling changing the color scale.
    zlim(ax_num,density_plot_range);
    clim(ax_num,density_plot_range);
    zlim(ax_exact,density_plot_range);
    clim(ax_exact,density_plot_range);
    zlim(ax_diff,diff_plot_range);
    clim(ax_diff,diff_plot_range);

    sgtitle(fig,sprintf('circular dromion, t = %.3f',target_t));
    drawnow;

    saveas(fig,fullfile(outDir,['frame_',num2str(si,'%03d'),'.png']));
    end

    fprintf('saved %3d/%3d, t = %+7.3f, E_I1 = %.3e, E_inf = %.3e\n', ...
        si,nSave,target_t,E_I1(si),E_inf(si));
end

%% 9. Error history
% Plot three time-dependent diagnostics; this section saves an image, not a separate MAT file.
%   E_I1 : relative error in the conserved quantity I1;
%   E_inf: maximum absolute intensity difference;
%   E_2  : discrete L2 norm of the intensity difference.
fig_error = figure('Color','w','Visible','on');
semilogy(t_save,max(E_I1,eps),'LineWidth',1.5); hold on;
semilogy(t_save,max(E_inf,eps),'LineWidth',1.5);
semilogy(t_save,max(E_2,eps),'LineWidth',1.5);
grid on;
xlabel('t'); ylabel('error');
legend('E_M','E_\infty','E_2','Location','best');
title('circular dromion error history');
drawnow;
saveas(fig_error,fullfile(outDir,'error_history.png'));

%% 10. One-row, two-column numerical validation figure for the manuscript
% The left panel overlays numerical intensity contours at four quarter-period times;
% the right panel shows the relative error E_I1 and the intensity infinity-norm error E_inf.
if ~all(snapshot_found)
    error('Not all four quarter-period snapshots were stored.');
end

snapshot_colors = [
    0.00 0.00 0.00
    0.85 0.10 0.10
    0.00 0.35 0.85
    0.00 0.55 0.20];
snapshot_labels = {'$t=0$','$t=2\pi$','$t=4\pi$','$t=6\pi$'};
contour_levels = [0.20 0.50 0.80];
label_positions = [
     10.8   1.8
      1.8  10.5
    -10.8  -1.8
      1.8 -10.5];

% Match the standalone replotting script: a text width of 16.6 cm and two plotting areas of equal width and height.
figure_width = 16.6;
panel_width = 6.45;
left_margin = 1.20;
column_gap = 2.00;
bottom_margin = 1.10;
figure_height = bottom_margin + panel_width + 0.25;
fig_validation = figure('Color','w','Units','centimeters', ...
    'Position',[2,2,figure_width,figure_height],'Visible','off');
positions = [left_margin bottom_margin panel_width panel_width; ...
    left_margin+panel_width+column_gap bottom_margin panel_width panel_width];

ax_contour = axes(fig_validation,'Units','centimeters','Position',positions(1,:));
hold(ax_contour,'on');
orbit_angle = linspace(0,2*pi,600);
plot(ax_contour,R*cos(orbit_angle),R*sin(orbit_angle),'--', ...
    'Color',[0.62 0.62 0.62],'LineWidth',0.70,'HandleVisibility','off');

% Cubic splines are used only to display smooth contours; the discrete time-stepping data are unchanged.
x_plot = linspace(x0,x1,401);
y_plot = linspace(y0,y1,401);
[x_plot_grid,y_plot_grid] = meshgrid(x_plot,y_plot);
legend_handles = gobjects(1,numel(snapshot_times));
for js = 1:numel(snapshot_times)
    density_plot = interp2(x,y,snapshot_density(:,:,js), ...
        x_plot_grid,y_plot_grid,'spline');
    density_plot = max(real(density_plot),0);
    contour(ax_contour,x_plot_grid,y_plot_grid,density_plot,contour_levels, ...
        'Color',snapshot_colors(js,:),'LineWidth',1.10, ...
        'HandleVisibility','off');
    legend_handles(js) = plot(ax_contour,nan,nan,'-', ...
        'Color',snapshot_colors(js,:),'LineWidth',1.10);

    text(ax_contour,label_positions(js,1),label_positions(js,2), ...
        snapshot_labels{js},'Interpreter','latex', ...
        'Color',snapshot_colors(js,:),'FontName','Times New Roman', ...
        'FontSize',10.0,'HorizontalAlignment','center', ...
        'VerticalAlignment','middle');
end
axis(ax_contour,'equal');
xlim(ax_contour,[x0 x1]); ylim(ax_contour,[y0 y1]);
xticks(ax_contour,[-15,-8,0,8,15]);
yticks(ax_contour,[-15,-8,0,8,15]);
xlabel(ax_contour,'$x$','Interpreter','latex','FontSize',13.0);
ylabel(ax_contour,'$y$','Interpreter','latex','FontSize',13.0);
set(ax_contour,'FontName','Times New Roman','FontSize',10.5, ...
    'TickLabelInterpreter','latex','LineWidth',0.95, ...
    'TickDir','out','Layer','top','Box','on');
grid(ax_contour,'on');
set(ax_contour,'GridColor',[0.85 0.85 0.85],'GridAlpha',0.55);
legend(ax_contour,legend_handles,snapshot_labels,'Interpreter','latex', ...
    'Location','northwest','FontSize',9.5,'Box','on');

ax_error = axes(fig_validation,'Units','centimeters','Position',positions(2,:));
semilogy(ax_error,t_save,max(E_I1,eps),'-','Color',[0.00 0.20 0.65], ...
    'LineWidth',1.15); hold(ax_error,'on');
semilogy(ax_error,t_save,max(E_inf,eps),'-','Color',[0.85 0.10 0.10], ...
    'LineWidth',1.15);
xlim(ax_error,[t0 t_end]);
xticks(ax_error,[0,2*pi,4*pi,6*pi,8*pi]);
xticklabels(ax_error,{'$0$','$2\pi$','$4\pi$','$6\pi$','$8\pi$'});
xlabel(ax_error,'$t$','Interpreter','latex','FontSize',13.0);
ylabel(ax_error,'Error','Interpreter','latex','FontSize',13.0);
set(ax_error,'FontName','Times New Roman','FontSize',10.5, ...
    'TickLabelInterpreter','latex','LineWidth',0.95, ...
    'TickDir','out','Layer','top','Box','on');
grid(ax_error,'on');
legend(ax_error,{'$E_M$','$E_{\infty}$'},'Interpreter','latex', ...
    'Location','best','FontSize',10.0,'Box','on');

validation_file = fullfile(outDir, ...
    'circular_dromion_numerical_validation_600dpi.png');
pbaspect(ax_contour,[1 1 1]);
pbaspect(ax_error,[1 1 1]);
set(ax_contour,'Position',positions(1,:));
set(ax_error,'Position',positions(2,:));
drawnow;
set(fig_validation,'PaperUnits','centimeters', ...
    'PaperPosition',[0,0,figure_width,figure_height], ...
    'PaperSize',[figure_width,figure_height]);
print(fig_validation,validation_file,'-dpng','-r600');
fprintf('Saved: %s\n',validation_file);
close(fig_validation);

save(fullfile(outDir,'circular_dromion_numerical_data.mat'), ...
    'X','Y','t_save','snapshot_times','snapshot_density','E_I1','E_inf');

toc
