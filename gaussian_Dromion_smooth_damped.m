clc; clear; close all;
tic

%% 1. Computational domain, grid, and time parameters
x0 = -15;  x1 = 15;  Nx = 128;
y0 = -15;  y1 = 15;  Ny = 128;
t0 = 0;    t_end = 8*pi*10;
period_steps = 10000*10;
dt = (t_end-t0)/period_steps;

% Frame output interval: display and save one evolution plot every ND RK4 steps.
% Change ND to control the output frequency; the final frame is always saved, even if ND does not divide the total step count.
ND = 1000;
frame_dpi = 300;
% Figure windows can be hidden during automated background checks; normal runs display frames every ND steps.
figure_visibility = 'on';
if strcmpi(getenv('GAUSSIAN_FIGURE_VISIBLE'),'off')
    figure_visibility = 'off';
end

% Production runs use t_end as specified above. For automated checks only,
% GAUSSIAN_TEST_TEND can specify an earlier final time; the test step count
% is rounded down so that the actual final time does not exceed the requested value.
test_t_end_text = getenv('GAUSSIAN_TEST_TEND');
is_short_test = ~isempty(test_t_end_text);
if is_short_test
    requested_test_end = str2double(test_t_end_text);
    if ~isfinite(requested_test_end) || requested_test_end <= t0
        error('GAUSSIAN_TEST_TEND 必须大于初始时刻。');
    end
    Nt = max(1,floor((requested_test_end-t0)/dt));
    t_end = t0+Nt*dt;
    fprintf('短程测试模式：0 <= t <= %.8f，共 %d 步。\n',t_end,Nt);
else
    Nt = period_steps;
end

if ~isscalar(ND) || ~isnumeric(ND) || ~isfinite(ND) || ND < 1 || ND ~= round(ND)
    error('ND 必须是正整数。');
end
ND = double(ND);
if ~isscalar(frame_dpi) || ~isnumeric(frame_dpi) || ~isfinite(frame_dpi) ...
        || frame_dpi < 1 || frame_dpi ~= round(frame_dpi)
    error('frame_dpi 必须是正整数。');
end
diagnostic_step = 100;

Lx = x1-x0;  dx = Lx/Nx;
Ly = y1-y0;  dy = Ly/Ny;
X = (x0+(0:Nx-1)*dx).';
Y = (y0+(0:Ny-1)*dy).';
[x,y] = meshgrid(X,Y);
Kx = (2*pi/Lx)*[0:Nx/2-1 -Nx/2:-1].';
Ky = (2*pi/Ly)*[0:Ny/2-1 -Ny/2:-1].';
[kx,ky] = meshgrid(Kx,Ky);
k2 = kx.^2+ky.^2;
y_mid = Y(1:end-1)+dy/2;
x_mid = X(1:end-1)+dx/2;

%% 2. Mean flows for circular motion and ten-layer smooth quintic damping
R = 8;
Omega = 1/4;
c = (2+sqrt(2))/2;
n_damp = 10;
gamma_max = 4;
Gamma = zeros(Ny,Nx);
for layer = 1:n_damp
    s = (n_damp-layer+1)/n_damp;
    smooth_ramp = 10*s^3-15*s^4+6*s^5;
    gamma_layer = gamma_max*smooth_ramp;
    Gamma(layer,:) = max(Gamma(layer,:),gamma_layer);
    Gamma(end-layer+1,:) = max(Gamma(end-layer+1,:),gamma_layer);
    Gamma(:,layer) = max(Gamma(:,layer),gamma_layer);
    Gamma(:,end-layer+1) = max(Gamma(:,end-layer+1),gamma_layer);
end

%% 3. Gaussian initial condition for the main field
% The full-plane integral of |A_G|^2 is 8*log(2), and its peak is 2.
% Relative to the exact dromion, the mass is four times as large and the
% peak intensity is twice as large. The center and local carrier phase
% remain (8,0) and exp(i*y), respectively.
A0 = sqrt(2)*exp(-pi*((x-R).^2+y.^2)/(8*log(2))).*exp(1i*y);
fut = fft2(A0);

target_mass = 8*log(2);
discrete_mass = sum(abs(A0).^2,'all')*dx*dy;
fprintf('target mass = %.15g, discrete mass = %.15g\n', ...
    target_mass,discrete_mass);
fprintf('target peak = 2, grid peak = %.15g\n',max(abs(A0(:)).^2));

%% 4. Output times and diagnostics
if is_short_test
    snapshot_steps = [0,Nt];
else
    orbit_period = 2*pi/Omega;
    requested_snapshot_times = [t0+(0:3)*orbit_period/4, ...
        t_end-(3:-1:0)*orbit_period/4];
    snapshot_steps = round((requested_snapshot_times-t0)/dt);
    if any(snapshot_steps < 0) || any(snapshot_steps > Nt) ...
            || max(abs(t0+snapshot_steps*dt-requested_snapshot_times)) > 1e-10
        error('八个四分之一周期快照必须位于计算区间内并与时间网格对齐。');
    end
    if numel(unique(snapshot_steps)) ~= 8
        error('前后两组快照必须是八个不同的时刻。');
    end
end
snapshot_times = t0+snapshot_steps*dt;
snapshot_density = zeros(Ny,Nx,numel(snapshot_steps));
snapshot_found = false(size(snapshot_steps));
frame_steps = unique([0:ND:Nt,Nt]);

% Diagnostics are recorded every 100 steps by default. Include all frame updates,
% quarter-period snapshots, and period endpoints in the saved steps so that no required data are missed for any ND.
save_steps = unique([0:diagnostic_step:Nt,frame_steps,snapshot_steps,Nt]);
t_save = t0+save_steps*dt;
nSave = numel(t_save);

M = zeros(nSave,1);
P = zeros(nSave,1);
peak_x = zeros(nSave,1);
peak_y = zeros(nSave,1);

if is_short_test
    outDir = fullfile(fileparts(mfilename('fullpath')),'gaussian_dromion_smooth_test_results');
else
    outDir = fullfile(fileparts(mfilename('fullpath')),'gaussian_dromion_smooth_results');
end
if ~exist(outDir,'dir')
    mkdir(outDir);
end
frameDir = fullfile(outDir,'frames');
if ~exist(frameDir,'dir')
    mkdir(frameDir);
end

%% 5. Main-field display updated every ND steps
fig_live = figure('Color','w','Units','centimeters', ...
    'Position',[2,2,25,10.8],'Visible',figure_visibility);
orbit_angle = linspace(0,2*pi,600);

%% 6. Fourier pseudospectral discretization and RK4 time stepping
current_step = 0;
for si = 1:nSave
    target_step = save_steps(si);
    for step = current_step+1:target_step
        t_now = t0+(step-1)*dt;
        fut0 = fut;

        K1 = gaussian_rhs(fut0,t_now,k2,kx,ky,X,Y,dx,dy, ...
            x_mid,y_mid,R,Omega,c,Gamma);
        K2 = gaussian_rhs(fut0+dt*K1/2,t_now+dt/2,k2,kx,ky,X,Y,dx,dy, ...
            x_mid,y_mid,R,Omega,c,Gamma);
        K3 = gaussian_rhs(fut0+dt*K2/2,t_now+dt/2,k2,kx,ky,X,Y,dx,dy, ...
            x_mid,y_mid,R,Omega,c,Gamma);
        K4 = gaussian_rhs(fut0+dt*K3,t_now+dt,k2,kx,ky,X,Y,dx,dy, ...
            x_mid,y_mid,R,Omega,c,Gamma);
        fut = fut0+dt*(K1+2*K2+2*K3+K4)/6;
    end
    current_step = target_step;

    A_num = ifft2(fut);
    density = abs(A_num).^2;
    M(si) = sum(density,'all')*dx*dy;
    [P(si),linear_index] = max(density(:));
    [row_index,column_index] = ind2sub(size(density),linear_index);
    peak_x(si) = X(column_index);
    peak_y(si) = Y(row_index);

    snapshot_index = find(target_step==snapshot_steps,1);
    if ~isempty(snapshot_index)
        snapshot_density(:,:,snapshot_index) = density;
        snapshot_found(snapshot_index) = true;
    end

    % Update every ND time steps: the left panel shows the squared modulus in 3D,
    % and the right panel shows contours with a dashed circle of radius 8 centered at the origin.
    if any(target_step==frame_steps)
        clf(fig_live);
        tiledlayout(fig_live,1,2,'Padding','compact','TileSpacing','compact');
        current_t = t_save(si);
        color_max = max(density(:));

        ax_surface = nexttile;
        surf(ax_surface,x,y,density,'EdgeColor','none');
        shading(ax_surface,'interp');
        view(ax_surface,42,28);
        xlim(ax_surface,[x0 x1]); ylim(ax_surface,[y0 y1]);
        zlim(ax_surface,[0,1.05*color_max]);
        clim(ax_surface,[0,color_max]);
        xlabel(ax_surface,'$x$','Interpreter','latex');
        ylabel(ax_surface,'$y$','Interpreter','latex');
        zlabel(ax_surface,'$|A|^2$','Interpreter','latex');
        title(ax_surface,sprintf('(a) $|A|^2$, $t=%.4f$',current_t), ...
            'Interpreter','latex','FontWeight','normal');
        set(ax_surface,'FontName','Times New Roman','FontSize',14, ...
            'TickLabelInterpreter','latex','LineWidth',1.1,'Box','on');
        grid(ax_surface,'on');
        colorbar(ax_surface);

        ax_contour_live = nexttile;
        contourf(ax_contour_live,x,y,density,30,'LineStyle','none');
        hold(ax_contour_live,'on');
        contour(ax_contour_live,x,y,density,8,'LineColor','k','LineWidth',0.55);
        plot(ax_contour_live,R*cos(orbit_angle),R*sin(orbit_angle),'k--', ...
            'LineWidth',1.65);
        axis(ax_contour_live,'equal');
        xlim(ax_contour_live,[x0 x1]); ylim(ax_contour_live,[y0 y1]);
        clim(ax_contour_live,[0,color_max]);
        xlabel(ax_contour_live,'$x$','Interpreter','latex');
        ylabel(ax_contour_live,'$y$','Interpreter','latex');
        title(ax_contour_live,sprintf('(b) $|A|^2$ contours, $t=%.4f$',current_t), ...
            'Interpreter','latex','FontWeight','normal');
        set(ax_contour_live,'FontName','Times New Roman','FontSize',14, ...
            'TickLabelInterpreter','latex','LineWidth',1.1, ...
            'Layer','top','Box','on');
        grid(ax_contour_live,'on');
        colorbar(ax_contour_live);
        colormap(fig_live,turbo(256));
        drawnow;

        frame_file = fullfile(frameDir,sprintf( ...
            'gaussian_frame_step_%06d_t_%010.4f.png',target_step,current_t));
        set(fig_live,'PaperUnits','centimeters', ...
            'PaperPosition',[0,0,25,10.8],'PaperSize',[25,10.8]);
        print(fig_live,frame_file,'-dpng',sprintf('-r%d',frame_dpi));
    end

    fprintf('saved %3d/%3d, t=%7.3f, M/M(0)=%.8f, peak=%.8f, at (%+.3f,%+.3f)\n', ...
        si,nSave,t_save(si),M(si)/M(1),P(si),peak_x(si),peak_y(si));
end

if ~all(snapshot_found)
    error('Not all quarter-period snapshots were stored.');
end

%% 7. Save numerical data, then generate the two-by-two summary figure and two standalone curves
% M is the mass, and P is the peak intensity; legacy field names are retained only for compatibility with earlier data-reading scripts.
I1 = M;
peak_density = P;
data_file = fullfile(outDir,'gaussian_dromion_smooth_damped_data.mat');
save(data_file,'X','Y','t_save','snapshot_times','snapshot_steps', ...
    'snapshot_density','M','P','I1','peak_density','peak_x','peak_y', ...
    'Gamma','n_damp','gamma_max','ND','frame_dpi','frame_steps', ...
    'x0','x1','y0','y1','Nx','Ny','dx','dy','dt','Nt','t0','t_end', ...
    'R','Omega','diagnostic_step','is_short_test','A_num');
[figure_file,peak_file,mass_file] = plot_gaussian_Dromion_diagnostics( ...
    data_file,figure_visibility);
fprintf('final: M/M(0)=%.10f, P/P(0)=%.10f, peak=(%+.6f,%+.6f)\n', ...
    M(end)/M(1),P(end)/P(1),peak_x(end),peak_y(end));
if ~isempty(figure_file)
    fprintf('Saved: %s\n',figure_file);
end
fprintf('Saved: %s\n',peak_file);
fprintf('Saved: %s\n',mass_file);
toc

function K = gaussian_rhs(fut,t,k2,kx,ky,X,Y,dx,dy, ...
    x_mid,y_mid,R,Omega,c,Gamma)
    A = ifft2(fut);
    density = abs(A).^2;
    density_hat = fft2(density);
    density_x = real(ifft2(1i*kx.*density_hat));
    density_y = real(ifft2(1i*ky.*density_hat));

    U_bottom = (3+2*sqrt(2))*sech(c*(X.'-R*cos(Omega*t))-log(2)/4).^2 ...
        -(X.'/4)*cos(Omega*t);
    V_left = (3+2*sqrt(2))*sech(c*(Y-R*sin(Omega*t))-log(2)/4).^2 ...
        -(Y/4)*sin(Omega*t);

    density_x_mid = interp1(Y,density_x,y_mid,'spline');
    U_step = dy*(density_x(1:end-1,:)+4*density_x_mid+density_x(2:end,:))/6;
    U = [U_bottom;U_bottom+cumsum(U_step,1)];

    density_y_mid = interp1(X,density_y.',x_mid,'spline').';
    V_step = dx*(density_y(:,1:end-1)+4*density_y_mid+density_y(:,2:end))/6;
    V = [V_left,V_left+cumsum(V_step,2)];

    K = 1i*(-k2.*fut+fft2((U+V).*A))-fft2(Gamma.*A);
end
