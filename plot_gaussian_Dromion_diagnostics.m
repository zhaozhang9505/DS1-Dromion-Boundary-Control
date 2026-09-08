function [figure_file,peak_file,mass_file] = plot_gaussian_Dromion_diagnostics(data_file,figure_visibility,export_standalone)
% Replot Figure 5 from saved numerical data: a two-by-two summary and standalone P(t) and M(t) curves.
% This function performs no time stepping and does not replace numerical snapshots with the analytical solution.
% When the third argument is false, only the summary is replotted; the standalone curves are left unchanged.
if nargin < 1 || isempty(data_file)
    data_file = fullfile(fileparts(mfilename('fullpath')), ...
        'gaussian_dromion_smooth_results','gaussian_dromion_smooth_damped_data.mat');
end
if nargin < 2
    figure_visibility = 'on';
end
if nargin < 3
    export_standalone = true;
end
d = load(data_file);
outDir = fileparts(data_file);
required_fields = {'M','P','snapshot_steps','Omega','t0','t_end','is_short_test'};
if ~all(isfield(d,required_fields))
    error('旧数据缺少所需字段或八时刻快照，请先运行更新后的主程序。');
end
figure_file = '';
peak_file = '';
mass_file = '';
if ~d.is_short_test
    period = 2*pi/d.Omega;
    expected_times = [d.t0+(0:3)*period/4,d.t_end-(3:-1:0)*period/4];
    if numel(d.snapshot_times) ~= 8 || ...
            max(abs(d.snapshot_times(:).'-expected_times)) > 1e-10
        error('快照时间与第一周期、最后周期的指定时刻不一致，不能生成正式总图。');
    end
    % Use the same text width, column positions, and square plotting-area dimensions as Figure 4.
    figure_width = 16.6;
    panel_width = 6.45;
    panel_height = 6.45;
    column_gap = 2.00;
    left_margin = 1.20;
    bottom_margin = 1.60;
    row_gap = 1.80;
    top_margin = 0.75;
    upper_row = bottom_margin+panel_height+row_gap;
    right_column = left_margin+panel_width+column_gap;
    figure_height = upper_row+panel_height+top_margin;
    fig = figure('Color','w','Units','centimeters', ...
        'Position',[2,2,figure_width,figure_height],'Visible',figure_visibility);
    positions = [left_margin upper_row panel_width panel_height; ...
                 right_column upper_row panel_width panel_height; ...
                 left_margin bottom_margin panel_width panel_height; ...
                 right_column bottom_margin panel_width panel_height];
    ax_first = axes(fig,'Units','centimeters','Position',positions(1,:));
    draw_contours(ax_first,d,1:4,'(a) First cycle');
    ax_last = axes(fig,'Units','centimeters','Position',positions(2,:));
    draw_contours(ax_last,d,5:8,'(b) Last cycle');
    ax_peak = axes(fig,'Units','centimeters','Position',positions(3,:));
    draw_history(ax_peak,d,d.P,'$P(t)$','(c) Peak intensity',[0.85 0.10 0.10]);
    ax_mass = axes(fig,'Units','centimeters','Position',positions(4,:));
    draw_history(ax_mass,d,d.M,'$M(t)$','(d) Mass',[0.00 0.20 0.65]);
    panel_axes = [ax_first ax_last ax_peak ax_mass];
    for panel_index = 1:4
        ax = panel_axes(panel_index);
        set(ax,'FontSize',10.5,'LineWidth',0.95);
        set(get(ax,'XLabel'),'FontSize',13.0);
        set(get(ax,'YLabel'),'FontSize',13.0);
        set(get(ax,'Title'),'FontSize',12.0);
        if panel_index >= 3
            ax.XAxis.FontSize = 9.0;
            set(findobj(ax,'Type','line'),'LineWidth',1.15);
        end
        pbaspect(ax,[1 1 1]);
        set(ax,'Position',positions(panel_index,:));
    end
    drawnow;
    for panel_index = 1:4
        actual_position = get(panel_axes(panel_index),'Position');
        assert(max(abs(actual_position(3:4)-[panel_width panel_height])) < 1e-8, ...
            'Each Figure 5 plot box must match the Figure 4 plot box.');
    end
    contour_handles = [findall(ax_first,'Type','contour');findall(ax_last,'Type','contour')];
    assert(numel(contour_handles) == 8,'Expected eight numerical contour snapshots.');
    assert(all(abs(cell2mat(get(contour_handles,'LineWidth'))-1.10) < 1e-8), ...
        'The top-row contour line widths must match Figure 4.');
    figure_file = fullfile(outDir,'gaussian_dromion_smooth_damped_evolution_600dpi.png');
    set(fig,'PaperUnits','centimeters', ...
        'PaperPosition',[0,0,figure_width,figure_height], ...
        'PaperSize',[figure_width,figure_height]);
    print(fig,figure_file,'-dpng','-r600');
    fprintf('Figure 5: %.2f x %.2f cm; four plot boxes: %.2f x %.2f cm; contour width: 1.10 pt.\n', ...
        figure_width,figure_height,panel_width,panel_height);
    fprintf('Saved: %s\n',figure_file);
    if strcmpi(figure_visibility,'off'), close(fig); end
end

if ~export_standalone
    return;
end

fig_peak = figure('Color','w','Units','centimeters', ...
    'Position',[2,2,13.5,10.5],'Visible',figure_visibility);
ax_peak = axes(fig_peak,'Position',[0.15 0.23 0.80 0.66]);
draw_history(ax_peak,d,d.P,'$P(t)$','Peak intensity versus time',[0.85 0.10 0.10]);
peak_file = fullfile(outDir,'gaussian_peak_intensity_history_600dpi.png');
set(fig_peak,'PaperUnits','centimeters','PaperPosition',[0,0,13.5,10.5], ...
    'PaperSize',[13.5,10.5]);
print(fig_peak,peak_file,'-dpng','-r600');

fig_mass = figure('Color','w','Units','centimeters', ...
    'Position',[2,2,13.5,10.5],'Visible',figure_visibility);
ax_mass = axes(fig_mass,'Position',[0.15 0.23 0.80 0.66]);
draw_history(ax_mass,d,d.M,'$M(t)$','Mass versus time',[0.00 0.20 0.65]);
mass_file = fullfile(outDir,'gaussian_integrated_intensity_history_600dpi.png');
set(fig_mass,'PaperUnits','centimeters','PaperPosition',[0,0,13.5,10.5], ...
    'PaperSize',[13.5,10.5]);
print(fig_mass,mass_file,'-dpng','-r600');
if strcmpi(figure_visibility,'off'), close(fig_peak); close(fig_mass); end
end

function draw_contours(ax,d,indices,panel_title)
    % Match the colors, contour line widths, and annotation fonts of Figure 4, retaining the numerical snapshots and intensity levels.
    colors = [0 0 0; 0.85 0.10 0.10; 0 0.35 0.85; 0 0.55 0.20];
    levels = [0.20 0.50 0.80];
    position_by_quadrant = [10.8 3.0; 1.8 11.0; -10.8 -3.0; 1.8 -11.0];
    hold(ax,'on');
    orbit_angle = linspace(0,2*pi,600);
    plot(ax,d.R*cos(orbit_angle),d.R*sin(orbit_angle),'--', ...
        'Color',[0.62 0.62 0.62],'LineWidth',0.70,'HandleVisibility','off');
    % Append the right and upper periodic endpoints to avoid spline extrapolation at the plot boundaries.
    X_closed = [d.X;d.x1]; Y_closed = [d.Y;d.y1];
    [xp,yp] = meshgrid(linspace(d.x0,d.x1,401),linspace(d.y0,d.y1,401));
    legend_handles = gobjects(1,4);
    labels = cell(1,4);
    for j = 1:4
        js = indices(j);
        z = d.snapshot_density(:,:,js);
        z_closed = [z,z(:,1);z(1,:),z(1,1)];
        z_plot = max(real(interp2(X_closed,Y_closed,z_closed,xp,yp,'spline')),0);
        contour(ax,xp,yp,z_plot,levels,'LineColor',colors(j,:), ...
            'LineWidth',1.10,'HandleVisibility','off');
        legend_handles(j) = plot(ax,nan,nan,'-','Color',colors(j,:),'LineWidth',1.10);
        labels{j} = pi_label(d.snapshot_times(js),true);
        quadrant = mod(round(d.snapshot_times(js)*d.Omega/(pi/2)),4)+1;
        pos = position_by_quadrant(quadrant,:);
        text(ax,pos(1),pos(2),labels{j},'Interpreter','latex', ...
            'Color',colors(j,:),'FontName','Times New Roman','FontSize',10.0, ...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    axis(ax,'equal');
    xlim(ax,[d.x0 d.x1]); ylim(ax,[d.y0 d.y1]);
    xticks(ax,[d.x0 -8 0 8 d.x1]); yticks(ax,[d.y0 -8 0 8 d.y1]);
    style_axes(ax);
    xlabel(ax,'$x$','Interpreter','latex','FontSize',20.5);
    ylabel(ax,'$y$','Interpreter','latex','FontSize',20.5);
    title(ax,panel_title,'Interpreter','latex','FontSize',18,'FontWeight','normal');
    legend(ax,legend_handles,labels,'Interpreter','latex', ...
        'Location','northwest','FontSize',9.5,'Box','on');
end

function draw_history(ax,d,values,y_label,panel_title,color)
    plot(ax,d.t_save,values,'-','Color',color,'LineWidth',1.75);
    xlim(ax,[d.t0 d.t_end]);
    style_axes(ax);
    if d.is_short_test
        xticks(ax,linspace(d.t0,d.t_end,6));
    else
        period = 2*pi/d.Omega;
        ticks = d.t0:period:d.t_end+1e-10;
        xticks(ax,ticks);
        tick_labels = arrayfun(@(t) pi_label(t,false),ticks,'UniformOutput',false);
        xticklabels(ax,tick_labels);
        xtickangle(ax,45);
        ax.XAxis.FontSize = 14;
    end
    ylim(ax,[0,1.08*max(values)]);
    xlabel(ax,'$t$','Interpreter','latex','FontSize',20.5);
    ylabel(ax,y_label,'Interpreter','latex','FontSize',20.5);
    title(ax,panel_title,'Interpreter','latex','FontSize',18,'FontWeight','normal');
end

function style_axes(ax)
    set(ax,'FontName','Times New Roman','FontSize',16.5, ...
        'TickLabelInterpreter','latex','LineWidth',1.40,'TickDir','out', ...
        'Layer','top','Box','on');
    grid(ax,'on');
    set(ax,'GridColor',[0.85 0.85 0.85],'GridAlpha',0.55);
end

function label = pi_label(t,include_time)
    k = t/pi;
    if abs(k) < 1e-10
        value = '0';
    elseif abs(k-round(k)) < 1e-9
        value = sprintf('%d\\pi',round(k));
    else
        value = sprintf('%.4g\\pi',k);
    end
    if include_time
        label = ['$t=',value,'$'];
    else
        label = ['$',value,'$'];
    end
end
