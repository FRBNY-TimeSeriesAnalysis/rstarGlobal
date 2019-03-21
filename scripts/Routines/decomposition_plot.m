function decomposition_plot(Year, Data, colorScheme, title_label, legend_label, axis_scale)

% Add global and idiosyncratic inflation terms.

Area                 = Data;
posArea              = Area;
negArea              = Area;
posArea(posArea < 0) = 0;
negArea(negArea > 0) = 0;


%f = figure;
box on;
hold on;

% Plot areas
a1              = area(Year, posArea);
a2              = area(Year, negArea);

for iPlot = 1:length(a1)  % Display settings
    a1(iPlot).FaceColor = colorScheme(iPlot, :);
    a2(iPlot).FaceColor = colorScheme(iPlot, :);
    a1(iPlot).LineStyle    = 'none';
    a2(iPlot).LineStyle    = 'none';
end


p1 = plot(Year, sum(Data, 2), 'Color', 'k', 'LineWidth', 2);



hline       = refline(0);
hline.Color = 'k';

xlim([Year(1), Year(end)])
%ylabel('Percent')
%xlabel('Year')

if numel(legend_label) > 0
    legend([a1, p1], legend_label, ...
        'Location', 'southoutside',...
        'Orientation', 'horizontal',...
        'Interpreter', 'latex', ...
        'boxoff')
    legend boxoff
end

if numel(title_label) > 0
    title(title_label, 'Interpreter', 'latex')
end

yl = get(gca, 'Ylim');
ylim(yl * axis_scale)

%printpdf(gcf, outfile)
end