function contributionPlot_realrate(Year, Rshort_X_50, Rshort_bar_X_50, Rshort_bar_X_idio_50, colorScheme, title_label, legend_label, outfile)

Area                 = [Rshort_X_50, Rshort_bar_X_idio_50];
posArea              = Area;
negArea              = Area;
posArea(posArea < 0) = 0;
negArea(negArea > 0) = 0;


f = figure;
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


p1 = plot(Year, Rshort_bar_X_50, 'Color', 'k', 'LineWidth', 2);



hline       = refline(0);
hline.Color = 'k';

xlim([Year(1), Year(end)])
%ylabel('Percent')
%xlabel('Year')
legend([a1, p1], legend_label, ...
                 'Location', 'southoutside',...
                  'Orientation', 'horizontal',...
                  'Interpreter', 'latex', ...
                  'boxoff')
legend boxoff
title(title_label, 'Interpreter', 'latex')

yl = get(gca, 'Ylim');
ylim(yl * 2)

printpdf(gcf, outfile)
end