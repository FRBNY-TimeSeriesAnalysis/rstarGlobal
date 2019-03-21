function contributionPlot_baanom(Year, M_bar, Ts_bar, Ts_bar_idio, ...
    ql_Pi_bar_us, qPi_bar_us_idio, colorScheme, title_label, legend_label, outfile)

% Add global and idiosyncratic inflation terms.

Area                 = [M_bar, Ts_bar, Ts_bar_idio, ql_Pi_bar_us qPi_bar_us_idio];
posArea              = Area;
negArea              = Area;
posArea(posArea < 0) = 0;
negArea(negArea > 0) = 0;


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


p1 = plot(Year, M_bar + Ts_bar + Ts_bar_idio + ql_Pi_bar_us + qPi_bar_us_idio, 'Color', 'k', 'LineWidth', 2);



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
printpdf(gcf, outfile)
end