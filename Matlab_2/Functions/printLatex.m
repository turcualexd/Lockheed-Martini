function printLatex(L)

figure('menubar','none') ;
t = text(1, 1, ['$$' latex(L) '$$'], 'interpreter', 'latex', 'FontSize', 20);
set(gca, 'visible', 'off', 'xlim', [0 2], 'ylim', [0 2], 'Position', [0 0 1 1]);
set(t, 'visible', 'on', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');