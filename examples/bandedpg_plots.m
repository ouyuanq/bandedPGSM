%% Example 1
A = readmatrix("ex1_time.txt");
n = A(:, 1);

figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 2), '-ok', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 3), '-*k', 'LineWidth', 1, 'MarkerSize', 8)
loglog(n, A(:, 4), '-+k', 'LineWidth', 1, 'MarkerSize', 8)

legend('banded PG', 'recurrence', 'G-NI', 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 12)
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('construction time (sec)', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)

loglog(n(end-3:end), n(end-3:end) .^ 2 / n(end) ^ 2 * (A(end, 4) * 1.2), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
loglog(n(end-5:end), n(end-5:end) / n(end) * (A(end, 3) / 1.5), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
loglog(n(end-5:end), n(end-5:end) / n(end) * (A(end, 2) / 1.5), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
text(1e3, 1.5e-1, '$\mathcal{O}(n^2)$', 'Interpreter', 'latex', 'FontSize', 12)
text(1.5e3, 5e-3, '$\mathcal{O}(n)$', 'Interpreter', 'latex', 'FontSize', 12)
text(2e3, 3e-4, '$\mathcal{O}(n)$', 'Interpreter', 'latex', 'FontSize', 12)

xlim([n(1) / 1.1, n(end)*1.1])
ylim([min(A(:, 2:end), [], 'all') / 1.5, max(A(:, 2:end), [], 'all') * 2])

exportgraphics(gcf, 'ex1_construction.png', 'Resolution', 200)
%%
A = readmatrix("ex1_accuracy.txt");
n = A(:, 1);

figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 2), '-ok', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 3), '-*k', 'LineWidth', 1, 'MarkerSize', 8)
loglog(n, A(:, 4), '-+k', 'LineWidth', 1, 'MarkerSize', 8)

legend('banded PG', 'recurrence', 'G-NI', 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 12)
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$L^2$ error', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)


xlim([n(1) / 1.1, n(end)*1.1])
ylim([min(A(:, 2:end), [], 'all') / 1.5, max(A(:, 2:end), [], 'all') * 2])

exportgraphics(gcf, 'ex1_accuracy.png')
%% Example 2
A = readmatrix("ex2_time.txt");
n = A(:, 1);

figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 2), '-ok', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 3), '-*k', 'LineWidth', 1, 'MarkerSize', 8)
loglog(n, A(:, 4), '-+k', 'LineWidth', 1, 'MarkerSize', 8)
loglog(n, A(:, 5), '-sk', 'LineWidth', 1, 'MarkerSize', 8)

legend('banded PG', 'recurrence', 'G-NI', 'US', 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 12)
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('construction time (sec)', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)

loglog(n(end-3:end), n(end-3:end) .^ 2 / n(end) ^ 2 * (A(end, 4) * 2), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
loglog(n(end-5:end), n(end-5:end) / n(end) * (A(end, 2) / 2), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
text(1e4, 5e1, '$\mathcal{O}(n^2)$', 'Interpreter', 'latex', 'FontSize', 12)
text(3e4, 2e-3, '$\mathcal{O}(n)$', 'Interpreter', 'latex', 'FontSize', 12)

xlim([n(1) / 1.1, n(end)*1.1])
ylim([min(A(:, 2:end), [], 'all') / 1.5, max(A(:, 2:end), [], 'all') * 2])

exportgraphics(gcf, 'ex2_time.png', 'Resolution', 200)
%%
A = readmatrix("ex2_accuracy.txt");
n = A(:, 1);

figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 2), '-ok', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 3), '-*k', 'LineWidth', 1, 'MarkerSize', 8)
loglog(n, A(:, 4), '-+k', 'LineWidth', 1, 'MarkerSize', 8)
loglog(n, A(:, 5), '-sk', 'LineWidth', 1, 'MarkerSize', 8)

legend('banded PG', 'recurrence', 'G-NI', 'US', 'Interpreter', 'latex', 'Location', 'southwest', 'FontSize', 12)
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$L^2$ error', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)

xlim([10 1e5])
ylim([1e-13, max(A(:, 2:end), [], 'all') * 5])
yticks(10 .^ (-13:4:-1))

exportgraphics(gcf, 'ex2_accuracy.png', 'Resolution', 300)

%% Example 3
A = readmatrix("ex3_time.txt");
n = A(:, 1);

figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 2), '-*k', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 3), '-ok', 'LineWidth', 1, 'MarkerSize', 8)

legend('original US', 'accelerated US', 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 12)
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('construction time (sec)', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)

loglog(n(end-4:end), n(end-4:end) / n(end) * (A(end, 3) / 1.7), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
text(4e5, 3e-1, '$\mathcal{O}(n)$', 'Interpreter', 'latex', 'FontSize', 12)

xlim([n(1) / 1.1, n(end)*1.1])
ylim([min(A(:, 2:end), [], 'all') / 1.5, max(A(:, 2:end), [], 'all') * 2])

exportgraphics(gcf, 'ex3_con.png', 'Resolution', 200)

figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 4), '-*k', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 5), '-ok', 'LineWidth', 1, 'MarkerSize', 8)

xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('solving time (sec)', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)

loglog(n(end-4:end), n(end-4:end) / n(end) * (A(end, 5) / 1.5), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
text(4e5, 7e-2, '$\mathcal{O}(n)$', 'Interpreter', 'latex', 'FontSize', 12)

xlim([n(1) / 1.1, n(end)*1.1])
ylim([min(A(:, 2:end), [], 'all') / 1.5, max(A(:, 2:end), [], 'all') * 2])

exportgraphics(gcf, 'ex3_sol.png', 'Resolution', 200)