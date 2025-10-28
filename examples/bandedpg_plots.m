%% construction of multiplication matrix
A = readmatrix("Multiplication_variable_C10.txt");
n = A(:, 1);

figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 3), '-ok', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 2), '-*k', 'LineWidth', 1, 'MarkerSize', 8)

legend('new', '3-term recurrence', 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 12)
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('execution time (sec)', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)

loglog(n(end-4:end), (n(end-4:end) / n(end)) .^ 1 * (A(end, 3) / 1.5), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
loglog(n(end-4:end), (n(end-4:end) / n(end)) .^ 2 * (A(end, 2) * 1.5), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
text(1e3, 1e-1, '$\mathcal{O}(m)$', 'Interpreter', 'latex', 'FontSize', 12)
text(5e2, 3e2, '$\mathcal{O}(m^2)$', 'Interpreter', 'latex', 'FontSize', 12)

xlim([25, n(end)*1.1])
ylim([A(1, 3) / 1.5, A(end, 2) * 2])
yticks(10 .^ (-2:3))

exportgraphics(gcf, 'Multiplication_variable_C10.png')

A = readmatrix("Multiplication_low_C2.txt");
n = A(:, 1);

figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 3), '-ok', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 2), '-*k', 'LineWidth', 1, 'MarkerSize', 8)

% legend('new', '3-term recurrence', 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 12)
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('execution time (sec)', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)

loglog(n(end-5:end), (n(end-5:end) / n(end)) .^ 1 * (A(end, 3) / 1.8), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
loglog(n(end-5:end), (n(end-5:end) / n(end)) .^ 1 * (A(end, 2) * 1.5), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
text(5e5, 1.5e-2, '$\mathcal{O}(n)$', 'Interpreter', 'latex', 'FontSize', 12)
text(1e5, 1e1, '$\mathcal{O}(n)$', 'Interpreter', 'latex', 'FontSize', 12)

xlim([25, n(end)*1.1])
ylim([A(1, 3) / 1.5, A(end, 2) * 2])
yticks(10 .^ (-5:1))

exportgraphics(gcf, 'Multiplication_low_C2.png')

%% third order ODE with variable coefficients
A = readmatrix("taylor_time.txt");
n = A(:, 1);

figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 2), '-ok', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 3), '-*k', 'LineWidth', 1, 'MarkerSize', 8)
loglog(n, A(:, 4), '-+k', 'LineWidth', 1, 'MarkerSize', 8)

legend('new', 'MPG(R)', 'MPG(NI)', 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 12)
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('construction time (sec)', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)

loglog(n(end-2:end), n(end-2:end) .^ 3 / n(end) ^ 3 * (A(end, 4) / 1.5), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
loglog(n(end-3:end), n(end-3:end) / n(end) * (A(end, 3) / 1.5), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
loglog(n(end-3:end), n(end-3:end) / n(end) * (A(end, 2) / 1.5), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
text(4e3, 5e-1, '$\mathcal{O}(n^3)$', 'Interpreter', 'latex', 'FontSize', 12)
text(3e3, 6e-4, '$\mathcal{O}(n)$', 'Interpreter', 'latex', 'FontSize', 12)
% text(3.5e3, 2e-1, '$\mathcal{O}(n)$', 'Interpreter', 'latex', 'FontSize', 12)

xlim([n(1) / 1.1, n(end)*1.1])
ylim([min(A(:, 2:end), [], 'all') / 1.5, max(A(:, 2:end), [], 'all') * 2])

exportgraphics(gcf, 'taylor_time.png', 'Resolution', 200)

%%
A = readmatrix("taylor_accuracy.txt");
n = A(:, 1);

figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 2), '-ok', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 3), '-*k', 'LineWidth', 1, 'MarkerSize', 8)
loglog(n, A(:, 4), '-+k', 'LineWidth', 1, 'MarkerSize', 8)
% loglog(n, A(:, 5), '-sk', 'LineWidth', 1, 'MarkerSize', 8)

% legend('new(dual)', 'MPG(R)', 'MPG(NI)', 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 12)
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$L^2$ error', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)


xlim([n(1) / 1.1, n(end)*1.1])
ylim([min(A(:, 2:4), [], 'all') / 2, max(A(:, 2:4), [], 'all') * 2])

exportgraphics(gcf, 'taylor_err.png')

%% Airy equation
A = readmatrix("airy_time.txt");
n = A(:, 1);

figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 2), '-ok', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 3), '-*k', 'LineWidth', 1, 'MarkerSize', 8)
loglog(n, A(:, 4), '-+k', 'LineWidth', 1, 'MarkerSize', 8)
loglog(n, A(:, 5), '-sk', 'LineWidth', 1, 'MarkerSize', 8)

legend('new', 'MPG(R)', 'MPG(NI)', 'US', 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 12)
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('execution time (sec)', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)

A4 = A(:, 4); A4 = A4(A4 ~= 0); lenA4 = length(A4);
loglog(n(lenA4-3:lenA4), n(lenA4-3:lenA4) .^ 3 / n(lenA4) ^ 3 * (A4(end) / 2), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
loglog(n(end-4:end), n(end-4:end) / n(end) * (A(end, 2) / 1.5), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
text(3e4, 2e-3, '$\mathcal{O}(n)$', 'Interpreter', 'latex', 'FontSize', 12)
text(2e4, 6e1, '$\mathcal{O}(n^3)$', 'Interpreter', 'latex', 'FontSize', 12)

xlim([n(1) / 1.2, n(end)*1.2])
ylim([min(A(:, 2), [], 'all') / 2, max(A(:, 2:end), [], 'all') * 2])
yticks([1e-4 1e-2 1e0 1e2 1e4])

exportgraphics(gcf, 'airy_time.png', 'Resolution', 200)

%%
A = readmatrix("airy_accuracy.txt");
n = A(:, 1);

figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 2), '-ok', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 3), '-*k', 'LineWidth', 1, 'MarkerSize', 8)
loglog(n, A(:, 4), '-+k', 'LineWidth', 1, 'MarkerSize', 8)
loglog(n, A(:, 5), '-sk', 'LineWidth', 1, 'MarkerSize', 8)

% legend('new', 'MPG(R)', 'MPG(NI)', 'US', 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 12)
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$L^2$ error', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)

xlim([10 1e5])
ylim([1e-13, max(A(:, 2:end), [], 'all') * 5])
yticks(10 .^ (-13:4:-1))

exportgraphics(gcf, 'airy_err.png', 'Resolution', 300)

%% fourth order ODE with non-polynomial coefficient
A = readmatrix("oscillatory_time.txt");
n = A(:, 1);

figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 2), '-ok', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 3), '-*k', 'LineWidth', 1, 'MarkerSize', 8)
loglog(n, A(:, 4), '-+k', 'LineWidth', 1, 'MarkerSize', 8)
legend('new', 'MPG(R)', 'MPG(NI)', 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 12)

xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('construction time (sec)', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)

loglog(n(end-3:end), n(end-3:end) / n(end) * (A(end, 2) / 1.5), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
text(3e3, 3e-3, '$\mathcal{O}(n)$', 'Interpreter', 'latex', 'FontSize', 12)
loglog(n(end-3:end), n(end-3:end) .^ 3 / n(end)^3 * (A(end, 4) / 1.5), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
text(2e3, 6e-2, '$\mathcal{O}(n^3)$', 'Interpreter', 'latex', 'FontSize', 12)

xlim([n(1) / 1.1, n(end)*1.1])
ylim([min(A(:, 2:end), [], 'all') / 1.5, max(A(:, 2:4), [], 'all') * 8])

exportgraphics(gcf, 'oscillatory_time.png', 'Resolution', 200)

%%
A = readmatrix("oscillatory_accuracy.txt");
n = A(:, 1);

figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 2), '-ok', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 3), '-*k', 'LineWidth', 1, 'MarkerSize', 8)
loglog(n, A(:, 4), '-+k', 'LineWidth', 1, 'MarkerSize', 8)

% legend('new', 'MPG(R)', 'MPG(NI)', 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 12)
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$L^2$ error', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)

xlim([n(1) / 1.1, n(end)*1.1])
ylim([min(A(:, 2:end), [], 'all') / 4, max(A(:, 2:end), [], 'all') * 4])

exportgraphics(gcf, 'oscillatory_err.png', 'Resolution', 200)

%% tenth order equation
A = readmatrix("tenth_time.txt");
n = A(:, 1);

figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 3), '-ok', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 2), '-*k', 'LineWidth', 1, 'MarkerSize', 8)

legend('US (accelerated)', 'US (original)', 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', 12)
xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('construction time (sec)', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)

loglog(n(end-4:end), n(end-4:end) / n(end) * (A(end, 3) / 1.7), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
text(4e5, 3e-1, '$\mathcal{O}(n)$', 'Interpreter', 'latex', 'FontSize', 12)

xlim([n(1) / 1.1, n(end)*1.1])
ylim([min(A(:, 2:end), [], 'all') / 1.5, max(A(:, 2:end), [], 'all') * 2])

exportgraphics(gcf, 'tenth_con.png', 'Resolution', 200)
%%
figure
set(gcf, 'Position', [200 200 600 350])
loglog(n, A(:, 4), '-*k', 'LineWidth', 1, 'MarkerSize', 8)
hold on
loglog(n, A(:, 5), '-ok', 'LineWidth', 1, 'MarkerSize', 8)

xlabel('$n$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('solution time (sec)', 'Interpreter', 'latex', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)

loglog(n(end-4:end), n(end-4:end) / n(end) * (A(end, 5) / 1.5), '--k', 'LineWidth', 1.2, 'HandleVisibility', 'off')
text(4e5, 7e-2, '$\mathcal{O}(n)$', 'Interpreter', 'latex', 'FontSize', 12)

xlim([n(1) / 1.1, n(end)*1.1])
ylim([min(A(:, 2:end), [], 'all') / 1.5, max(A(:, 2:end), [], 'all') * 2])

exportgraphics(gcf, 'tenth_sol.png', 'Resolution', 200)