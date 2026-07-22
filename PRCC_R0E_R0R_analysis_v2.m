function PRCC_R0E_R0R_analysis_v2
clc;
clear;
close all;
rng(1);
%% 
% SETTINGS
% Number of Latin Hypercube samples
N = 5000;

% Save figures automatically?
saveFigures = false;

% Figure export resolution
exportResolution = 600;

% GLOBAL PLOTTING STYLE
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesLineWidth', 1.0);
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultAxesBox', 'on');
set(groot, 'defaultAxesTickDir', 'out');

% Professional figure colors
colorR0E = [0.10 0.35 0.65];   % deep academic blue
colorR0R = [0.70 0.25 0.12];   % deep academic reddish-brown

% BASELINE PARAMETER VALUES
base.mu     = 0.02;

% Normalization condition Lambda = eta = mu
base.Lambda = base.mu;
base.eta    = base.mu;

base.beta   = 0.70;
base.gamma  = 0.25;
base.psi    = 0.30;
base.alpha  = 0.50;
base.omega  = 0.05;
base.sigma  = 0.40;

base.lambda = 0.25;
base.phi    = 0.60;
base.delta  = 0.15;
base.theta  = 0.10;

baseline_R0E = base.beta * ((base.omega + base.mu) + base.sigma * base.psi) / ...
               ((base.omega + base.psi + base.mu) * (base.gamma + base.mu));

baseline_R0R = base.lambda / (base.delta + base.mu);

fprintf('\n==============================================================\n');
fprintf('NORMALIZED THEORETICAL BASELINE VALUES\n');
fprintf('All rates are per normalized model time unit.\n');
fprintf('==============================================================\n');
fprintf('Lambda  = %.4f   (= mu)\n', base.Lambda);
fprintf('mu      = %.4f\n', base.mu);
fprintf('eta     = %.4f   (= mu)\n', base.eta);
fprintf('beta    = %.4f\n', base.beta);
fprintf('gamma   = %.4f\n', base.gamma);
fprintf('psi     = %.4f\n', base.psi);
fprintf('alpha   = %.4f\n', base.alpha);
fprintf('omega   = %.4f\n', base.omega);
fprintf('sigma   = %.4f\n', base.sigma);
fprintf('lambda  = %.4f\n', base.lambda);
fprintf('phi     = %.4f\n', base.phi);
fprintf('delta   = %.4f\n', base.delta);
fprintf('theta   = %.4f\n', base.theta);
fprintf('--------------------------------------------------------------\n');
fprintf('Baseline R0E = %.6f\n', baseline_R0E);
fprintf('Baseline R0R = %.6f\n', baseline_R0R);
fprintf('LHS sample size N = %d\n', N);
fprintf('==============================================================\n\n');

fprintf('FIGURE DESCRIPTION\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Figure 1: PRCC bar plot for R0E. Positive bars indicate parameters\n');
fprintf('          that increase the epidemic reproduction number, while\n');
fprintf('          negative bars indicate parameters that reduce it.\n\n');
fprintf('Figure 2: Scatter plots of sampled epidemic parameters against R0E.\n');
fprintf('          All panels use the same y-axis scale and deep blue markers.\n\n');
fprintf('Figure 3: PRCC bar plot for R0R. Positive bars indicate parameters\n');
fprintf('          that increase the rumor reproduction number, while\n');
fprintf('          negative bars indicate parameters that reduce it.\n\n');
fprintf('Figure 4: Scatter plots of sampled rumor parameters against R0R.\n');
fprintf('          All panels use the same y-axis scale and reddish-brown markers.\n');
fprintf('--------------------------------------------------------------\n\n');

%PRCC FOR EPIDEMIC REPRODUCTION NUMBER R0E

paramNamesE = {'\beta','\omega','\mu','\psi','\sigma','\gamma'};
plainNamesE = {'beta','omega','mu','psi','sigma','gamma'};

mE = length(paramNamesE);

lowE  = zeros(1,mE);
highE = zeros(1,mE);

% Parameter uncertainty ranges: 50% to 150% of normalized baseline
lowE(1)  = 0.5 * base.beta;   highE(1) = 1.5 * base.beta;    % beta
lowE(2)  = 0.5 * base.omega;  highE(2) = 1.5 * base.omega;   % omega
lowE(3)  = 0.5 * base.mu;     highE(3) = 1.5 * base.mu;      % mu
lowE(4)  = 0.5 * base.psi;    highE(4) = 1.5 * base.psi;     % psi
lowE(5)  = 0.5 * base.sigma;  highE(5) = 1.5 * base.sigma;   % sigma
lowE(6)  = 0.5 * base.gamma;  highE(6) = 1.5 * base.gamma;   % gamma

XE = manual_lhs(N, lowE, highE);

betaE  = XE(:,1);
omegaE = XE(:,2);
muE    = XE(:,3);
psiE   = XE(:,4);
sigmaE = XE(:,5);
gammaE = XE(:,6);

% Normalization inside every sample: Lambda = mu
LambdaE = muE; %#ok<NASGU>

R0E = betaE .* ((omegaE + muE) + sigmaE .* psiE) ./ ...
      ((omegaE + psiE + muE) .* (gammaE + muE));

[PRCC_E, pval_E] = compute_prcc(XE, R0E);

fprintf('\n==============================================================\n');
fprintf('PRCC ANALYSIS FOR EPIDEMIC REPRODUCTION NUMBER R0E\n');
fprintf('==============================================================\n');
fprintf('R0E = beta * ((omega+mu) + sigma*psi) / ((omega+psi+mu)*(gamma+mu))\n\n');

print_distribution_summary('R0E', R0E);

fprintf('\n%-10s   %-10s   %-12s   %-6s\n','Parameter','PRCC','p-value','Sig.');
fprintf('%s\n', repmat('-',1,50));
for i = 1:mE
    starText = significance_stars(pval_E(i));
    fprintf('%-10s   %+0.6f    %0.3e     %-6s\n', ...
        plainNamesE{i}, PRCC_E(i), pval_E(i), starText);
end

print_ranking('R0E', plainNamesE, PRCC_E, pval_E);

%Figure 1: PRCC bar plot for R0E
fig1 = figure('Color','w','Position',[100 100 780 540]);
plot_prcc_bar(PRCC_E, pval_E, paramNamesE);

title('Partial Rank Correlation Coefficients for R_{0E}', ...
    'Interpreter','tex', ...
    'FontName','Times New Roman', ...
    'FontWeight','bold');

ylabel('PRCC', 'FontName','Times New Roman');
apply_academic_axes(gca);

if saveFigures
    save_publication_figure(fig1, 'Figure_PRCC_R0E_bar.png', exportResolution);
end

%Figure 2: Scatter plots for R0E
fig2 = figure('Color','w','Position',[80 80 1120 720]);
ylE = padded_limits(R0E);

for i = 1:mE
    subplot(2,3,i);
    
    scatter(XE(:,i), R0E, 13, ...
        'MarkerFaceColor', colorR0E, ...
        'MarkerEdgeColor', colorR0E, ...
        'MarkerFaceAlpha', 0.35, ...
        'MarkerEdgeAlpha', 0.35);
    
    xlabel(paramNamesE{i}, ...
        'Interpreter','tex', ...
        'FontName','Times New Roman');
    
    ylabel('R_{0E}', ...
        'Interpreter','tex', ...
        'FontName','Times New Roman');
    
title(sprintf('PRCC = %+0.3f', PRCC_E(i)), ...
        'FontWeight','bold', ...
        'FontName','Times New Roman');
    
    ylim(ylE);
    grid on;
    grid minor;
    apply_academic_axes(gca);
end

add_super_title('Scatter Plots of Sampled Epidemic Parameters versus R_{0E}');

if saveFigures
    save_publication_figure(fig2, 'Figure_PRCC_R0E_scatter.png', exportResolution);
end

%PRCC FOR RUMOR REPRODUCTION NUMBER R0R

paramNamesR = {'\lambda','\delta','\mu'};
plainNamesR = {'lambda','delta','mu'};

mR = length(paramNamesR);

lowR  = zeros(1,mR);
highR = zeros(1,mR);

% Parameter uncertainty ranges: 50% to 150% of normalized baseline
lowR(1)  = 0.5 * base.lambda;   highR(1) = 1.5 * base.lambda;   % lambda
lowR(2)  = 0.5 * base.delta;    highR(2) = 1.5 * base.delta;    % delta
lowR(3)  = 0.5 * base.mu;       highR(3) = 1.5 * base.mu;       % mu

XR = manual_lhs(N, lowR, highR);

lambdaR = XR(:,1);
deltaR  = XR(:,2);
muR     = XR(:,3);

% Normalization inside every sample: eta = mu
etaR = muR; %#ok<NASGU>

R0R = lambdaR ./ (deltaR + muR);

[PRCC_R, pval_R] = compute_prcc(XR, R0R);

fprintf('\n==============================================================\n');
fprintf('PRCC ANALYSIS FOR RUMOR REPRODUCTION NUMBER R0R\n');
fprintf('==============================================================\n');
fprintf('R0R = lambda / (delta + mu)\n\n');

print_distribution_summary('R0R', R0R);

fprintf('\n%-10s   %-10s   %-12s   %-6s\n','Parameter','PRCC','p-value','Sig.');
fprintf('%s\n', repmat('-',1,50));
for i = 1:mR
    starText = significance_stars(pval_R(i));
    fprintf('%-10s   %+0.6f    %0.3e     %-6s\n', ...
        plainNamesR{i}, PRCC_R(i), pval_R(i), starText);
end

print_ranking('R0R', plainNamesR, PRCC_R, pval_R);

%Figure 3: PRCC bar plot for R0R
fig3 = figure('Color','w','Position',[100 100 680 520]);
plot_prcc_bar(PRCC_R, pval_R, paramNamesR);

title('Partial Rank Correlation Coefficients for R_{0R}', ...
    'Interpreter','tex', ...
    'FontName','Times New Roman', ...
    'FontWeight','bold');

ylabel('PRCC', 'FontName','Times New Roman');
apply_academic_axes(gca);

if saveFigures
    save_publication_figure(fig3, 'Figure_PRCC_R0R_bar.png', exportResolution);
end

% Figure 4: Scatter plots for R0R
fig4 = figure('Color','w','Position',[90 90 1000 400]);
ylR = padded_limits(R0R);

for i = 1:mR
    subplot(1,3,i);
    
    scatter(XR(:,i), R0R, 15, ...
        'MarkerFaceColor', colorR0R, ...
        'MarkerEdgeColor', colorR0R, ...
        'MarkerFaceAlpha', 0.40, ...
        'MarkerEdgeAlpha', 0.40);
    
    xlabel(paramNamesR{i}, ...
        'Interpreter','tex', ...
        'FontName','Times New Roman');
    
    ylabel('R_{0R}', ...
        'Interpreter','tex', ...
        'FontName','Times New Roman');
    
 title(sprintf('PRCC = %+0.3f', PRCC_R(i)), ...
        'FontWeight','bold', ...
        'FontName','Times New Roman');
    
    ylim(ylR);
    grid on;
    grid minor;
    apply_academic_axes(gca);
end

add_super_title('Scatter Plots of Sampled Rumor Parameters versus R_{0R}');

if saveFigures
    save_publication_figure(fig4, 'Figure_PRCC_R0R_scatter.png', exportResolution);
end

fprintf('\n==============================================================\n');
fprintf('ANALYSIS COMPLETE\n');
if saveFigures
    fprintf('Generated figures:\n');
    fprintf('  Figure_PRCC_R0E_bar.png\n');
    fprintf('  Figure_PRCC_R0E_scatter.png\n');
    fprintf('  Figure_PRCC_R0R_bar.png\n');
    fprintf('  Figure_PRCC_R0R_scatter.png\n');
else
    fprintf('Figures were generated but not saved. Set saveFigures = true to export them.\n');
end
fprintf('==============================================================\n\n');

end

% Manual Latin Hypercube Sampling

function X = manual_lhs(N, low, high)
m = length(low);
X = zeros(N,m);

for j = 1:m
    perm = randperm(N)';
    u = rand(N,1);
    lhs_col = (perm - u) / N;
    X(:,j) = low(j) + lhs_col .* (high(j) - low(j));
end

end

% PRCC calculation

function [PRCC, pval] = compute_prcc(X, Y)
[N, m] = size(X);

RankX = zeros(N,m);
for j = 1:m
    RankX(:,j) = rank_vector_avg_ties(X(:,j));
end

RankY = rank_vector_avg_ties(Y);

PRCC = zeros(m,1);
pval = zeros(m,1);

for i = 1:m
    otherIndex = setdiff(1:m,i);
    Z = [ones(N,1), RankX(:,otherIndex)];

    bx = Z \ RankX(:,i);
    residualX = RankX(:,i) - Z * bx;

    by = Z \ RankY;
    residualY = RankY - Z * by;

    C = corrcoef(residualX, residualY);
    r = C(1,2);
    PRCC(i) = r;

    df = N - m - 1;
    if abs(r) >= 1
        pval(i) = 0;
    else
        tstat = abs(r) * sqrt(df / (1 - r^2));
        x = df / (df + tstat^2);
        pval(i) = betainc(x, df/2, 0.5);
    end
end

end

% Rank transformation with average ranks for ties


function ranks = rank_vector_avg_ties(x)
N = length(x);
[xs, idx] = sort(x);

naive = zeros(N,1);
naive(idx) = 1:N;

ranks = naive;

i = 1;
while i <= N
    j = i;
    while j < N && xs(j+1) == xs(j)
        j = j + 1;
    end
    
    if j > i
        avgRank = (i + j) / 2;
        ranks(idx(i:j)) = avgRank;
    end
    
    i = j + 1;
end

end


% Significance stars

function starText = significance_stars(p)

if p < 0.001
    starText = '***';
elseif p < 0.01
    starText = '**';
elseif p < 0.05
    starText = '*';
else
    starText = 'ns';
end

end


% Distribution summary printout


function print_distribution_summary(name, Y)

fprintf('Distribution of %s across %d LHS samples:\n', name, length(Y));
fprintf('  mean    = %0.4f\n', mean(Y));
fprintf('  std     = %0.4f\n', std(Y));
fprintf('  min     = %0.4f\n', min(Y));
fprintf('  median  = %0.4f\n', median(Y));
fprintf('  max     = %0.4f\n', max(Y));
fprintf('  P(.<1)  = %0.4f   (fraction of samples with %s < 1)\n', ...
    mean(Y < 1), name);
fprintf('  P(.>1)  = %0.4f\n', mean(Y > 1));

end


% Print PRCC ranking by absolute value

function print_ranking(name, paramNames, PRCC, pval)

[~, ord] = sort(abs(PRCC), 'descend');

fprintf('\n%s sensitivity ranking by |PRCC|:\n', name);

for k = 1:length(ord)
    i = ord(k);
    fprintf('  %d. %-8s  PRCC = %+0.4f  (p = %0.2e %s)\n', ...
        k, paramNames{i}, PRCC(i), pval(i), significance_stars(pval(i)));
end

end


% PRCC bar plot helper

function plot_prcc_bar(PRCC, pval, names)

m = length(PRCC);
hold on;

posMask = PRCC >= 0;
negMask = PRCC < 0;

% Muted journal-style colors for positive and negative PRCC values
positiveColor = [0.30 0.55 0.80];
negativeColor = [0.80 0.45 0.35];

if any(posMask)
    b1 = bar(find(posMask), PRCC(posMask), 0.65);
    set(b1, ...
        'FaceColor', positiveColor, ...
        'EdgeColor', 'k', ...
        'LineWidth', 0.8);
end

if any(negMask)
    b2 = bar(find(negMask), PRCC(negMask), 0.65);
    set(b2, ...
        'FaceColor', negativeColor, ...
        'EdgeColor', 'k', ...
        'LineWidth', 0.8);
end



plot([0 m+1], [0 0], 'k-', 'LineWidth', 1.0);

set(gca, 'XTick', 1:m);
set(gca, 'XTickLabel', names);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 12);

ylabel('PRCC', 'FontName','Times New Roman');

ylim([-1 1]);
xlim([0 m+1]);

grid on;
grid minor;
apply_academic_axes(gca);

hold off;

end


%Padded y-axis limits


function yl = padded_limits(Y)

ymin = min(Y);
ymax = max(Y);

pad = 0.06 * (ymax - ymin);

if pad == 0
    pad = 0.1;
end

yl = [ymin - pad, ymax + pad];

end


%  Academic axis formatting


function apply_academic_axes(ax)

set(ax, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 12, ...
    'LineWidth', 1.0, ...
    'Box', 'on', ...
    'TickDir', 'out', ...
    'TickLength', [0.015 0.015], ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'Layer', 'top');

ax.GridAlpha = 0.20;
ax.MinorGridAlpha = 0.08;

end


function add_super_title(titleText)

if exist('sgtitle', 'file') == 2 || exist('sgtitle', 'builtin') == 5
    sgtitle(titleText, ...
        'FontName', 'Times New Roman', ...
        'FontWeight', 'bold', ...
        'FontSize', 14, ...
        'Interpreter', 'tex');
else
    annotation('textbox', [0 0.95 1 0.04], ...
        'String', titleText, ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'FontName', 'Times New Roman', ...
        'FontWeight', 'bold', ...
        'FontSize', 14);
end

end


function save_publication_figure(figHandle, fileName, resolution)

if exist('exportgraphics', 'file') == 2
    exportgraphics(figHandle, fileName, 'Resolution', resolution);
else
    set(figHandle, 'PaperPositionMode', 'auto');
    print(figHandle, fileName, '-dpng', ['-r', num2str(resolution)]);
end

end