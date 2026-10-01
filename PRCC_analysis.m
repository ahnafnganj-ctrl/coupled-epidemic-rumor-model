function results = PRCC_analysis
% PRCC_R0E_R0R_ANALYSIS  Invasion-threshold sensitivities and figures.
% MATLAB R2020a; no Statistics and Machine Learning Toolbox required.
% Run: results = PRCC_R0E_R0R_analysis;
% The original R0E/R0R LHS calls and their order are retained.
clc;
clear;
close all;
rng(1);
%% 
% SETTINGS
% Number of Latin Hypercube samples
N = 5000;

% Default: vary the six R0E parameters plus alpha for R_E|1.
% lambda, delta and theta are held at baseline; A1 still changes with mu.
% Set true only to also vary lambda, delta and theta in the R_E|1 analysis.
varyRumourParametersInRE1 = false;

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
colorRE1 = [0.00 0.50 0.35];   % deep green

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

baseline_A1 = (base.theta + base.mu) * ...
    (base.lambda - base.delta - base.mu) / ...
    (base.lambda * (base.delta + base.theta + base.mu));
baseline_q1 = base.psi * (1 - base.alpha * baseline_A1);
baseline_RE1 = base.beta * ...
    (base.omega + base.mu + base.sigma * baseline_q1) / ...
    ((baseline_q1 + base.omega + base.mu) * (base.gamma + base.mu));
assert(baseline_R0R > 1 && baseline_q1 >= 0, ...
    'The baseline must admit the persistent-rumour background E1.');

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
fprintf('Baseline R_E|1 = %.6f\n', baseline_RE1);
fprintf('Baseline A1 = %.6f, q1 = %.6f\n', baseline_A1, baseline_q1);
fprintf('LHS sample size N = %d\n', N);
fprintf('==============================================================\n\n');

fprintf('FIGURE DESCRIPTION\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Figure 1: Combined PRCC bar plots for R_E|1, R0E and R0R.\n');
fprintf('          Positive bars are blue; negative bars are red.\n\n');
fprintf('Figure 2: Scatter plots of sampled epidemic parameters against R0E.\n');
fprintf('          All panels use the same y-axis scale and deep blue markers.\n\n');
fprintf('Figure 3: Scatter plots of sampled rumor parameters against R0R.\n');
fprintf('          All panels use the same y-axis scale and reddish-brown markers.\n');
fprintf('Figure 4: Scatter plots of varied parameters against R_E|1.\n');
fprintf('          All panels use the same y-axis scale and green markers.\n');
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
local_E = local_threshold_indices(base, plainNamesE, 1);

fprintf('\n==============================================================\n');
fprintf('PRCC ANALYSIS FOR EPIDEMIC REPRODUCTION NUMBER R0E\n');
fprintf('==============================================================\n');
fprintf('R0E = beta * ((omega+mu) + sigma*psi) / ((omega+psi+mu)*(gamma+mu))\n\n');

print_distribution_summary('R0E', R0E);

fprintf('\n%-10s   %-12s   %-10s   %-12s   %-6s\n', ...
    'Parameter','Local index','PRCC','p-value','Sig.');
fprintf('%s\n', repmat('-',1,70));
for i = 1:mE
    starText = significance_stars(pval_E(i));
    fprintf('%-10s   %+0.6f      %+0.6f    %0.3e     %-6s\n', ...
        plainNamesE{i}, local_E(i), PRCC_E(i), pval_E(i), starText);
end

print_ranking('R0E', plainNamesE, PRCC_E, pval_E);

%Figure 1: Combined PRCC bar plot in the style of main_thesis Figure 11
fig1 = figure('Color','w','Position',[60 100 1500 500], ...
    'Renderer','painters');

co = [31 119 180; 214 39 40; 44 160 44; ...
      148 103 189; 255 127 14]/255;

labels_E = arrayfun(@(x) sprintf('%.2f',x), ...
    PRCC_E,'UniformOutput',false);

axE = axes('Parent',fig1,'Position',[.46 .16 .29 .73]);
bars_fig11(axE,PRCC_E, ...
    {'\bf\beta','\bf\omega','\bf\mu','\bf\psi','\bf\sigma','\bf\gamma'}, ...
    '\bf(b) \rm{\itR}_{0E}',co,labels_E);

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
local_R = local_threshold_indices(base, plainNamesR, 2);

fprintf('\n==============================================================\n');
fprintf('PRCC ANALYSIS FOR RUMOR REPRODUCTION NUMBER R0R\n');
fprintf('==============================================================\n');
fprintf('R0R = lambda / (delta + mu)\n\n');

print_distribution_summary('R0R', R0R);

fprintf('\n%-10s   %-12s   %-10s   %-12s   %-6s\n', ...
    'Parameter','Local index','PRCC','p-value','Sig.');
fprintf('%s\n', repmat('-',1,70));
for i = 1:mR
    starText = significance_stars(pval_R(i));
    fprintf('%-10s   %+0.6f      %+0.6f    %0.3e     %-6s\n', ...
        plainNamesR{i}, local_R(i), PRCC_R(i), pval_R(i), starText);
end

print_ranking('R0R', plainNamesR, PRCC_R, pval_R);

% Complete Figure 1 with the rumour PRCC panel
labels_R = arrayfun(@(x) sprintf('%.2f',x), ...
    PRCC_R,'UniformOutput',false);

axR = axes('Parent',fig1,'Position',[.82 .16 .16 .73]);
bars_fig11(axR,PRCC_R, ...
    {'\bf\lambda','\bf\delta','\bf\mu'}, ...
    '\bf(c) \rm{\itR}_{0R}',co,labels_R);

add_super_title('Scatter Plots of Sampled Rumor Parameters versus R_{0R}');

if saveFigures
    save_publication_figure(fig3, 'Figure_PRCC_R0R_scatter.png', exportResolution);
end

%% PRCC FOR DISEASE INVASION AT THE PERSISTENT-RUMOUR BACKGROUND
% Sample this block AFTER the original R0E and R0R analyses. This preserves
% their random-number draws and their PRCC results for the same MATLAB RNG.
paramNames1 = [paramNamesE, {'\alpha'}];
plainNames1 = [plainNamesE, {'alpha'}];
if varyRumourParametersInRE1
    paramNames1 = [paramNames1, {'\lambda','\delta','\theta'}];
    plainNames1 = [plainNames1, {'lambda','delta','theta'}];
end
m1 = numel(plainNames1);
base1 = zeros(1,m1);
for i = 1:m1
    base1(i) = base.(plainNames1{i});
end
X1all = manual_lhs(N, 0.5*base1, 1.5*base1);

beta1  = X1all(:,1);
omega1 = X1all(:,2);
mu1    = X1all(:,3);
psi1   = X1all(:,4);
sigma1 = X1all(:,5);
gamma1 = X1all(:,6);
alpha1 = X1all(:,7);
if varyRumourParametersInRE1
    lambda1 = X1all(:,8);
    delta1  = X1all(:,9);
    theta1  = X1all(:,10);
else
    lambda1 = base.lambda * ones(N,1);
    delta1  = base.delta * ones(N,1);
    theta1  = base.theta * ones(N,1);
end

% Unit layer totals: Lambda = eta = mu in every sample.
% A1 is recomputed, including its dependence on sampled mu.
R0R_1 = lambda1 ./ (delta1 + mu1);
A1all = (theta1 + mu1) .* (lambda1 - delta1 - mu1) ./ ...
    (lambda1 .* (delta1 + theta1 + mu1));
q1all = psi1 .* (1 - alpha1 .* A1all);
RE1all = beta1 .* (omega1 + mu1 + sigma1 .* q1all) ./ ...
    ((q1all + omega1 + mu1) .* (gamma1 + mu1));

% R_E|1 is used only where the persistent-rumour background exists.
% All 5000 samples pass under the default seven-parameter settings.
valid1 = R0R_1 > 1 & A1all > 0 & A1all <= 1 & q1all >= 0 ...
    & isfinite(RE1all);
N1 = sum(valid1);
assert(N1 > m1 + 1, 'Too few admissible samples for R_E|1 PRCC.');
X1 = X1all(valid1,:);
RE1 = RE1all(valid1);
[PRCC_1, pval_1] = compute_prcc(X1, RE1);
local_1 = local_threshold_indices(base, plainNames1, 3);

fprintf('\n==============================================================\n');
fprintf('PRCC ANALYSIS FOR DISEASE INVASION NUMBER R_E|1\n');
fprintf('==============================================================\n');
fprintf('A1 = (theta+mu)*(lambda-delta-mu)/(lambda*(delta+theta+mu))\n');
fprintf('q1 = psi*(1-alpha*A1)\n');
fprintf('R_E|1 = beta*(omega+mu+sigma*q1)/((q1+omega+mu)*(gamma+mu))\n');
if varyRumourParametersInRE1
    fprintf('Ten parameters vary, including lambda, delta and theta.\n');
else
    fprintf('Seven parameters vary: beta, omega, mu, psi, sigma, gamma, alpha.\n');
    fprintf('lambda=%.4f, delta=%.4f, theta=%.4f are held at baseline.\n', ...
        base.lambda, base.delta, base.theta);
end
fprintf('A1 is recomputed for every sample; Lambda = eta = mu.\n');
fprintf('Admissible R_E|1 samples: %d of %d (excluded: %d).\n\n', ...
    N1, N, N-N1);
print_distribution_summary('R_E|1', RE1);
fprintf('\n%-10s   %-12s   %-10s   %-12s   %-6s\n', ...
    'Parameter','Local index','PRCC','p-value','Sig.');
fprintf('%s\n', repmat('-',1,70));
for i = 1:m1
    fprintf('%-10s   %+0.6f      %+0.6f    %0.3e     %-6s\n', ...
        plainNames1{i}, local_1(i), PRCC_1(i), pval_1(i), ...
        significance_stars(pval_1(i)));
end
print_ranking('R_E|1', plainNames1, PRCC_1, pval_1);

% Complete Figure 1 with the additional R_E|1 panel.
labels_1 = arrayfun(@(x) sprintf('%.2f',x), PRCC_1, ...
    'UniformOutput',false);
boldNames1 = cellfun(@(s) ['\bf' s], paramNames1, ...
    'UniformOutput',false);
ax1 = axes('Parent',fig1,'Position',[.055 .16 .34 .73]);
bars_fig11(ax1, PRCC_1, boldNames1, ...
    '\bf(a) \rm{\itR}_{E|1}', co, labels_1);
if saveFigures
    save_publication_figure(fig1, ...
        'Figure_PRCC_RE1_R0E_R0R_combined.png', exportResolution);
end

% Consolidated table, in the same output order as the combined figure.
outputNames = [repmat({'R_E|1'},m1,1); repmat({'R0E'},mE,1); ...
    repmat({'R0R'},mR,1)];
parameterNames = [plainNames1(:); plainNamesE(:); plainNamesR(:)];
summaryTable = table(outputNames, parameterNames, ...
    [local_1; local_E; local_R], [PRCC_1; PRCC_E; PRCC_R], ...
    [pval_1; pval_E; pval_R], ...
    'VariableNames',{'Output','Parameter','LocalIndex','PRCC','PValue'});
fprintf('\nCOMBINED LOCAL-INDEX AND PRCC TABLE\n');
fprintf('%-8s %-10s %12s %10s %12s\n', ...
    'Output','Parameter','Local index','PRCC','p-value');
for i = 1:height(summaryTable)
    fprintf('%-8s %-10s %+12.4f %+10.2f %12.3e\n', ...
        summaryTable.Output{i}, summaryTable.Parameter{i}, ...
        summaryTable.LocalIndex(i), summaryTable.PRCC(i), ...
        summaryTable.PValue(i));
end

% Return numerical results without saving data or figures automatically.
results.baselineParameters = base;
results.baseline = struct('R0E',baseline_R0E,'R0R',baseline_R0R, ...
    'RE1',baseline_RE1,'A1',baseline_A1,'q1',baseline_q1);
results.settings = struct('N',N,'seed',1,'rangeFactors',[0.5 1.5], ...
    'varyRumourParametersInRE1',varyRumourParametersInRE1, ...
    'saveFigures',saveFigures);
results.settings.randomGenerator = rng;
results.summaryTable = summaryTable;
results.R0E = struct('parameters',{plainNamesE},'samples',XE, ...
    'values',R0E,'localIndices',local_E,'PRCC',PRCC_E,'pValues',pval_E);
results.R0R = struct('parameters',{plainNamesR},'samples',XR, ...
    'values',R0R,'localIndices',local_R,'PRCC',PRCC_R,'pValues',pval_R);
results.RE1 = struct('parameters',{plainNames1},'samples',X1, ...
    'values',RE1,'localIndices',local_1,'PRCC',PRCC_1,'pValues',pval_1, ...
    'candidateSamples',X1all,'validMask',valid1,'nRetained',N1, ...
    'A1',A1all(valid1),'q1',q1all(valid1));
results.figures = fig1;
drawnow;

fprintf('\n==============================================================\n');
fprintf('ANALYSIS COMPLETE\n');
if saveFigures
    fprintf('Generated figures:\n');
    fprintf('  Figure_PRCC_RE1_R0E_R0R_combined.png\n');


else
    fprintf('Figures were generated but not saved. Set saveFigures = true to export them.\n');
end
fprintf('==============================================================\n\n');

end

% Baseline local elasticities: complex-step derivatives of the formulas.
% Output indices: 1 = R0E, 2 = R0R, 3 = R_E|1.
function indices = local_threshold_indices(base, names, outputIndex)
Y0 = threshold_values(base);
h = 1e-20;
indices = zeros(numel(names),1);
for j = 1:numel(names)
    p = base;
    name = names{j};
    p.(name) = base.(name) + 1i*h;
    if strcmp(name,'mu')
        p.Lambda = p.mu;
        p.eta = p.mu;
    end
    Y = threshold_values(p);
    indices(j) = base.(name) * imag(Y(outputIndex)) / ...
        (h * Y0(outputIndex));
end
end

function Y = threshold_values(p)
A1 = (p.theta+p.mu)*(p.lambda-p.delta-p.mu) / ...
    (p.lambda*(p.delta+p.theta+p.mu));
q1 = p.psi*(1-p.alpha*A1);
R0E = p.beta*(p.omega+p.mu+p.sigma*p.psi) / ...
    ((p.psi+p.omega+p.mu)*(p.gamma+p.mu));
R0R = p.lambda/(p.delta+p.mu);
RE1 = p.beta*(p.omega+p.mu+p.sigma*q1) / ...
    ((q1+p.omega+p.mu)*(p.gamma+p.mu));
Y = [R0E, R0R, RE1];
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


% Combined PRCC plotting helper: identical styling to main_thesis Fig. 11
function bars_fig11(ax,vals,names,titleText,co_palette,manual_labels)

hold(ax,'on');
hb = bar(ax,1:numel(vals),vals,0.65, ...
    'FaceColor','flat','EdgeColor','none');

colors = repmat(co_palette(1,:),numel(vals),1);
colors(vals<0,:) = repmat(co_palette(2,:),sum(vals<0),1);
hb.CData = colors;
yline(ax,0,'Color',[.4 .4 .4],'LineWidth',1.2);

for k = 1:numel(vals)
    if vals(k)>=0
        offset=.06;
        va='bottom';
    else
        offset=-.06;
        va='top';
    end
    text(ax,k,vals(k)+offset,manual_labels{k}, ...
        'HorizontalAlignment','center','VerticalAlignment',va, ...
        'FontName','Arial','FontSize',11.5);
end

set(ax,'XTick',1:numel(vals),'XTickLabel',names, ...
    'TickLabelInterpreter','tex');
ylim(ax,[-1.15 1.15]);
xlim(ax,[.4 numel(vals)+.6]);
ylabel(ax,'PRCC','FontName','Arial','FontSize',14);

grid(ax,'on');
set(ax,'FontName','Arial','FontSize',12,'FontWeight','normal', ...
    'TickLabelInterpreter','tex','Box','off','TickDir','out', ...
    'GridColor',[.8 .8 .8],'GridAlpha',.5,'LineWidth',1.0);

title(ax,titleText,'FontName','Arial','FontSize',15, ...
    'FontWeight','normal','Interpreter','tex');
hold(ax,'off');

end
