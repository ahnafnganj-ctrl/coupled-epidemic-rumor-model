%% Scenario Plots for I(t), A(t), and V(t)

clear; clc; close all;

% BASELINE PARAMETER VALUES
% ============================================================

p = struct();

% Demographic / recruitment parameters
p.Lambda = 0.020000;
p.mu     = 0.020000;
p.eta    = 0.020000;

% Epidemic-layer parameters
p.beta   = 0.700000;
p.gamma  = 0.250000;
p.psi    = 0.300000;
p.alpha  = 0.500000;
p.omega  = 0.050000;
p.sigma  = 0.400000;

% Rumor-layer parameters
p.lambda_r = 0.250000;
p.phi      = 0.600000;
p.delta    = 0.150000;
p.theta    = 0.100000;

% Simulation horizon
p.tEnd = 10000;

% Initial condition
% State order: y = [S I V R U A C]
p.y0 = [0.89; 0.10; 0.01; 0.00; 0.90; 0.10; 0.00];


% SCENARIO DEFINITIONS
% ============================================================

scenarioNames = { ...
    'Baseline', ...
    'Low rumor transmission', ...
    'High rumor transmission', ...
    'High vaccination', ...
    'High leakiness'};

sc = repmat(p,1,5);

% Scenario 1: Baseline
sc(1) = p;

% Scenario 2: Low rumor transmission
sc(2).lambda_r = 0.05;

% Scenario 3: High rumor transmission
sc(3).lambda_r = 0.90;

% Scenario 4: High vaccination pressure
sc(4).psi = 0.70;

% Scenario 5: High vaccine leakiness
sc(5).sigma = 0.70;

% PRINT BASELINE AND SCENARIO INFORMATION
% ============================================================

fprintf('\n============================================================\n');
fprintf('Standalone scenario plots for I(t), A(t), and V(t)\n');
fprintf('No Code_Library dependency\n');
fprintf('============================================================\n\n');

fprintf('Baseline parameter values:\n');
fprintf('Lambda   = %.4f\n',p.Lambda);
fprintf('mu       = %.4f\n',p.mu);
fprintf('eta      = %.4f\n',p.eta);
fprintf('beta     = %.4f\n',p.beta);
fprintf('gamma    = %.4f\n',p.gamma);
fprintf('psi      = %.4f\n',p.psi);
fprintf('alpha    = %.4f\n',p.alpha);
fprintf('omega    = %.4f\n',p.omega);
fprintf('sigma    = %.4f\n',p.sigma);
fprintf('lambda_r = %.4f\n',p.lambda_r);
fprintf('phi      = %.4f\n',p.phi);
fprintf('delta    = %.4f\n',p.delta);
fprintf('theta    = %.4f\n',p.theta);
fprintf('tEnd     = %.0f\n\n',p.tEnd);

fprintf('Scenario definitions:\n');
for j = 1:5
    fprintf('%d. %-25s  lambda_r = %.3f, psi = %.3f, sigma = %.3f\n', ...
        j, scenarioNames{j}, sc(j).lambda_r, sc(j).psi, sc(j).sigma);
end
fprintf('\n');


% SIMULATE ALL SCENARIOS ONCE
% ============================================================

fprintf('Simulating all scenarios...\n');

Tcell = cell(1,5);
Ycell = cell(1,5);

for j = 1:5
    [t,y] = simulateModel(sc(j),p.y0,p.tEnd);
    Tcell{j} = t;
    Ycell{j} = y;

    fprintf('%-25s final I = %.6f, final A = %.6f, final V = %.6f\n', ...
        scenarioNames{j}, y(end,2), y(end,6), y(end,3));
end

fprintf('\nSimulation complete.\n\n');


% FIGURE SETTINGS
% ============================================================

set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',12);
set(groot,'defaultTextFontSize',12);
set(groot,'defaultAxesLineWidth',1.05);
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesBox','on');

C = thesisColors();

lineWidth = 2.30;

% Shifted log-time axis.
% MATLAB cannot display t = 0 on a true logarithmic axis, so x = t + 1 is used.
logTimeTicksOriginal = [0 10 100 1000 10000];
logTimeTicksOriginal = logTimeTicksOriginal(logTimeTicksOriginal <= p.tEnd);
logTimeTickPositions = logTimeTicksOriginal + 1;

logTimeTickLabels = cell(size(logTimeTicksOriginal));
for k = 1:length(logTimeTicksOriginal)
    if logTimeTicksOriginal(k) == 0
        logTimeTickLabels{k} = '0';
    else
        logTimeTickLabels{k} = sprintf('10^{%d}',round(log10(logTimeTicksOriginal(k))));
    end
end

logTimeXLim = [1 p.tEnd + 1];

% FIGURES: I(t), A(t), V(t)
% ============================================================

varIndex = [2 6 3];

varTitles = { ...
    'Infected fraction I(t)', ...
    'Active rumor fraction A(t)', ...
    'Vaccinated fraction V(t)'};

varYLabels = { ...
    'Infected fraction, I(t)', ...
    'Active-rumor fraction, A(t)', ...
    'Vaccinated fraction, V(t)'};

saveNames = { ...
    'Scenario_Infected_I', ...
    'Scenario_ActiveRumor_A', ...
    'Scenario_Vaccinated_V'};

for v = 1:3

    fig = figure( ...
        'Name',varTitles{v}, ...
        'Color','w', ...
        'Position',[120 100 980 620]);

    ax = axes('Parent',fig);
    hold(ax,'on');

    for j = 1:5

        t = Tcell{j};
        y = Ycell{j};

        tPlot = t + 1;

        plot(ax,tPlot,y(:,varIndex(v)), ...
            'Color',C(j,:), ...
            'LineWidth',lineWidth);

    end

    xlim(ax,logTimeXLim);
    ylim(ax,[0 1]);

    set(ax, ...
        'XScale','log', ...
        'XTick',logTimeTickPositions, ...
        'XTickLabel',logTimeTickLabels, ...
        'FontName','Times New Roman', ...
        'FontSize',12.5, ...
        'LineWidth',1.10, ...
        'TickDir','out', ...
        'Layer','top', ...
        'XMinorTick','on', ...
        'YMinorTick','on', ...
        'Box','on');

    grid(ax,'on');
    ax.GridAlpha = 0.24;
    ax.MinorGridAlpha = 0.10;

    xlabel(ax,'Normalized model time, t', ...
        'FontName','Times New Roman', ...
        'FontSize',13);

    ylabel(ax,varYLabels{v}, ...
        'FontName','Times New Roman', ...
        'FontSize',13);

    title(ax,varTitles{v}, ...
        'FontName','Times New Roman', ...
        'FontWeight','bold', ...
        'FontSize',15);

    legend(ax,scenarioNames, ...
        'Location','eastoutside', ...
        'FontName','Times New Roman', ...
        'FontSize',11, ...
        'Box','on');

    annotation(fig,'textbox',[0.16 0.015 0.68 0.045], ...
        'String','Scenario comparison under low/high rumor transmission, high vaccination pressure, and high vaccine leakiness.', ...
        'EdgeColor','none', ...
        'HorizontalAlignment','center', ...
        'FontName','Times New Roman', ...
        'FontSize',11);

    hold(ax,'off');

    saveFigureStandalone(fig,saveNames{v});

end

fprintf('Three figures saved successfully in the folder named Figures.\n\n');


% LOCAL FUNCTIONS
% =====================================================================

function dydt = coupledModelODE(~,y,p)

    S = y(1);
    I = y(2);
    V = y(3);
    R = y(4);
    U = y(5);
    A = y(6);
    C = y(7);

    dS = p.Lambda ...
        - p.beta*S*I ...
        - p.psi*(1 - p.alpha*A)*S ...
        + p.omega*V ...
        - p.mu*S;

    dI = p.beta*S*I ...
        + p.sigma*p.beta*V*I ...
        - (p.gamma + p.mu)*I;

    dV = p.psi*(1 - p.alpha*A)*S ...
        - p.sigma*p.beta*V*I ...
        - (p.omega + p.mu)*V;

    dR = p.gamma*I ...
        - p.mu*R;

    dU = p.eta ...
        - p.lambda_r*U*A ...
        - p.phi*I*U ...
        + p.theta*C ...
        - p.mu*U;

    dA = p.lambda_r*U*A ...
        + p.phi*I*U ...
        - (p.delta + p.mu)*A;

    dC = p.delta*A ...
        - (p.theta + p.mu)*C;

    dydt = [dS; dI; dV; dR; dU; dA; dC];

end

function [t,y] = simulateModel(p,y0,tEnd)

    % Dense log-focused output grid for smooth early and long-term trajectories.
    tEval = unique([ ...
        0, ...
        logspace(-5,log10(tEnd),6000)]);

    % Guard against tiny floating-point endpoint overshoot.
    tEval = tEval(tEval <= tEnd);
    if tEval(end) < tEnd
        tEval = [tEval, tEnd];
    end

    opts = odeset( ...
        'RelTol',1e-10, ...
        'AbsTol',1e-12, ...
        'NonNegative',1:7, ...
        'MaxStep',2.0);

    [t,y] = ode45(@(t,y) coupledModelODE(t,y,p),tEval,y0,opts);

    % Clean very small numerical round-off values.
    y(y < 0 & y > -1e-10) = 0;

end

function C = thesisColors()

    C = [ ...
        0.0000 0.4470 0.7410;   % blue
        0.8500 0.3250 0.0980;   % orange
        0.9290 0.6940 0.1250;   % yellow
        0.4940 0.1840 0.5560;   % purple
        0.4660 0.6740 0.1880;   % green
        0.3010 0.7450 0.9330;   % cyan
        0.6350 0.0780 0.1840;   % dark red
        0.2500 0.2500 0.2500];  % gray

end

function saveFigureStandalone(fig,baseName)

    scriptDir = fileparts(mfilename('fullpath'));

    if isempty(scriptDir)
        scriptDir = pwd;
    end

    outDir = fullfile(scriptDir,'Figures');

    if ~exist(outDir,'dir')
        mkdir(outDir);
    end

    savefig(fig,fullfile(outDir,[baseName '.fig']));

    try
        exportgraphics(fig,fullfile(outDir,[baseName '.png']),'Resolution',600);
    catch
        print(fig,fullfile(outDir,[baseName '.png']),'-dpng','-r600');
    end

    try
        exportgraphics(fig,fullfile(outDir,[baseName '.eps']),'ContentType','vector');
    catch
        print(fig,fullfile(outDir,[baseName '.eps']),'-depsc','-painters');
    end

end