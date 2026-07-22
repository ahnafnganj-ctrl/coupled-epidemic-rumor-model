%% Three Theoretical Regimes 
clear; clc; close all;

p = struct();
p.Lambda   = 0.020000; p.mu       = 0.020000; p.eta      = 0.020000;
p.beta     = 0.700000; p.gamma    = 0.250000; p.psi      = 0.300000;
p.alpha    = 0.500000; p.omega    = 0.050000; p.sigma    = 0.400000;
p.lambda_r = 0.250000; p.phi      = 0.600000; p.delta    = 0.150000;
p.theta    = 0.100000; p.tEnd     = 10000;

p.y0 = [0.89; 0.10; 0.01; 0.00; 0.90; 0.10; 0.00];

fprintf('\n============================================================\n');
fprintf('Objective 2: Three theoretical regimes\n');
fprintf('Standalone version: no Code_Library dependency\n');
fprintf('============================================================\n\n');

C = thesisColors();

caseNames = { ...
    'Case I: disease-free and rumor-free', ...
    'Case II: disease-free, rumor persists', ...
    'Case III: endemic coexistence'};

params = repmat(p,1,3);
params(1).beta     = 0.40;
params(1).lambda_r = 0.10;

params(2).beta     = 0.40;
params(2).lambda_r = 0.25;

params(3) = p;

y0cases = repmat(p.y0,1,3);

logTimeTicksOriginal = [0 10 100 1000 10000];
logTimeTickPositions = logTimeTicksOriginal + 1;
logTimeTickLabels    = {'0','10^1','10^2','10^3','10^4'};
logTimeXLim          = [1 10001];


for j = 1:3
    [t,y] = simulateModel(params(j),y0cases(:,j),params(j).tEnd);
    T = thresholds(params(j));

    fig = figure('Name',caseNames{j},'Color','w','Position',[100 100 850 500]);
    hold on

    tPlot = t + 1;

    h1 = plot(tPlot, y(:,2), 'Color', C(2,:), 'LineWidth', 2.6); % I
    h2 = plot(tPlot, y(:,6), 'Color', C(6,:), 'LineWidth', 2.6); % A
    h3 = plot(tPlot, y(:,3), 'Color', C(3,:), 'LineWidth', 2.3, 'LineStyle','--'); % V

    xlim(logTimeXLim); ylim([0 1]);
    set(gca, 'XScale','log', ...
        'XTick',logTimeTickPositions, ...
        'XTickLabel',logTimeTickLabels, ...
        'FontName','Times New Roman', ...
        'FontSize',11, ...
        'LineWidth',1.10, ...
        'TickDir','out', ...
        'Layer','top', ...
        'XMinorTick','on', ...
        'YMinorTick','on', ...
        'Box','on');
    grid on; ax = gca; ax.GridAlpha = 0.24; ax.MinorGridAlpha = 0.10;
    pbaspect(ax,[1 1 1]);

    title(ax, sprintf('%s\nR_{0E}=%.2f, R_{0R}=%.2f', ...
        caseNames{j}, T.R0E, T.R0R), ...
        'FontName','Times New Roman','FontWeight','bold','FontSize',13,'Interpreter','tex');

    xlabel('Normalized model time, t', 'FontName','Times New Roman','FontSize',12.5);
    ylabel('Normalized population', 'FontName','Times New Roman','FontSize',12.5);

    legend([h1 h2 h3], {'I infected','A active rumor','V vaccinated'}, ...
        'Location','best','Box','off','FontName','Times New Roman','FontSize',11.5);

    hold off

    fprintf('%s: R0E=%.4f, R0R=%.4f, final I=%.8f, final A=%.8f, final V=%.8f.\n', ...
        caseNames{j}, T.R0E, T.R0R, y(end,2), y(end,6), y(end,3));

    saveFigureStandalone(fig,sprintf('MS_objective2_case%d',j));
end

fprintf('\nFigures saved successfully in the folder named Figures.\n\n');



function dydt = coupledModelODE(~,y,p)
    S = max(y(1),0); I = max(y(2),0); V = max(y(3),0); R = max(y(4),0);
    U = max(y(5),0); A = max(y(6),0); C = max(y(7),0);

    dS = p.Lambda - p.beta*S*I - p.psi*(1 - p.alpha*A)*S + p.omega*V - p.mu*S;
    dI = p.beta*S*I + p.sigma*p.beta*V*I - (p.gamma + p.mu)*I;
    dV = p.psi*(1 - p.alpha*A)*S - p.sigma*p.beta*V*I - (p.omega + p.mu)*V;
    dR = p.gamma*I - p.mu*R;
    dU = p.eta - p.lambda_r*U*A + p.theta*C - p.phi*I*U - p.mu*U;
    dA = p.lambda_r*U*A + p.phi*I*U - (p.delta + p.mu)*A;
    dC = p.delta*A - (p.theta + p.mu)*C;

    dydt = [dS; dI; dV; dR; dU; dA; dC];
end

function [t,y] = simulateModel(p,y0,tEnd)
    tEval = unique([0, logspace(-5,log10(tEnd),9000)]);
    opts = odeset('RelTol',1e-10,'AbsTol',1e-12,'NonNegative',1:7,'MaxStep',2.0);
    sol = ode45(@(t,y) coupledModelODE(t,y,p),[0 tEnd],y0,opts);
    t = tEval(:); y = deval(sol,t).';
    y(y < 0 & y > -1e-10) = 0;
end

function T = thresholds(p)
    S0 = p.Lambda*(p.omega + p.mu) / (p.mu*(p.psi + p.omega + p.mu));
    V0 = p.Lambda*p.psi / (p.mu*(p.psi + p.omega + p.mu));
    U0 = p.eta/p.mu;
    T.R0E = p.beta*(S0 + p.sigma*V0)/(p.gamma + p.mu);
    T.R0R = p.lambda_r*U0/(p.delta + p.mu);
    T.beta_c0 = p.mu*(p.gamma + p.mu)*(p.psi + p.omega + p.mu) / ...
                (p.Lambda*((p.omega + p.mu) + p.sigma*p.psi));
    T.lambda_c = p.mu*(p.delta + p.mu)/p.eta;
end

function C = thesisColors()
    C = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250;
         0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330;
         0.6350 0.0780 0.1840; 0.25 0.25 0.25];
end

function saveFigureStandalone(fig,baseName)
    scriptDir = fileparts(mfilename('fullpath'));
    if isempty(scriptDir), scriptDir = pwd; end
    outDir = fullfile(scriptDir,'Figures');
    if ~exist(outDir,'dir'), mkdir(outDir); end
    savefig(fig,fullfile(outDir,[baseName '.fig']));
    try exportgraphics(fig,fullfile(outDir,[baseName '.png']),'Resolution',600); catch, print(fig,fullfile(outDir,[baseName '.png']),'-dpng','-r600'); end
    try exportgraphics(fig,fullfile(outDir,[baseName '.eps']),'ContentType','vector'); catch, print(fig,fullfile(outDir,[baseName '.eps']),'-depsc','-painters'); end
end