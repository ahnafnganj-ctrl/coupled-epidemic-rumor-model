%% Three-panel coupling heatmaps for I_eq, V_eq, and A_eq
clc;
clear;
close all;

% Baseline normalized theoretical parameter set
p.Lambda  = 0.02;
p.mu      = 0.02;
p.eta     = 0.02;
p.beta    = 0.70;
p.gamma   = 0.25;
p.psi     = 0.30;
p.alpha   = 0.50;
p.omega   = 0.05;
p.sigma   = 0.40;
p.lambda_r = 0.25;
p.phi      = 0.60;
p.delta    = 0.15;
p.theta    = 0.10;

% Initial condition
% y = [S I V R U A C]^T
y0 = [0.89; 0.10; 0.01; 0.00; 0.90; 0.10; 0.00];

% Numerical settings
nGrid = 55;
tspan_eq = [0 1000];
odeOptions = odeset('RelTol',1e-7,'AbsTol',1e-9,'NonNegative',1:7);

% Figure style
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultAxesFontSize',12);
set(groot,'defaultTextFontSize',12);
set(groot,'defaultAxesLineWidth',1.10);
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesBox','on');

% Parameter ranges
alpha_values = linspace(0,1,nGrid);
phi_values   = linspace(0,1,nGrid);
sigma_values = linspace(0,1,nGrid);

fprintf('\nComputing three-panel map 1: (alpha, phi) -> I_eq, V_eq, A_eq...\n');
M_I_alpha_phi = computeEquilibriumMatrix(p,y0,tspan_eq,odeOptions,alpha_values,phi_values,'alpha','phi','I');
M_V_alpha_phi = computeEquilibriumMatrix(p,y0,tspan_eq,odeOptions,alpha_values,phi_values,'alpha','phi','V');
M_A_alpha_phi = computeEquilibriumMatrix(p,y0,tspan_eq,odeOptions,alpha_values,phi_values,'alpha','phi','A');

drawThreePanelHeatmap(alpha_values,phi_values, ...
    {M_I_alpha_phi,M_V_alpha_phi,M_A_alpha_phi}, ...
    'Rumor effect on vaccination, \alpha','Infection-induced rumor activation, \phi', ...
    {'I_{eq}: (\alpha,\phi)','V_{eq}: (\alpha,\phi)','A_{eq}: (\alpha,\phi)'}, ...
    'Coupling Heatmaps in the (\alpha,\phi) Plane');

fprintf('Computing three-panel map 2: (sigma, phi) -> I_eq, V_eq, A_eq...\n');
M_I_sigma_phi = computeEquilibriumMatrix(p,y0,tspan_eq,odeOptions,sigma_values,phi_values,'sigma','phi','I');
M_V_sigma_phi = computeEquilibriumMatrix(p,y0,tspan_eq,odeOptions,sigma_values,phi_values,'sigma','phi','V');
M_A_sigma_phi = computeEquilibriumMatrix(p,y0,tspan_eq,odeOptions,sigma_values,phi_values,'sigma','phi','A');

drawThreePanelHeatmap(sigma_values,phi_values, ...
    {M_I_sigma_phi,M_V_sigma_phi,M_A_sigma_phi}, ...
    'Vaccine leakiness, \sigma','Infection-induced rumor activation, \phi', ...
    {'I_{eq}: (\sigma,\phi)','V_{eq}: (\sigma,\phi)','A_{eq}: (\sigma,\phi)'}, ...
    'Coupling Heatmaps in the (\sigma,\phi) Plane');

fprintf('Three-panel coupling heatmaps complete.\n');

function drawThreePanelHeatmap(xValues,yValues,matrixCell,xLabelText,yLabelText,titleCell,mainTitle)
fig = figure('Color','w','Position',[80 120 1450 520]);
colormap(fig,parula(256));
for k = 1:3
    ax = subplot(1,3,k);
    contourf(ax,xValues,yValues,matrixCell{k},36,'LineColor','none');
    colorbar(ax);
    xlabel(ax,xLabelText,'FontName','Times New Roman');
    ylabel(ax,yLabelText,'FontName','Times New Roman');
    title(ax,titleCell{k},'FontName','Times New Roman','FontWeight','bold');
    set(ax,'FontName','Times New Roman','LineWidth',1.10,'TickDir','out','Layer','top','Box','on');
    grid(ax,'off');
    pbaspect(ax,[1 1 1]);
end
sgtitle(mainTitle,'FontName','Times New Roman','FontWeight','bold','FontSize',16);
end

function M = computeEquilibriumMatrix(p,y0,tspan,odeOptions,xValues,yValues,xName,yName,targetName)
M = zeros(length(yValues),length(xValues));
for ix = 1:length(xValues)
    for iy = 1:length(yValues)
        pp = p;
        pp = setParameterValue(pp,xName,xValues(ix));
        pp = setParameterValue(pp,yName,yValues(iy));
        try
            [~,Y] = ode45(@(t,y) coupledModelRHS(t,y,pp),tspan,y0,odeOptions);
            M(iy,ix) = getTargetValue(Y(end,:),targetName);
        catch
            M(iy,ix) = NaN;
        end
    end
end
end

function pp = setParameterValue(pp,name,value)
switch lower(name)
    case 'beta'
        pp.beta = value;
    case 'lambda'
        pp.lambda_r = value;
    case 'delta'
        pp.delta = value;
    case 'phi'
        pp.phi = value;
    case 'sigma'
        pp.sigma = value;
    case 'psi'
        pp.psi = value;
    case 'alpha'
        pp.alpha = value;
    otherwise
        error('Unknown parameter name: %s',name);
end
end

function value = getTargetValue(y,targetName)
switch upper(targetName)
    case 'I'
        value = y(2);
    case 'V'
        value = y(3);
    case 'A'
        value = y(6);
    otherwise
        error('Unknown target name: %s',targetName);
end
end

function dydt = coupledModelRHS(~,y,p)
S = y(1); I = y(2); V = y(3); R = y(4); U = y(5); A = y(6); C = y(7);

dS = p.Lambda - p.beta*S*I - p.psi*(1-p.alpha*A)*S + p.omega*V - p.mu*S;
dI = p.beta*S*I + p.sigma*p.beta*V*I - (p.gamma+p.mu)*I;
dV = p.psi*(1-p.alpha*A)*S - p.sigma*p.beta*V*I - (p.omega+p.mu)*V;
dR = p.gamma*I - p.mu*R;
dU = p.eta - p.lambda_r*U*A - p.phi*I*U + p.theta*C - p.mu*U;
dA = p.lambda_r*U*A + p.phi*I*U - (p.delta+p.mu)*A;
dC = p.delta*A - (p.theta+p.mu)*C;

dydt = [dS; dI; dV; dR; dU; dA; dC];
end

function drawSingleHeatmap(xValues,yValues,M,xLabelText,yLabelText,titleText)
figure('Color','w','Position',[120 80 760 640]);
contourf(xValues,yValues,M,36,'LineColor','none');
colormap(parula(256));
cb = colorbar;
cb.FontName = 'Times New Roman';
xlabel(xLabelText,'FontName','Times New Roman');
ylabel(yLabelText,'FontName','Times New Roman');
title(titleText,'FontName','Times New Roman','FontWeight','bold');
set(gca,'FontName','Times New Roman','LineWidth',1.10,'TickDir','out','Layer','top','Box','on');
grid off;
pbaspect([1 1 1]);
end
