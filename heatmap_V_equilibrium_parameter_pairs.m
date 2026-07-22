%% Heatmaps for vaccinated equilibrium V_eq under selected parameter pairs
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
psi_values   = linspace(0,1,nGrid);
delta_values = linspace(0,1,nGrid);
alpha_values = linspace(0,1,nGrid);
sigma_values = linspace(0,1,nGrid);
phi_values   = linspace(0,1,nGrid);

fprintf('\nComputing V_eq heatmap: (psi, delta)...\n');
M1 = computeEquilibriumMatrix(p,y0,tspan_eq,odeOptions,psi_values,delta_values,'psi','delta','V');
drawSingleHeatmap(psi_values,delta_values,M1, ...
    'Vaccination pressure, \psi','Counter-rumor transition rate, \delta','V_{eq}: (\psi,\delta)');

fprintf('Computing V_eq heatmap: (alpha, psi)...\n');
M2 = computeEquilibriumMatrix(p,y0,tspan_eq,odeOptions,alpha_values,psi_values,'alpha','psi','V');
drawSingleHeatmap(alpha_values,psi_values,M2, ...
    'Rumor effect on vaccination, \alpha','Vaccination pressure, \psi','V_{eq}: (\alpha,\psi)');

fprintf('Computing V_eq heatmap: (sigma, phi)...\n');
M3 = computeEquilibriumMatrix(p,y0,tspan_eq,odeOptions,sigma_values,phi_values,'sigma','phi','V');
drawSingleHeatmap(sigma_values,phi_values,M3, ...
    'Vaccine leakiness, \sigma','Infection-induced rumor activation, \phi','V_{eq}: (\sigma,\phi)');

fprintf('Computing V_eq heatmap: (psi, phi)...\n');
M4 = computeEquilibriumMatrix(p,y0,tspan_eq,odeOptions,psi_values,phi_values,'psi','phi','V');
drawSingleHeatmap(psi_values,phi_values,M4, ...
    'Vaccination pressure, \psi','Infection-induced rumor activation, \phi','V_{eq}: (\psi,\phi)');

fprintf('Computing V_eq heatmap: (alpha, sigma)...\n');
M5 = computeEquilibriumMatrix(p,y0,tspan_eq,odeOptions,alpha_values,sigma_values,'alpha','sigma','V');
drawSingleHeatmap(alpha_values,sigma_values,M5, ...
    'Rumor effect on vaccination, \alpha','Vaccine leakiness, \sigma','V_{eq}: (\alpha,\sigma)');

fprintf('V_eq heatmaps complete.\n');

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
