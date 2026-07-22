%% Chapter 5 baseline stability – baseline trajectories only (Figures 1A & 1B)
%  Self‑contained version: all helper functions are included as local
%  functions at the end of this script.  No external .m files are required.

clc;
clear;
close all;

% ============================================================
% OUTPUT SETTINGS
% ============================================================

saveFigures = true;
exportResolution = 600;

outDir = fullfile(pwd, 'chapter5_baseline_stability_output');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% ============================================================
% GLOBAL FIGURE STYLE: RESEARCH-PAPER STANDARD
% ============================================================

set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 11);
set(groot, 'defaultTextFontSize', 11);
set(groot, 'defaultAxesLineWidth', 1.05);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'on');
set(groot, 'defaultLineLineWidth', 1.75);

lineColors = [
    0.0000 0.2000 0.9000   % S: saturated blue
    0.9000 0.0500 0.0500   % I: saturated red
    0.0000 0.6500 0.1500   % V: saturated green
    1.0000 0.7500 0.0000   % R: strong yellow/golden yellow
    0.4660 0.6740 0.1880
    0.3010 0.7450 0.9330
    0.6350 0.0780 0.1840
];

% ============================================================
% LATEST NORMALIZED THEORETICAL BASELINE
% ============================================================

p.Lambda = 0.02;
p.mu     = 0.02;
p.eta    = 0.02;

p.beta   = 0.70;
p.gamma  = 0.25;
p.psi    = 0.30;
p.alpha  = 0.50;
p.omega  = 0.05;
p.sigma  = 0.40;

p.lambda = 0.25;
p.phi    = 0.60;
p.delta  = 0.15;
p.theta  = 0.10;

% ============================================================
% INITIAL CONDITION USED IN THE THESIS
% ============================================================

S0 = 0.89;
I0 = 0.10;
V0 = 0.01;
R0 = 0.00;

U0 = 0.90;
A0 = 0.10;
C0 = 0.00;

y0 = [S0; I0; V0; R0; U0; A0; C0];

% ============================================================
% NUMERICAL ENDEMIC EQUILIBRIUM FROM THE LATEST CHAPTER 4
% ============================================================

E_star = [
    0.15899202
    0.02031128
    0.56680567
    0.25389103
    0.54736463
    0.20117127
    0.25146409
];

% ============================================================
% EQUILIBRIA, REPRODUCTION NUMBERS, AND THRESHOLDS
% ============================================================

EQ = computeEq(p);

R0E      = EQ.R0E;
R0R      = EQ.R0R;
lambda_c = EQ.lambda_c;
beta_c0  = EQ.beta_c0;
beta_cr  = EQ.beta_cr;

fprintf('\n==============================================================\n');
fprintf('CHAPTER 5: BASELINE STABILITY – TRAJECTORY OUTPUT\n');
fprintf('==============================================================\n');
fprintf('All rates are per normalized model time unit.\n');
fprintf('The model is normalized so that S+I+V+R -> 1 and U+A+C -> 1.\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Latest baseline parameter set:\n');
fprintf('Lambda  = %.6f\n', p.Lambda);
fprintf('mu      = %.6f\n', p.mu);
fprintf('eta     = %.6f\n', p.eta);
fprintf('beta    = %.6f\n', p.beta);
fprintf('gamma   = %.6f\n', p.gamma);
fprintf('psi     = %.6f\n', p.psi);
fprintf('alpha   = %.6f\n', p.alpha);
fprintf('omega   = %.6f\n', p.omega);
fprintf('sigma   = %.6f\n', p.sigma);
fprintf('lambda  = %.6f\n', p.lambda);
fprintf('phi     = %.6f\n', p.phi);
fprintf('delta   = %.6f\n', p.delta);
fprintf('theta   = %.6f\n', p.theta);
fprintf('--------------------------------------------------------------\n');
fprintf('Baseline reproduction numbers:\n');
fprintf('R0E = %.8f\n', R0E);
fprintf('R0R = %.8f\n', R0R);
fprintf('--------------------------------------------------------------\n');
fprintf('Analytical threshold values:\n');
fprintf('lambda_c   = %.8f\n', lambda_c);
fprintf('beta_c^(0) = %.8f\n', beta_c0);
fprintf('beta_c^(r) = %.8f\n', beta_cr);
fprintf('--------------------------------------------------------------\n');
fprintf('Rumor-present disease-free equilibrium E_r:\n');
fprintf('S_r = %.8f, I_r = %.8f, V_r = %.8f, R_r = %.8f\n', ...
    EQ.Er(1), EQ.Er(2), EQ.Er(3), EQ.Er(4));
fprintf('U_r = %.8f, A_r = %.8f, C_r = %.8f\n', ...
    EQ.Er(5), EQ.Er(6), EQ.Er(7));
fprintf('psi_r = %.8f\n', EQ.psi_r);
fprintf('==============================================================\n\n');

% ============================================================
% BASELINE TRAJECTORY SIMULATION
% ============================================================

tEnd = 10000;

tEval = unique([0, logspace(-3, log10(tEnd), 1400)]);

odeOptions = odeset( ...
    'RelTol', 1e-9, ...
    'AbsTol', 1e-11, ...
    'NonNegative', 1:7);

fprintf('Solving baseline trajectory up to t = %.0f...\n', tEnd);

[T, Y] = ode45(@(t,y) rhsModel(t, y, p), tEval, y0, odeOptions);

S = Y(:,1);
I = Y(:,2);
V = Y(:,3);
R = Y(:,4);
U = Y(:,5);
A = Y(:,6);
C = Y(:,7);

idxLog = T >= 0;
Tplot = T(idxLog) + 1;

logTimeTicksOriginal = [0 10 100 1000 10000];
logTimeTickPositions = logTimeTicksOriginal + 1;
logTimeTickLabels = {'0','10^1','10^2','10^3','10^4'};
logTimeXLim = [1 10001];

N_disease = S + I + V + R;
N_rumor   = U + A + C;

massErrorDisease = abs(N_disease - 1);
massErrorRumor   = abs(N_rumor - 1);

yEnd = Y(end,:)';
resEnd = rhsModel(0, yEnd, p);

residualInfNorm = max(abs(resEnd));
distanceFinal = max(abs(yEnd - E_star));

lastWindowIndex = T >= 0.90*tEnd;
lastWindowDrift = max(max(abs(Y(lastWindowIndex,:) - yEnd')));

vaccinationFactor = 1 - p.alpha*A;
minVaccinationFactor = min(vaccinationFactor);

maxMassErrorDisease = max(massErrorDisease);
maxMassErrorRumor   = max(massErrorRumor);

[peakI, peakIIndex] = max(I);
tPeakI = T(peakIIndex);

[peakA, peakAIndex] = max(A);
tPeakA = T(peakAIndex);

fprintf('\n==============================================================\n');
fprintf('BASELINE TRAJECTORY NUMERICAL DIAGNOSTICS\n');
fprintf('==============================================================\n');
fprintf('Final state y(tEnd):\n');
fprintf('S(tEnd) = %.10f\n', yEnd(1));
fprintf('I(tEnd) = %.10f\n', yEnd(2));
fprintf('V(tEnd) = %.10f\n', yEnd(3));
fprintf('R(tEnd) = %.10f\n', yEnd(4));
fprintf('U(tEnd) = %.10f\n', yEnd(5));
fprintf('A(tEnd) = %.10f\n', yEnd(6));
fprintf('C(tEnd) = %.10f\n', yEnd(7));
fprintf('--------------------------------------------------------------\n');
fprintf('Endemic equilibrium E* used for Chapter 4 local stability:\n');
fprintf('S* = %.10f\n', E_star(1));
fprintf('I* = %.10f\n', E_star(2));
fprintf('V* = %.10f\n', E_star(3));
fprintf('R* = %.10f\n', E_star(4));
fprintf('U* = %.10f\n', E_star(5));
fprintf('A* = %.10f\n', E_star(6));
fprintf('C* = %.10f\n', E_star(7));
fprintf('--------------------------------------------------------------\n');
fprintf('Residual and convergence checks:\n');
fprintf('||f(y(tEnd))||_inf       = %.4e\n', residualInfNorm);
fprintf('||y(tEnd)-E*||_inf       = %.4e\n', distanceFinal);
fprintf('Maximum final-window drift = %.4e\n', lastWindowDrift);
fprintf('--------------------------------------------------------------\n');
fprintf('Mass conservation checks:\n');
fprintf('max |S+I+V+R-1|          = %.4e\n', maxMassErrorDisease);
fprintf('max |U+A+C-1|            = %.4e\n', maxMassErrorRumor);
fprintf('--------------------------------------------------------------\n');
fprintf('Admissibility check:\n');
fprintf('min(1-alpha*A(t))        = %.8f\n', minVaccinationFactor);
fprintf('--------------------------------------------------------------\n');
fprintf('Transient peak values:\n');
fprintf('Peak infected fraction    = %.8f at t = %.4f\n', peakI, tPeakI);
fprintf('Peak active-rumor fraction= %.8f at t = %.4f\n', peakA, tPeakA);
fprintf('==============================================================\n\n');

% ============================================================
% JACOBIAN AND LOCAL STABILITY CHECK AT E*
% ============================================================

Jstar = jacobianAt(p, E_star);
eigStar = eig(Jstar);
maxRealEig = max(real(eigStar));

fprintf('\n==============================================================\n');
fprintf('LOCAL STABILITY CHECK AT THE ENDEMIC EQUILIBRIUM\n');
fprintf('==============================================================\n');
fprintf('Eigenvalues of the full 7-dimensional Jacobian J(E*):\n');
for k = 1:length(eigStar)
    fprintf('xi_%d = %+.10f %+.10fi\n', k, real(eigStar(k)), imag(eigStar(k)));
end
fprintf('--------------------------------------------------------------\n');
fprintf('Maximum real part = %.10e\n', maxRealEig);
fprintf('Local stability conclusion: %s\n', bool2str(maxRealEig < 0));
fprintf('==============================================================\n\n');

% ============================================================
% FIGURE 1A: DISEASE-LAYER BASELINE TRAJECTORIES
% ============================================================

fig1A = figure('Color','w','Position',[80 60 820 620]);

hold on;
plot(Tplot, S(idxLog), 'Color', lineColors(1,:), 'LineWidth', 2.0);
plot(Tplot, I(idxLog), 'Color', lineColors(2,:), 'LineWidth', 2.0);
plot(Tplot, V(idxLog), 'Color', lineColors(3,:), 'LineWidth', 2.0);
plot(Tplot, R(idxLog), 'Color', lineColors(4,:), 'LineWidth', 2.0);
yline(E_star(1),'--','Color',lineColors(1,:),'LineWidth',1.0);
yline(E_star(2),'--','Color',lineColors(2,:),'LineWidth',1.0);
yline(E_star(3),'--','Color',lineColors(3,:),'LineWidth',1.0);
yline(E_star(4),'--','Color',lineColors(4,:),'LineWidth',1.0);
xlabel('Normalized model time, t');
ylabel('Population fraction');
title('Disease-layer baseline trajectories');
legend({'S(t)','I(t)','V(t)','R(t)'}, 'Location','eastoutside', 'Box','off');
ylim([0 1]);
xlim(logTimeXLim);
grid on;
axis square;
box on;
set(gca, ...
    'FontName','Times New Roman', ...
    'LineWidth',1.05, ...
    'TickDir','out', ...
    'Layer','top', ...
    'XScale','log', ...
    'XTick',logTimeTickPositions, ...
    'XTickLabel',logTimeTickLabels, ...
    'XMinorTick','on', ...
    'YMinorTick','on');

if saveFigures
    if exist('exportgraphics','file') == 2
        exportgraphics(fig1A, fullfile(outDir,'Figure_1A_Disease_Layer_Baseline_Trajectories.png'), ...
            'Resolution', exportResolution);
    else
        set(fig1A,'PaperPositionMode','auto');
        print(fig1A, fullfile(outDir,'Figure_1A_Disease_Layer_Baseline_Trajectories.png'), ...
            '-dpng', ['-r',num2str(exportResolution)]);
    end
end

% ============================================================
% FIGURE 1B: RUMOR-LAYER BASELINE TRAJECTORIES
% ============================================================

fig1B = figure('Color','w','Position',[120 90 820 620]);

hold on;
plot(Tplot, U(idxLog), 'Color', lineColors(5,:), 'LineWidth', 2.0);
plot(Tplot, A(idxLog), 'Color', lineColors(6,:), 'LineWidth', 2.0);
plot(Tplot, C(idxLog), 'Color', lineColors(7,:), 'LineWidth', 2.0);
yline(E_star(5),'--','Color',lineColors(5,:),'LineWidth',1.0);
yline(E_star(6),'--','Color',lineColors(6,:),'LineWidth',1.0);
yline(E_star(7),'--','Color',lineColors(7,:),'LineWidth',1.0);
xlabel('Normalized model time, t');
ylabel('Population fraction');
title('Rumor-layer baseline trajectories');
legend({'U(t)','A(t)','C(t)'}, 'Location','eastoutside', 'Box','off');
ylim([0 1]);
xlim(logTimeXLim);
grid on;
axis square;
box on;
set(gca, ...
    'FontName','Times New Roman', ...
    'LineWidth',1.05, ...
    'TickDir','out', ...
    'Layer','top', ...
    'XScale','log', ...
    'XTick',logTimeTickPositions, ...
    'XTickLabel',logTimeTickLabels, ...
    'XMinorTick','on', ...
    'YMinorTick','on');

if saveFigures
    if exist('exportgraphics','file') == 2
        exportgraphics(fig1B, fullfile(outDir,'Figure_1B_Rumor_Layer_Baseline_Trajectories.png'), ...
            'Resolution', exportResolution);
    else
        set(fig1B,'PaperPositionMode','auto');
        print(fig1B, fullfile(outDir,'Figure_1B_Rumor_Layer_Baseline_Trajectories.png'), ...
            '-dpng', ['-r',num2str(exportResolution)]);
    end
end

% ============================================================
% FINAL COMMAND-WINDOW DESCRIPTION
% ============================================================

fprintf('==============================================================\n');
fprintf('FIGURE DESCRIPTIONS\n');
fprintf('==============================================================\n');

fprintf('Figure 1A: Disease-layer baseline trajectories.\n');
fprintf('           A shifted logarithmic time axis x=t+1 is used so that\n');
fprintf('           the displayed tick labels are 0, 10^1, 10^2, 10^3, 10^4.\n');
fprintf('           Dashed horizontal levels mark the disease-layer components\n');
fprintf('           of the endemic equilibrium E*.\n\n');

fprintf('Figure 1B: Rumor-layer baseline trajectories.\n');
fprintf('           A shifted logarithmic time axis x=t+1 is used so that\n');
fprintf('           the displayed tick labels are 0, 10^1, 10^2, 10^3, 10^4.\n');
fprintf('           Dashed horizontal levels mark the rumor-layer components\n');
fprintf('           of the endemic equilibrium E*.\n\n');

fprintf('Stability information has been printed in the command window.\n');
fprintf('==============================================================\n');

if saveFigures
    fprintf('Both figures were saved at %d dpi in:\n%s\n', exportResolution, outDir);
else
    fprintf('Figures were generated but not saved. Set saveFigures = true to save them.\n');
end
fprintf('==============================================================\n');

% ============================================================
% END OF MAIN SCRIPT – LOCAL HELPER FUNCTIONS FOLLOW
% ============================================================

% -------------------------------------------------------------------------
function s = bool2str(b)
% BOOL2STR  Convert logical scalar into a readable conclusion.
if b
    s = 'STABLE / HOLDS';
else
    s = 'UNSTABLE / FAILS';
end
end

% -------------------------------------------------------------------------
function EQ = computeEq(p)
% COMPUTEEQ  Equilibria, reproduction numbers, and thresholds.
%
% Output structure EQ contains:
%   E0        rumor-free disease-free equilibrium
%   Er        rumor-present disease-free equilibrium
%   R0E       epidemic reproduction number at E0
%   R0R       rumor reproduction number
%   lambda_c  rumor threshold
%   beta_c0   disease threshold at E0
%   beta_cr   disease threshold at Er
%   psi_r     rumor-modified vaccination rate at Er
%
% If the rumor-present DFE does not exist, Er and beta_cr are NaN.

% Rumor-free disease-free equilibrium
S0 = p.Lambda*(p.omega + p.mu)/(p.mu*(p.psi + p.omega + p.mu));
V0 = p.Lambda*p.psi/(p.mu*(p.psi + p.omega + p.mu));
U0 = p.eta/p.mu;

E0 = [S0; 0; V0; 0; U0; 0; 0];

% Reproduction numbers
R0E = p.beta*(S0 + p.sigma*V0)/(p.gamma + p.mu);
R0R = p.lambda*U0/(p.delta + p.mu);

% Thresholds
lambda_c = p.mu*(p.delta + p.mu)/p.eta;

beta_c0 = p.mu*(p.gamma + p.mu)*(p.omega + p.mu + p.psi) ...
        /(p.Lambda*((p.omega + p.mu) + p.sigma*p.psi));

% Rumor-present disease-free equilibrium
if R0R > 1
    Ur = (p.delta + p.mu)/p.lambda;

    Ar = ((p.theta + p.mu)*(p.eta*p.lambda - p.mu*(p.delta + p.mu))) ...
       /(p.mu*p.lambda*(p.delta + p.theta + p.mu));

    Cr = p.delta*Ar/(p.theta + p.mu);

    psi_r = p.psi*(1 - p.alpha*Ar);

    Sr = p.Lambda*(p.omega + p.mu)/(p.mu*(psi_r + p.omega + p.mu));
    Vr = p.Lambda*psi_r/(p.mu*(psi_r + p.omega + p.mu));

    Er = [Sr; 0; Vr; 0; Ur; Ar; Cr];

    beta_cr = p.mu*(p.gamma + p.mu)*(p.omega + p.mu + psi_r) ...
            /(p.Lambda*((p.omega + p.mu) + p.sigma*psi_r));
else
    Er = NaN(7,1);
    psi_r = NaN;
    beta_cr = NaN;
end

EQ.E0 = E0;
EQ.Er = Er;
EQ.R0E = R0E;
EQ.R0R = R0R;
EQ.lambda_c = lambda_c;
EQ.beta_c0 = beta_c0;
EQ.beta_cr = beta_cr;
EQ.psi_r = psi_r;
end

% -------------------------------------------------------------------------
function J = jacobianAt(p, x)
% JACOBIANAT  Full 7x7 Jacobian of the coupled SVIR + rumor model.
% State order: [S; I; V; R; U; A; C]

S = x(1);
I = x(2);
V = x(3);
U = x(5);
A = x(6);

J = zeros(7,7);

% Row 1: dS/dt
J(1,1) = -p.beta*I - p.psi*(1 - p.alpha*A) - p.mu;
J(1,2) = -p.beta*S;
J(1,3) =  p.omega;
J(1,6) =  p.psi*p.alpha*S;

% Row 2: dI/dt
J(2,1) =  p.beta*I;
J(2,2) =  p.beta*S + p.sigma*p.beta*V - (p.gamma + p.mu);
J(2,3) =  p.sigma*p.beta*I;

% Row 3: dV/dt
J(3,1) =  p.psi*(1 - p.alpha*A);
J(3,2) = -p.sigma*p.beta*V;
J(3,3) = -p.sigma*p.beta*I - (p.omega + p.mu);
J(3,6) = -p.psi*p.alpha*S;

% Row 4: dR/dt
J(4,2) =  p.gamma;
J(4,4) = -p.mu;

% Row 5: dU/dt
J(5,2) = -p.phi*U;
J(5,5) = -p.lambda*A - p.phi*I - p.mu;
J(5,6) = -p.lambda*U;
J(5,7) =  p.theta;

% Row 6: dA/dt
J(6,2) =  p.phi*U;
J(6,5) =  p.lambda*A + p.phi*I;
J(6,6) =  p.lambda*U - (p.delta + p.mu);

% Row 7: dC/dt
J(7,6) =  p.delta;
J(7,7) = -(p.theta + p.mu);
end

% -------------------------------------------------------------------------
function dy = rhsModel(~, y, p)
% RHSMODEL  Right-hand side of the coupled SVIR + rumor ODE system.
% State order: [S; I; V; R; U; A; C]

S = y(1);
I = y(2);
V = y(3);
R = y(4); %#ok<NASGU>
U = y(5);
A = y(6);
C = y(7);

dy = zeros(7,1);

dy(1) = p.Lambda ...
      - p.beta*S*I ...
      - p.psi*(1 - p.alpha*A)*S ...
      + p.omega*V ...
      - p.mu*S;

dy(2) = p.beta*S*I ...
      + p.sigma*p.beta*V*I ...
      - (p.gamma + p.mu)*I;

dy(3) = p.psi*(1 - p.alpha*A)*S ...
      - p.sigma*p.beta*V*I ...
      - (p.omega + p.mu)*V;

dy(4) = p.gamma*I ...
      - p.mu*y(4);

dy(5) = p.eta ...
      - p.lambda*U*A ...
      - p.phi*I*U ...
      + p.theta*C ...
      - p.mu*U;

dy(6) = p.lambda*U*A ...
      + p.phi*I*U ...
      - (p.delta + p.mu)*A;

dy(7) = p.delta*A ...
      - (p.theta + p.mu)*C;
end