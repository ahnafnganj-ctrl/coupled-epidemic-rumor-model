function out = coupling_comparison()
% COUPLING_COMPARISON_R2020A Reproduce Table 3 and plot its four cases.
% MATLAB R2020a; base MATLAB only (no additional toolboxes).
% Save as coupling_comparison_R2020a.m and run:
%     out = coupling_comparison_R2020a;
%
% Equilibrium fractions are computed from the steady-state equations.
% Imax is computed over [0,1000], using peak events and both endpoints.
% Integration starts at t=0; the logarithmic plots start at t=0.01.
% This function displays one figure and does not save any files.

%% Baseline parameters and initial condition
p0 = struct('Lambda',0.02, 'eta',0.02, 'mu',0.02, ...
    'beta',0.70, 'gamma',0.25, 'psi',0.30, 'alpha',0.50, ...
    'omega',0.05, 'sigma',0.40, 'lambda',0.25, 'phi',0.60, ...
    'delta',0.15, 'theta',0.10);

% State order: S, I, V, R, U, A, C.
X0 = [0.89; 0.10; 0.01; 0; 0.90; 0.10; 0];
T = 1000;
tPlot = logspace(-2,log10(T),1800);

caseNames = {'Full baseline'; 'alpha = 0'; ...
             'phi = 0'; 'alpha = phi = 0'};

% Unicode Greek letters allow all title text to use Arial, without LaTeX.
caseTitles = {'Full baseline'; [char(945) ' = 0']; ...
    [char(966) ' = 0']; [char(945) ' = ' char(966) ' = 0']};

parameters = repmat(p0,4,1);
parameters(2).alpha = 0;
parameters(3).phi = 0;
parameters(4).alpha = 0;
parameters(4).phi = 0;

values = zeros(4,4);       % Columns: Istar, Astar, Vstar, Imax.
equilibria = zeros(7,4);
reproduction = zeros(4,3); % Columns: R0E, R0R, RE1.
peakTimes = zeros(4,1);
residuals = zeros(4,1);
terminalErrors = zeros(4,1);
balanceErrors = zeros(4,1);
solutions = cell(4,1);
trajectories = cell(4,1);

%% Compute equilibria, trajectories, peaks, and reproduction numbers
for k = 1:4
    p = parameters(k);

    Xstar = endemic_equilibrium(p);
    equilibria(:,k) = Xstar;
    residuals(k) = norm(model_rhs(0,Xstar,p),inf);

    opts = odeset('RelTol',1e-10, 'AbsTol',1e-12, ...
        'MaxStep',0.5, 'NonNegative',1:7, ...
        'Events',@(t,X) infection_peak_event(t,X,p));

    sol = ode45(@(t,X) model_rhs(t,X,p),[0 T],X0,opts);

    Xplot = deval(sol,tPlot);
    XT = deval(sol,T);
    solutions{k} = sol;
    trajectories{k} = Xplot;

    % Check both endpoints and all detected local infection maxima.
    candidateTimes = [0; T];
    candidatePeaks = [X0(2); XT(2)];

    if isfield(sol,'xe') && ~isempty(sol.xe)
        candidateTimes = [candidateTimes; sol.xe(:)];
        eventInfection = sol.ye(2,:);
        candidatePeaks = [candidatePeaks; eventInfection(:)];
    end

    [Imax,idx] = max(candidatePeaks);
    peakTimes(k) = candidateTimes(idx);

    values(k,:) = [Xstar(2),Xstar(6),Xstar(3),Imax];
    reproduction(k,:) = reproduction_numbers(p);

    terminalErrors(k) = norm(XT-Xstar,inf);
    balanceErrors(k) = max([ ...
        max(abs(sum(Xplot(1:4,:),1)-p.Lambda/p.mu)), ...
        max(abs(sum(Xplot(5:7,:),1)-p.eta/p.mu))]);
end

Table3 = table(caseNames,values(:,1),values(:,2),values(:,3), ...
    values(:,4),'VariableNames',{'Case','Istar','Astar','Vstar','Imax'});

ReproductionNumbers = table(caseNames,reproduction(:,1), ...
    reproduction(:,2),reproduction(:,3), ...
    'VariableNames',{'Case','R0E','R0R','RE1'});

%% Compute the manuscript percentages from unrounded values
stats.I_increase_vs_alpha0 = 100*(values(1,1)/values(2,1)-1);

stats.Imax_increase_vs_alpha0 = 100*(values(1,4)/values(2,4)-1);

stats.V_decrease_vs_alpha0 = 100*(1-values(1,3)/values(2,3));

stats.I_increase_vs_phi0 = 100*(values(1,1)/values(3,1)-1);

stats.A_increase_vs_phi0 = 100*(values(1,2)/values(3,2)-1);

stats.RE1_increase_over_R0E = ...
    100*(reproduction(1,3)/reproduction(1,1)-1);

%% Print the table and requested statistics
fprintf('\nTABLE 3: COUPLING COMPARISON\n');
fprintf('Baseline: beta=0.70, lambda=0.25, alpha=0.50, phi=0.60.\n');
fprintf('Istar, Astar and Vstar are equilibrium fractions.\n');
fprintf('Imax = max I(t) over 0 <= t <= %.0f.\n\n',T);

fprintf('%-20s %12s %12s %12s %12s\n', ...
    'Case','Istar','Astar','Vstar','Imax');

for k = 1:4
    fprintf('%-20s %12.6f %12.6f %12.6f %12.6f\n', ...
        caseNames{k},values(k,:));
end

fprintf('\nREPRODUCTION NUMBERS\n');
fprintf('%-20s %14s %14s %14s\n','Case','R0E','R0R','R_E|1');

for k = 1:4
    fprintf('%-20s %14.8f %14.8f %14.8f\n', ...
        caseNames{k},reproduction(k,:));
end

fprintf('\nFULL BASELINE RELATIVE TO alpha = 0\n');

fprintf('Endemic infection increase:    %.4f %% (%.2f %%)\n', ...
    stats.I_increase_vs_alpha0,stats.I_increase_vs_alpha0);

fprintf('Peak infection increase:       %.4f %% (%.2f %%)\n', ...
    stats.Imax_increase_vs_alpha0,stats.Imax_increase_vs_alpha0);

fprintf('Endemic vaccination decrease:  %.4f %% (%.2f %%)\n', ...
    stats.V_decrease_vs_alpha0,stats.V_decrease_vs_alpha0);

fprintf('\nFULL BASELINE RELATIVE TO phi = 0\n');

fprintf('Endemic infection increase:    %.4f %% (%.2f %%)\n', ...
    stats.I_increase_vs_phi0,stats.I_increase_vs_phi0);

fprintf('Active-rumour increase:        %.4f %% (%.2f %%)\n', ...
    stats.A_increase_vs_phi0,stats.A_increase_vs_phi0);

fprintf('\nDISEASE INVASION COMPARISON\n');

fprintf('R0E at E0:                    %.8f\n',reproduction(1,1));

fprintf('R_E|1 at E1, full baseline:    %.8f\n',reproduction(1,3));

fprintf('R_E|1 at E1, alpha = 0:        %.8f\n',reproduction(2,3));

fprintf('Increase from R0E to R_E|1:    %.4f %% (%.2f %%)\n', ...
    stats.RE1_increase_over_R0E,stats.RE1_increase_over_R0E);

fprintf('R0R > 1 in all four cases; R_E|1 uses the persistent-rumour E1.\n');

fprintf('\nNUMERICAL CHECKS\n');

fprintf('%-20s %11s %15s %17s\n', ...
    'Case','Peak time','Eq. residual','|X(T)-Xstar|_inf');

for k = 1:4
    fprintf('%-20s %11.6f %15.3e %17.3e\n', ...
        caseNames{k},peakTimes(k),residuals(k),terminalErrors(k));
end

fprintf('Largest population-balance error: %.3e\n',max(balanceErrors));

%% Publication-style 2-by-2 figure: I(t), A(t), and V(t) only
colours = [0.0000 0.4470 0.6980; ...  % I: blue
           0.8350 0.3690 0.0000; ...  % A: orange
           0.0000 0.6200 0.4510];    % V: green

styles = {'-','--','-.'};
stateIndices = [2 6 3];

fig = figure('Color','w','Units','centimeters', ...
    'Position',[2 2 34 24],'Renderer','painters', ...
    'Name','Coupling comparison','NumberTitle','off');

axesHandles = gobjects(4,1);

for k = 1:4
    ax = subplot(2,2,k,'Parent',fig);
    axesHandles(k) = ax;
    hold(ax,'on');

    h = gobjects(3,1);

    for j = 1:3
        h(j) = plot(ax,tPlot,trajectories{k}(stateIndices(j),:), ...
            'Color',colours(j,:), ...
            'LineStyle',styles{j}, ...
            'LineWidth',2.5);
    end

    set(ax, ...
        'XScale','log', ...
        'XLim',[tPlot(1) T], ...
        'XTick',10.^(-2:3), ...
        'YLim',[0 0.70], ...
        'YTick',0:0.1:0.7, ...
        'FontName','Arial', ...
        'FontSize',16, ...
        'FontWeight','normal', ...
        'LineWidth',1.2, ...
        'TickDir','out', ...
        'TickLength',[0.015 0.015], ...
        'Box','on', ...
        'Layer','top', ...
        'XMinorTick','on', ...
        'YMinorTick','off', ...
        'XGrid','on', ...
        'YGrid','on', ...
        'XMinorGrid','off', ...
        'YMinorGrid','off', ...
        'GridAlpha',0.12, ...
        'GridColor',[0.35 0.35 0.35]);

    title(ax,caseTitles{k}, ...
        'FontName','Arial', ...
        'FontSize',21, ...
        'FontWeight','bold', ...
        'Interpreter','none');

    xlabel(ax,'Time, t', ...
        'FontName','Arial', ...
        'FontSize',18, ...
        'FontWeight','bold', ...
        'Interpreter','none');

    ylabel(ax,'Population fraction', ...
        'FontName','Arial', ...
        'FontSize',18, ...
        'FontWeight','bold', ...
        'Interpreter','none');

    % A single vertical legend, inside the first subplot.
    if k == 1
        lgd = legend(ax,h,{'I(t)','A(t)','V(t)'}, ...
            'Location','best', ...
            'Orientation','vertical', ...
            'FontName','Arial', ...
            'FontSize',16, ...
            'Interpreter','none');

        set(lgd,'Box','on','Color','w', ...
            'EdgeColor',[0.75 0.75 0.75]);
    end

    hold(ax,'off');
end

linkaxes(axesHandles,'xy');
drawnow;

%% Return results without writing to disk
out.parameters = parameters;
out.initialState = X0;
out.Table3 = Table3;
out.ReproductionNumbers = ReproductionNumbers;
out.percentages = stats;
out.equilibria = equilibria;
out.peakTimes = peakTimes;
out.equilibriumResiduals = residuals;
out.terminalErrors = terminalErrors;
out.populationBalanceErrors = balanceErrors;
out.tPlot = tPlot;
out.trajectories = trajectories;
out.solutions = solutions;
out.figure = fig;

end

%% Model equations
function F = model_rhs(~,X,p)

S=X(1); I=X(2); V=X(3); R=X(4);
U=X(5); A=X(6); C=X(7);

q = p.psi*(1-p.alpha*A);
activation = (p.lambda*A+p.phi*I)*U;

F = [p.Lambda-p.beta*S*I-q*S+p.omega*V-p.mu*S; ...
     p.beta*(S+p.sigma*V)*I-(p.gamma+p.mu)*I; ...
     q*S-p.sigma*p.beta*V*I-(p.omega+p.mu)*V; ...
     p.gamma*I-p.mu*R; ...
     p.eta-activation+p.theta*C-p.mu*U; ...
     activation-(p.delta+p.mu)*A; ...
     p.delta*A-(p.theta+p.mu)*C];

end

%% Detect local maxima of infection
function [value,isterminal,direction] = infection_peak_event(t,X,p)

F = model_rhs(t,X,p);

value = F(2);        % dI/dt = 0.
isterminal = 0;      % Continue to t=1000.
direction = -1;     % Positive-to-negative crossing: local maximum.

end

%% Compute the endemic equilibrium
function Xstar = endemic_equilibrium(p)
% Scalar reduction of the original steady-state equations.
% Valid for all four manuscript cases, including alpha=0 and phi=0.

g = p.gamma+p.mu;
fun = @(I) infection_equilibrium_residual(I,p);
Iupper = p.Lambda/g;

assert(fun(0)>0 && fun(Iupper)<0, ...
    'The endemic-root bracket is not valid for these parameters.');

Istar = fzero(fun,[0 Iupper], ...
    optimset('TolX',1e-13,'Display','off'));

[S,V,U,A,C] = equilibrium_components(Istar,p);

Xstar = [S; Istar; V; p.gamma*Istar/p.mu; U; A; C];

assert(all(isfinite(Xstar)) && all(Xstar>=0), ...
    'The reconstructed equilibrium is not admissible.');

assert(norm(model_rhs(0,Xstar,p),inf)<1e-10, ...
    'The equilibrium residual is too large.');

end

function f = infection_equilibrium_residual(I,p)

[S,V] = equilibrium_components(I,p);

f = p.beta*(S+p.sigma*V)-(p.gamma+p.mu);

end

function [S,V,U,A,C] = equilibrium_components(I,p)
% At equilibrium: C=delta*A/r and U=M-kappa*A.
% Hence kappa*lambda*A^2 + B*A - phi*I*M = 0.

d = p.delta+p.mu;
r = p.theta+p.mu;
M = p.eta/p.mu;
kappa = 1+p.delta/r;

assert(p.lambda*M>d, ...
    'This Table 3 calculation uses the persistent-rumour regime.');

B = kappa*p.phi*I+d-p.lambda*M;
D = sqrt(B^2+4*kappa*p.lambda*p.phi*I*M);

if B >= 0
    A = 2*p.phi*I*M/(B+D); % Stable evaluation of the positive root.
else
    A = (D-B)/(2*kappa*p.lambda);
end

C = p.delta*A/r;
U = M-kappa*A;
q = p.psi*(1-p.alpha*A);

assert(q>=0 && U>=0,'Inadmissible rumour or vaccination state.');

% Solve dS/dt=dV/dt=0 for fixed I and A.
denV = p.sigma*p.beta*I+p.omega+p.mu;

den = (p.beta*I+p.mu)*denV ...
    + q*(p.sigma*p.beta*I+p.mu);

S = p.Lambda*denV/den;
V = p.Lambda*q/den;

end

%% Reproduction numbers at the two disease-free backgrounds
function nums = reproduction_numbers(p)

N = p.Lambda/p.mu;
M = p.eta/p.mu;
g = p.gamma+p.mu;
d = p.delta+p.mu;
r = p.theta+p.mu;

R0E = p.beta*N*(p.omega+p.mu+p.sigma*p.psi) ...
    /(g*(p.psi+p.omega+p.mu));

R0R = p.lambda*M/d;

assert(R0R>1, ...
    'The rumour-persistent equilibrium E1 requires R0R>1.');

A1 = r*(p.lambda*M-d)/(p.lambda*(d+p.theta));
q1 = p.psi*(1-p.alpha*A1);

RE1 = p.beta*N*(p.omega+p.mu+p.sigma*q1) ...
    /(g*(q1+p.omega+p.mu));

nums = [R0E R0R RE1];

end