%% FULL OPTIMAL-CONTROL STUDY FOR THE COUPLED EPIDEMIC--RUMOR MODEL
clear; close all; clc;
%% 1. Output folder, parameters, weights, bounds, and initial condition
outdir = fullfile(pwd,'full_optimal_control_matlab_results');
if ~exist(outdir,'dir'), mkdir(outdir); end

p.Lambda = 0.02; p.mu    = 0.02; p.eta   = 0.02;
p.beta   = 0.70; p.gamma = 0.25; p.psi   = 0.30;
p.alpha  = 0.50; p.omega = 0.05; p.sigma = 0.40;
p.lam    = 0.25; p.phi   = 0.60; p.delta = 0.15;
p.theta  = 0.10;

% w = [B_I, B_A, C_1, C_2, C_3].
w = [10, 4, 12, 4, 8];
umax = [0.50; 0.50; 0.50];
x0 = [0.89; 0.10; 0.01; 0.00; 0.90; 0.10; 0.00];
T = 100; N = 2001; t = linspace(0,T,N);
relax = 0.20; tol = 2e-7; maxIter = 180;

% The exact sufficient condition requested for psi*(1-alpha*A) >= 0.
positivityMargin = 1 - p.alpha*p.eta/p.mu;
assert(p.sigma >= 0 && p.sigma <= 1, ...
    'The vaccine leakiness sigma must lie in [0,1].');
assert(positivityMargin >= -100*eps, ...
    'The condition alpha*eta/mu <= 1 is violated.');
fprintf('Positivity check: alpha*eta/mu = %.8f <= 1; margin = %.8f.\n', ...
    p.alpha*p.eta/p.mu, positivityMargin);

%% 2. Main optimal solution and no-control comparator
[x,z,u,history] = solve_sweep(t,x0,p,w,umax,relax,tol,maxIter);
xNo = forward_states(t,x0,zeros(3,N),p);
raw = unconstrained_controls(x,z,w);
projected = project_controls(raw,umax);
H = hamiltonian_path(x,z,u,p,w);

m  = compute_metrics(t,x,u,p,w);
m0 = compute_metrics(t,xNo,zeros(3,N),p,w);

Nexact = p.Lambda/p.mu + (sum(x0(1:4))-p.Lambda/p.mu)*exp(-p.mu*t);
Mexact = p.eta/p.mu + (sum(x0(5:7))-p.eta/p.mu)*exp(-p.mu*t);
checks.iterations = numel(history);
checks.lastControlChange = history(end);
checks.projectionResidual = max(abs(u(:)-projected(:)));
checks.transversalityResidual = max(abs(z(:,end)));
checks.epidemicMassError = max(abs(sum(x(1:4,:),1)-Nexact));
checks.rumorMassError = max(abs(sum(x(5:7,:),1)-Mexact));
checks.minimumState = min(x(:));
checks.hamiltonianRange = max(H)-min(H);
checks.hamiltonianRelativeRange = checks.hamiltonianRange/max(abs(mean(H)),eps);

fprintf('\nMAIN OPTIMAL SOLUTION\n');
fprintf('Iterations                         : %d\n',checks.iterations);
fprintf('Final maximum control change       : %.9e\n',checks.lastControlChange);
fprintf('Projection residual                : %.9e\n',checks.projectionResidual);
fprintf('Terminal-costate residual          : %.9e\n',checks.transversalityResidual);
fprintf('Epidemic population-balance error  : %.9e\n',checks.epidemicMassError);
fprintf('Rumor population-balance error     : %.9e\n',checks.rumorMassError);
fprintf('Minimum computed state             : %.9e\n',checks.minimumState);
fprintf('Hamiltonian relative range         : %.9e\n',checks.hamiltonianRelativeRange);
fprintf('Objective J: no control %.9f -> optimal %.9f (%.4f%% reduction)\n', ...
    m0.J,m.J,100*(1-m.J/m0.J));
fprintf('Integral I: no control %.9f -> optimal %.9f (%.4f%% reduction)\n', ...
    m0.cumI,m.cumI,100*(1-m.cumI/m0.cumI));
fprintf('Integral A: no control %.9f -> optimal %.9f (%.4f%% reduction)\n', ...
    m0.cumA,m.cumA,100*(1-m.cumA/m0.cumA));

%% 3. Figure 1: all seven state trajectories
fig = figure('Color','w','Position',[80 80 1120 760]);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
diseaseNames = {'S','I','V','R'}; rumorNames = {'U','A','C'};
cc = [0.043 0.447 0.522; 0.608 0.110 0.110; 0.427 0.157 0.851; 0.851 0.467 0.024];
nexttile; hold on
for j=1:4, plot(t,xNo(j,:),'LineWidth',1.8,'Color',cc(j,:)); end
title('Disease layer: no control'); xlabel('normalized time'); ylabel('fraction');
legend(diseaseNames,'Location','best','Box','off'); grid on
nexttile; hold on
for j=1:4, plot(t,x(j,:),'LineWidth',1.8,'Color',cc(j,:)); end
title('Disease layer: optimal control'); xlabel('normalized time'); ylabel('fraction');
legend(diseaseNames,'Location','best','Box','off'); grid on
nexttile; hold on
for j=1:3, plot(t,xNo(4+j,:),'LineWidth',1.8,'Color',cc(j,:)); end
title('Rumor layer: no control'); xlabel('normalized time'); ylabel('fraction');
legend(rumorNames,'Location','best','Box','off'); grid on
nexttile; hold on
for j=1:3, plot(t,x(4+j,:),'LineWidth',1.8,'Color',cc(j,:)); end
title('Rumor layer: optimal control'); xlabel('normalized time'); ylabel('fraction');
legend(rumorNames,'Location','best','Box','off'); grid on
title(tl,'All state trajectories');
exportgraphics(fig,fullfile(outdir,'01_all_state_trajectories.png'),'Resolution',240);

%% 4. Figure 2: primary disease and rumor outcomes
fig = figure('Color','w','Position',[80 80 1120 430]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(t,xNo(2,:),'--','LineWidth',2.2,'Color',cc(2,:)); hold on
plot(t,x(2,:),'LineWidth',2.2,'Color',cc(1,:));
xlabel('normalized time'); ylabel('I(t)'); title('Infected fraction'); grid on
legend('no control','optimal control','Location','best','Box','off');
nexttile; plot(t,xNo(6,:),'--','LineWidth',2.2,'Color',cc(2,:)); hold on
plot(t,x(6,:),'LineWidth',2.2,'Color',cc(3,:));
xlabel('normalized time'); ylabel('A(t)'); title('Active-rumor fraction'); grid on
legend('no control','optimal control','Location','best','Box','off');
exportgraphics(fig,fullfile(outdir,'02_primary_outcomes.png'),'Resolution',240);

%% 5. Figure 3: optimal control profiles
fig = figure('Color','w','Position',[80 80 980 470]); hold on
plot(t,u(1,:),'LineWidth',2.2,'Color',cc(1,:));
plot(t,u(2,:),'LineWidth',2.2,'Color',cc(4,:));
plot(t,u(3,:),'LineWidth',2.2,'Color',cc(3,:));
xlabel('normalized time'); ylabel('control intensity'); ylim([-0.015 0.52]);
title('Optimal control profiles'); grid on
legend('u_1: vaccination','u_2: counter-rumor','u_3: treatment','Location','best','Box','off');
exportgraphics(fig,fullfile(outdir,'03_control_profiles.png'),'Resolution',240);

%% 6. Figure 4: unconstrained expressions and projected controls
fig = figure('Color','w','Position',[80 80 1160 410]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
controlNames = {'u_1: vaccination','u_2: counter-rumor','u_3: treatment'};
controlColors = [cc(1,:);cc(4,:);cc(3,:)];
for j=1:3
    nexttile; plot(t,raw(j,:),'LineWidth',1.6,'Color',[0.45 0.50 0.58]); hold on
    plot(t,u(j,:),'LineWidth',2.1,'Color',controlColors(j,:));
    yline(0,'k-'); yline(umax(j),'k--'); grid on
    xlabel('time'); ylabel('intensity'); title(controlNames{j});
end
legend('unconstrained expression','projected optimum','Location','best','Box','off');
exportgraphics(fig,fullfile(outdir,'04_projection_and_switching.png'),'Resolution',240);

%% 7. Figure 5: all costate trajectories
fig = figure('Color','w','Position',[80 80 1050 650]);
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile; hold on
for j=1:4, plot(t,z(j,:),'LineWidth',1.7,'Color',cc(j,:)); end
ylabel('costate value'); title('Disease-layer costates'); grid on
legend('\xi_S','\xi_I','\xi_V','\xi_R','Location','best','Box','off');
nexttile; hold on
for j=1:3, plot(t,z(4+j,:),'LineWidth',1.7,'Color',cc(j,:)); end
xlabel('normalized time'); ylabel('costate value'); title('Rumor-layer costates'); grid on
legend('\xi_U','\xi_A','\xi_C','Location','best','Box','off');
exportgraphics(fig,fullfile(outdir,'05_costate_trajectories.png'),'Resolution',240);

%% 8. Figure 6: mechanism diagnostics
gINo = p.beta*(xNo(1,:)+p.sigma*xNo(3,:))-(p.gamma+p.mu);
gI   = p.beta*(x(1,:)+p.sigma*x(3,:))-(p.gamma+u(3,:)+p.mu);
rLinearNo = p.lam*xNo(5,:)-(p.delta+p.mu);
rLinear   = p.lam*x(5,:)-(p.delta+u(2,:)+p.mu);
rSourceNo = p.phi*xNo(2,:).*xNo(5,:);
rSource   = p.phi*x(2,:).*x(5,:);
fig = figure('Color','w','Position',[80 80 1120 430]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(t,gINo,'--','LineWidth',2,'Color',cc(2,:)); hold on
plot(t,gI,'LineWidth',2,'Color',cc(1,:)); yline(0,'k-'); grid on
xlabel('time'); ylabel('rate'); title('\beta(S+\sigmaV)-(\gamma+u_3+\mu)');
legend('no control','optimal control','Location','best','Box','off');
nexttile; plot(t,rLinearNo,'--','LineWidth',2,'Color',cc(2,:)); hold on
plot(t,rLinear,'LineWidth',2,'Color',cc(4,:));
plot(t,rSourceNo,':','LineWidth',2,'Color',[0.45 0.50 0.58]);
plot(t,rSource,'LineWidth',1.8,'Color',cc(3,:)); yline(0,'k-'); grid on
xlabel('time'); ylabel('rate or source'); title('Active-rumor mechanisms');
legend('linear rate, no control','linear rate, controlled', ...
       '\phi IU, no control','\phi IU, controlled','Location','best','Box','off');
exportgraphics(fig,fullfile(outdir,'06_mechanism_diagnostics.png'),'Resolution',240);

%% 9. Figure 7: objective decomposition
BI=w(1); BA=w(2); C1=w(3); C2=w(4); C3=w(5);
inst = [BI*x(2,:); BA*x(6,:); 0.5*C1*u(1,:).^2; ...
        0.5*C2*u(2,:).^2; 0.5*C3*u(3,:).^2];
componentNames = {'infection burden','rumor burden','vaccination cost', ...
                  'counter-rumor cost','treatment cost'};
componentColors = [cc(2,:);cc(4,:);cc(1,:);0.15 0.39 0.92;cc(3,:)];
fig = figure('Color','w','Position',[80 80 1120 440]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on
for j=1:5
    plot(t,cumtrapz(t,inst(j,:)),'LineWidth',2,'Color',componentColors(j,:));
end
xlabel('time'); ylabel('accumulated contribution');
title('Cumulative objective components'); grid on
legend(componentNames,'Location','best','Box','off');

% -- Updated bar chart with distinct colors (Figure 7, right panel) --
nexttile; 
bh = barh(m.Jparts, 'FaceColor','flat');
% Assign each bar its own color (same order as componentColors)
for k = 1:5
    bh.CData(k,:) = componentColors(k,:);
end
set(gca,'YTick',1:5,'YTickLabel',componentNames);
xlabel('contribution to J'); title('Final objective decomposition'); grid on
% -----------------------------------------------------------------

exportgraphics(fig,fullfile(outdir,'07_objective_decomposition.png'),'Resolution',240);

%% 10. Figure 8: integrated transition flows
% -- Updated bar chart with distinct colors --
flowNames = {'baseline vaccination','additional vaccination', ...
             'counter-rumor conversion','treatment recovery'};
flowColors = [0.043 0.447 0.522; 0.851 0.467 0.024; 0.608 0.110 0.110; 0.427 0.157 0.851];
fig = figure('Color','w','Position',[80 80 920 440]);
bh = bar(m.flows, 'FaceColor','flat');
for k = 1:4
    bh.CData(k,:) = flowColors(k,:);
end
grid on
set(gca,'XTickLabel',flowNames);
xtickangle(15); ylabel('integrated flow over [0,T]');
title('Integrated controlled transition flows');
% --------------------------------------------
exportgraphics(fig,fullfile(outdir,'08_integrated_flows.png'),'Resolution',240);

%% 11. Figure 9: convergence and autonomous-Hamiltonian check
fig = figure('Color','w','Position',[80 80 1120 410]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; semilogy(1:numel(history),history,'LineWidth',2,'Color',cc(1,:)); hold on
yline(tol,'k--'); xlabel('iteration'); ylabel('maximum control change');
title('Forward-backward convergence'); grid on
nexttile; plot(t,H,'LineWidth',1.8,'Color',cc(3,:)); hold on
yline(mean(H),'k--'); xlabel('time'); ylabel('H(t)');
title('Hamiltonian along the computed extremal'); grid on
exportgraphics(fig,fullfile(outdir,'09_numerical_optimality_checks.png'),'Resolution',240);

%% 12. Strategy comparison: none, single, paired, and all controls
strategyNames = {'none','u1','u2','u3','u1+u2','u1+u3','u2+u3','all'};
strategyBounds = [0 0 0; .5 0 0; 0 .5 0; 0 0 .5; ...
                  .5 .5 0; .5 0 .5; 0 .5 .5; .5 .5 .5];
strategyMetrics = repmat(empty_metrics(),1,numel(strategyNames));
strategyIterations = zeros(1,numel(strategyNames));
for s=1:numel(strategyNames)
    if s==1
        xs=xNo; us=zeros(3,N); hs=0;
    elseif s==numel(strategyNames)
        xs=x; us=u; hs=history;
    else
        [xs,~,us,hs]=solve_sweep(t,x0,p,w,strategyBounds(s,:)', ...
                                  relax,tol,maxIter);
    end
    strategyMetrics(s)=compute_metrics(t,xs,us,p,w);
    strategyIterations(s)=numel(hs);
end
redI = 100*(1-[strategyMetrics.cumI]/m0.cumI);
redA = 100*(1-[strategyMetrics.cumA]/m0.cumA);
fig = figure('Color','w','Position',[80 80 1030 470]);
bar([redI(:),redA(:)]); grid on; ylim([min(-2,min([redI redA])-3),100]);
set(gca,'XTickLabel',strategyNames); ylabel('reduction relative to no control (%)');
title('Single, paired, and combined control strategies');
legend('cumulative infection','cumulative active rumor','Location','best','Box','off');
exportgraphics(fig,fullfile(outdir,'10_strategy_comparison.png'),'Resolution',240);

fprintf('\nSTRATEGY COMPARISON\n');
fprintf('%-8s %12s %12s %12s %12s %12s\n', ...
    'strategy','J','int I','int A','red I (%)','red A (%)');
for s=1:numel(strategyNames)
    fprintf('%-8s %12.6f %12.6f %12.6f %12.4f %12.4f\n', ...
        strategyNames{s},strategyMetrics(s).J,strategyMetrics(s).cumI, ...
        strategyMetrics(s).cumA,redI(s),redA(s));
end

%% 13. Two-dimensional objective-weight study
BIvals=[5 10 20]; BAvals=[2 4 8];
weightCumI=zeros(3); weightCumA=zeros(3); weightJ=zeros(3); weightEffort=zeros(3);
for ib=1:numel(BIvals)
    for ia=1:numel(BAvals)
        ws=[BIvals(ib),BAvals(ia),12,4,8];
        ts=linspace(0,T,601);
        [xs,~,us,~]=solve_sweep(ts,x0,p,ws,umax,relax,tol,maxIter);
        ms=compute_metrics(ts,xs,us,p,ws);
        weightCumI(ib,ia)=ms.cumI; weightCumA(ib,ia)=ms.cumA;
        weightJ(ib,ia)=ms.J; weightEffort(ib,ia)=sum(ms.L2sq);
    end
end
fig=figure('Color','w','Position',[80 80 1000 760]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
axh=nexttile; plot_heatmap(axh,weightCumI,BAvals,BIvals,'integral I dt');
axh=nexttile; plot_heatmap(axh,weightCumA,BAvals,BIvals,'integral A dt');
axh=nexttile; plot_heatmap(axh,weightJ,BAvals,BIvals,'objective J');
axh=nexttile; plot_heatmap(axh,weightEffort,BAvals,BIvals,'total integral u_j^2 dt');
title(tl,'Objective-weight sensitivity');
exportgraphics(fig,fullfile(outdir,'11_weight_grid.png'),'Resolution',240);

%% 14. Control-cost, common-bound, and horizon sensitivity
costScale=[0.5 1 2]; costI=zeros(size(costScale)); costA=zeros(size(costScale));
for j=1:numel(costScale)
    ws=[w(1:2),costScale(j)*w(3:5)]; ts=linspace(0,T,701);
    [xs,~,us,~]=solve_sweep(ts,x0,p,ws,umax,relax,tol,maxIter);
    ms=compute_metrics(ts,xs,us,p,ws); costI(j)=ms.cumI; costA(j)=ms.cumA;
end
commonBound=[.25 .50 .75]; boundI=zeros(size(commonBound)); boundA=zeros(size(commonBound));
for j=1:numel(commonBound)
    ts=linspace(0,T,701); ub=commonBound(j)*ones(3,1);
    [xs,~,us,~]=solve_sweep(ts,x0,p,w,ub,relax,tol,maxIter);
    ms=compute_metrics(ts,xs,us,p,w); boundI(j)=ms.cumI; boundA(j)=ms.cumA;
end
horizons=[50 100 150]; horizonI=zeros(size(horizons)); horizonA=zeros(size(horizons));
for j=1:numel(horizons)
    Ts=horizons(j); ts=linspace(0,Ts,10*Ts+1);
    [xs,~,us,~]=solve_sweep(ts,x0,p,w,umax,relax,tol,maxIter);
    ms=compute_metrics(ts,xs,us,p,w); horizonI(j)=ms.cumI; horizonA(j)=ms.cumA;
end
fig=figure('Color','w','Position',[80 80 1240 420]);
tl=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; plot(costScale,costI,'o-','LineWidth',2,'Color',cc(1,:)); hold on
plot(costScale,costA,'s-','LineWidth',2,'Color',cc(4,:)); grid on
xlabel('control-cost multiplier'); ylabel('cumulative burden'); title('Cost scale');
legend('integral I dt','integral A dt','Location','best','Box','off');
nexttile; plot(commonBound,boundI,'o-','LineWidth',2,'Color',cc(1,:)); hold on
plot(commonBound,boundA,'s-','LineWidth',2,'Color',cc(4,:)); grid on
xlabel('common upper bound'); ylabel('cumulative burden'); title('Bound');
legend('integral I dt','integral A dt','Location','best','Box','off');
nexttile; plot(horizons,horizonI,'o-','LineWidth',2,'Color',cc(1,:)); hold on
plot(horizons,horizonA,'s-','LineWidth',2,'Color',cc(4,:)); grid on
xlabel('time horizon'); ylabel('cumulative burden'); title('Horizon');
legend('integral I dt','integral A dt','Location','best','Box','off');
exportgraphics(fig,fullfile(outdir,'12_cost_bound_horizon_sensitivity.png'),'Resolution',240);

%% 15. Initial-condition robustness
% Changed 'manuscript baseline' to 'baseline'
initialNames={'low I / low A','high I / low A','low I / high A', ...
              'high I / high A','baseline'};
initialStates=[.93 .02 .05 0 .98 .02 0; ...
               .75 .20 .05 0 .98 .02 0; ...
               .93 .02 .05 0 .75 .25 0; ...
               .75 .20 .05 0 .75 .25 0; ...
               x0'];
initialRedI=zeros(1,5); initialRedA=zeros(1,5);
for j=1:5
    ts=linspace(0,T,701); xinit=initialStates(j,:)';
    [xs,~,us,~]=solve_sweep(ts,xinit,p,w,umax,relax,tol,maxIter);
    xn=forward_states(ts,xinit,zeros(3,numel(ts)),p);
    ms=compute_metrics(ts,xs,us,p,w); mn=compute_metrics(ts,xn,zeros(3,numel(ts)),p,w);
    initialRedI(j)=100*(1-ms.cumI/mn.cumI);
    initialRedA(j)=100*(1-ms.cumA/mn.cumA);
end
fig=figure('Color','w','Position',[80 80 1050 490]);
bar([initialRedI(:),initialRedA(:)]); grid on; ylim([0 100]);
set(gca,'XTickLabel',initialNames); xtickangle(15);
ylabel('reduction under optimal control (%)');
title('Robustness to initial infection and rumor prevalence');
legend('infection burden','active-rumor burden','Location','best','Box','off');
exportgraphics(fig,fullfile(outdir,'13_initial_condition_robustness.png'),'Resolution',240);

%% 16. Mesh-refinement study
meshN=[501 1001 2001 4001]; meshJ=zeros(size(meshN));
meshI=zeros(size(meshN)); meshA=zeros(size(meshN));
for j=1:numel(meshN)
    ts=linspace(0,T,meshN(j));
    [xs,~,us,~]=solve_sweep(ts,x0,p,w,umax,relax,tol,maxIter);
    ms=compute_metrics(ts,xs,us,p,w);
    meshJ(j)=ms.J; meshI(j)=ms.cumI; meshA(j)=ms.cumA;
end
errJ=abs(meshJ-meshJ(end))/meshJ(end);
errI=abs(meshI-meshI(end))/meshI(end);
errA=abs(meshA-meshA(end))/meshA(end);
meshH=T./(meshN-1);
fig=figure('Color','w','Position',[80 80 1100 430]);
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(meshN,meshJ,'o-','LineWidth',2,'Color',cc(1,:)); grid on
xlabel('grid points'); ylabel('objective'); title('Mesh convergence of objective');
nexttile; loglog(meshH(1:end-1),errJ(1:end-1),'o-','LineWidth',2,'Color',cc(3,:)); hold on
loglog(meshH(1:end-1),errI(1:end-1),'s-','LineWidth',2,'Color',cc(1,:));
loglog(meshH(1:end-1),errA(1:end-1),'d-','LineWidth',2,'Color',cc(4,:)); grid on
xlabel('step size h'); ylabel('relative error'); title('Relative error against N=4001');
legend('J','integral I','integral A','Location','best','Box','off');
exportgraphics(fig,fullfile(outdir,'14_mesh_refinement.png'),'Resolution',240);

%% 17. Save complete numerical workspace
save(fullfile(outdir,'full_optimal_control_results.mat'), ...
    'p','w','umax','x0','T','N','t','x','xNo','z','u','raw','H','history', ...
    'm','m0','checks','strategyNames','strategyBounds','strategyMetrics', ...
    'redI','redA','BIvals','BAvals','weightCumI','weightCumA','weightJ', ...
    'weightEffort','costScale','costI','costA','commonBound','boundI','boundA', ...
    'horizons','horizonI','horizonA','initialNames','initialStates', ...
    'initialRedI','initialRedA','meshN','meshH','meshJ','meshI','meshA', ...
    'errJ','errI','errA');

fprintf('\nAll 14 figures and the MAT data file were written to:\n%s\n',outdir);

%% LOCAL FUNCTIONS

function [x,z,u,history]=solve_sweep(t,x0,p,w,umax,relax,tol,maxIter)
%SOLVE_SWEEP Solve the coupled state-adjoint-projection fixed point.
% The fallback attempts reduce the new-candidate weight if the first fixed-
% point iteration has not reached the required tolerance.
attemptRelax=[relax 0.10 0.05];
attemptIter=[maxIter 2*maxIter 4*maxIter];
for attempt=1:numel(attemptRelax)
    rho=attemptRelax(attempt); limit=attemptIter(attempt);
    u=zeros(3,numel(t)); history=zeros(1,limit);
    for iteration=1:limit
        old=u;
        x=forward_states(t,x0,u,p);
        z=backward_adjoints(t,x,u,p,w);
        candidate=project_controls(unconstrained_controls(x,z,w),umax);
        u=rho*candidate+(1-rho)*old;
        history(iteration)=max(abs(u(:)-old(:)));
        if history(iteration)<tol, break; end
    end
    history=history(1:iteration);
    if history(end)<tol, break; end
end
x=forward_states(t,x0,u,p);
z=backward_adjoints(t,x,u,p,w);
end

function x=forward_states(t,x0,u,p)
%FORWARD_STATES Classical RK4 integration of the seven controlled states.
N=numel(t); h=t(2)-t(1); x=zeros(7,N); x(:,1)=x0;
for k=1:N-1
    um=0.5*(u(:,k)+u(:,k+1));
    k1=state_rhs(x(:,k),u(:,k),p);
    k2=state_rhs(x(:,k)+0.5*h*k1,um,p);
    k3=state_rhs(x(:,k)+0.5*h*k2,um,p);
    k4=state_rhs(x(:,k)+h*k3,u(:,k+1),p);
    x(:,k+1)=x(:,k)+h*(k1+2*k2+2*k3+k4)/6;
end
end

function z=backward_adjoints(t,x,u,p,w)
%BACKWARD_ADJOINTS Classical RK4, backward from z(T)=0.
N=numel(t); h=t(2)-t(1); z=zeros(7,N);
for k=N-1:-1:1
    xm=0.5*(x(:,k)+x(:,k+1)); um=0.5*(u(:,k)+u(:,k+1));
    k1=adjoint_rhs(x(:,k+1),z(:,k+1),u(:,k+1),p,w);
    k2=adjoint_rhs(xm,z(:,k+1)-0.5*h*k1,um,p,w);
    k3=adjoint_rhs(xm,z(:,k+1)-0.5*h*k2,um,p,w);
    k4=adjoint_rhs(x(:,k),z(:,k+1)-h*k3,u(:,k),p,w);
    z(:,k)=z(:,k+1)-h*(k1+2*k2+2*k3+k4)/6;
end
end

function dx=state_rhs(x,u,p)
%STATE_RHS Controlled coupled epidemic--rumor state equations.
S=x(1); I=x(2); V=x(3); R=x(4); %#ok<NASGU>
U=x(5); A=x(6); C=x(7);
u1=u(1); u2=u(2); u3=u(3);
q=p.psi*(1-p.alpha*A)+u1;
dx=[p.Lambda-p.beta*S*I-q*S+p.omega*V-p.mu*S;
    p.beta*S*I+p.sigma*p.beta*V*I-(p.gamma+u3+p.mu)*I;
    q*S-p.sigma*p.beta*V*I-(p.omega+p.mu)*V;
    (p.gamma+u3)*I-p.mu*R;
    p.eta-p.lam*U*A+p.theta*C-p.phi*I*U-p.mu*U;
    p.lam*U*A+p.phi*I*U-(p.delta+u2+p.mu)*A;
    (p.delta+u2)*A-(p.theta+p.mu)*C];
end

function dz=adjoint_rhs(x,z,u,p,w)
%ADJOINT_RHS The seven costate differential equations z'=-H_x.
% Rewritten without vertical concatenation of long expressions to avoid
% "dimensions not consistent" errors.
S=x(1); I=x(2); V=x(3); U=x(5); A=x(6);
zS=z(1); zI=z(2); zV=z(3); zR=z(4); zU=z(5); zA=z(6); zC=z(7);
u1=u(1); u2=u(2); u3=u(3);
BI=w(1); BA=w(2);
q=p.psi*(1-p.alpha*A)+u1;

dz1 = zS*(p.beta*I+q+p.mu) - zI*p.beta*I - zV*q;

dz2 = -BI + zS*p.beta*S ...
      - zI*(p.beta*S + p.sigma*p.beta*V - (p.gamma+u3+p.mu)) ...
      + zV*p.sigma*p.beta*V - zR*(p.gamma+u3) ...
      + (zU-zA)*p.phi*U;

dz3 = -zS*p.omega - zI*p.sigma*p.beta*I ...
      + zV*(p.sigma*p.beta*I + p.omega + p.mu);

dz4 = p.mu*zR;

dz5 = zU*(p.lam*A + p.phi*I + p.mu) - zA*(p.lam*A + p.phi*I);

dz6 = -BA - zS*p.alpha*p.psi*S + zV*p.alpha*p.psi*S ...
      + zU*p.lam*U - zA*(p.lam*U - (p.delta+u2+p.mu)) ...
      - zC*(p.delta+u2);

dz7 = -zU*p.theta + zC*(p.theta + p.mu);

dz = [dz1; dz2; dz3; dz4; dz5; dz6; dz7];
end

function raw=unconstrained_controls(x,z,w)
%UNCONSTRAINED_CONTROLS Stationary solutions of H_{u_j}=0.
raw=[(z(1,:)-z(3,:)).*x(1,:)/w(3);
     (z(6,:)-z(7,:)).*x(6,:)/w(4);
     (z(2,:)-z(4,:)).*x(2,:)/w(5)];
end

function u=project_controls(raw,umax)
%PROJECT_CONTROLS Euclidean projection onto [0,u1max]x[0,u2max]x[0,u3max].
u=max(raw,0);
u=min(u,repmat(umax,1,size(raw,2)));
end

function H=hamiltonian_path(x,z,u,p,w)
%HAMILTONIAN_PATH Evaluate H=L+z^T f at every time node.
N=size(x,2); H=zeros(1,N);
for k=1:N
    L=w(1)*x(2,k)+w(2)*x(6,k)+0.5*sum(w(3:5)'.*u(:,k).^2);
    H(k)=L+z(:,k)'*state_rhs(x(:,k),u(:,k),p);
end
end

function m=compute_metrics(t,x,u,p,w)
%COMPUTE_METRICS Objective, burdens, effort, peaks, and transition flows.
m=empty_metrics();
m.Jparts=[w(1)*trapz(t,x(2,:)), ...
          w(2)*trapz(t,x(6,:)), ...
          0.5*w(3)*trapz(t,u(1,:).^2), ...
          0.5*w(4)*trapz(t,u(2,:).^2), ...
          0.5*w(5)*trapz(t,u(3,:).^2)];
m.J=sum(m.Jparts);
m.cumI=trapz(t,x(2,:)); m.cumA=trapz(t,x(6,:));
m.peakI=max(x(2,:)); m.peakA=max(x(6,:)); m.finalState=x(:,end);
m.L1=[trapz(t,u(1,:)),trapz(t,u(2,:)),trapz(t,u(3,:))];
m.L2sq=[trapz(t,u(1,:).^2),trapz(t,u(2,:).^2),trapz(t,u(3,:).^2)];
m.maxControls=max(u,[],2)';
q0=p.psi*(1-p.alpha*x(6,:));
m.flows=[trapz(t,q0.*x(1,:)), ...
         trapz(t,u(1,:).*x(1,:)), ...
         trapz(t,u(2,:).*x(6,:)), ...
         trapz(t,u(3,:).*x(2,:))];
end

function m=empty_metrics()
%EMPTY_METRICS Fixed-layout structure used by preallocation and reporting.
m=struct('J',0,'cumI',0,'cumA',0,'peakI',0,'peakA',0, ...
    'finalState',zeros(7,1),'Jparts',zeros(1,5),'flows',zeros(1,4), ...
    'L1',zeros(1,3),'L2sq',zeros(1,3),'maxControls',zeros(1,3));
end

function plot_heatmap(ax,M,xvalues,yvalues,ttl)
%PLOT_HEATMAP Draw a labeled 3-by-3 numerical heatmap.
imagesc(ax,1:numel(xvalues),1:numel(yvalues),M); axis(ax,'xy');
colormap(ax,parula); colorbar(ax); title(ax,ttl);
set(ax,'XTick',1:numel(xvalues),'XTickLabel',xvalues, ...
       'YTick',1:numel(yvalues),'YTickLabel',yvalues);
xlabel(ax,'B_A'); ylabel(ax,'B_I');
span=max(M(:))-min(M(:));
for ii=1:size(M,1)
    for jj=1:size(M,2)
        if span==0 || M(ii,jj)<min(M(:))+0.55*span, color='k'; else, color='w'; end
        text(ax,jj,ii,sprintf('%.2f',M(ii,jj)), ...
            'HorizontalAlignment','center','Color',color,'FontSize',9);
    end
end
end