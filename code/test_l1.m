%%
% test for extended support of L1 norm regularization.

addpath('toolbox/');

%%
% Setting.

N = 100;

rep = ['results/l1/'];
if not(exist(rep))
    mkdir(rep);
end

% number of measurements
P = round(N/2);
% measurement matrix
Phi = randn(P,N);

%%
% helpers

mynorm = @(x)norm(x(:));
norml1 = @(x)sum(abs(x(:)));
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);
setfigname = @(name)set(gcf, 'Name', name, 'NumberTitle','off');

%%
% Compare CVX and FB.

% generate random sparse signal
r = 10;
randn('state', 123);
x0 = gen_sparse(N,r);
y = Phi*x0;

lambda=0.1;

% solving with CVX
[xcvx,p,eta] = perform_l1_reg_cvx(y,Phi,lambda);

% solving with BP
options.niter = 5000;
options.tau = 1.8/norm(Phi)^2;
[xfb,Elist] = perform_l1_reg_fb(y,Phi,lambda, options);
fprintf('Should be 0: %.3e\n', mynorm(xcvx-xfb)/mynorm(xcvx));

clf;
plot([xcvx xfb], '.-');
legend('CVX', 'BP');

%%
% Test for minimal norm certificate (aka limitting p when lambda->0)

r = 10;
x0 = gen_sparse(N,r);
y = Phi*x0(:);
[x,p,eta] = perform_l1_reg_cvx(y,Phi,0);
% should be 0 for exact recovery
fprintf('Should be 0: %.3e\n', mynorm(x-x0)/mynorm(x0) );
[p,eta] = compute_certificate_l1(x,Phi);

% compare support with extended support
clf;
plot([x/2 eta], '.-');
axis tight;

%%
% Phase transition plots.

rlist = 1:30;
rlist = 10;
nrep = 1e3; % number of replications
nrep = 150;

ExactRecovery = []; % recovery 0/1
SupportExcess = []; % exceeding support
SupportFB = [];
X0 = [];
thresh_exactrecov = 1e-4;
thresh_rankexcess = 1e-4;
for i=1:length(rlist)
    progressbar(i,length(rlist));
    r = rlist(i);
    for j=1:nrep
        % generate random instance
        x0 = gen_sparse(N,r);
        y = Phi*x0(:);
        X0(:,i,j) = x0;
        % solve
        [x,p,eta] = perform_l1_reg_cvx(y,Phi,0);
        ExactRecovery(i,j) = mynorm(x-x0)/mynorm(x0)<thresh_exactrecov;
        % compute exceeding sparsity
        [p,eta] = compute_certificate_l1(x,Phi);
        SupportExcess(i,j) = sum( abs(abs(eta)-1)<=thresh_rankexcess ) - r;
    end
end

% save for future use.
save([rep 'data-svg'], 'X0', 'SupportExcess', 'ExactRecovery');

lw = 2;
ms = 15;
ExactRecov = sum(ExactRecovery,2)/nrep;
SupportStab = sum(SupportExcess==0,2)/nrep;
ermax = max(SupportExcess(ExactRecovery==1));

%%
% Plot phase transition.

clf; hold on;
for er=0:ermax
    a = sum( (ExactRecovery==1) & (SupportExcess<=er) , 2 )/nrep;
    col = [er/ermax,0,1-er/ermax];
    plot(rlist, a, '.-', 'LineWidth', lw, 'MarkerSize', ms, 'color', col);
end
axis tight; box on; drawnow;
SetAR(2/3);
set(gca, 'FontSize', 20);
saveas(gcf, [rep 'phase-transition-progress.eps'], 'epsc');

%%
% Histogram of low-complexity index.

% target low-complexity index
r = 10;

% compute histograms of identified manifold dimensions
S = SupportExcess(rlist==r,:);
s = 0:min(20,max(S(:)));
h = hist(S(:), s); h(end)=0;
h = h/sum(h);
clf; hold on;
for i=1:length(s)
    t = (i-1)/(length(s)-1);
    col = [t 0 1-t];
    bar(s(i), h(i), 'Facecolor', col, 'EdgeColor', col);
end
set(gca, 'FontSize', 20);
SetAR(2/3);
box on;
axis tight;
saveas(gcf, [rep 'histogram-low-complexity.eps'], 'epsc');
saveas(gcf, [rep 'histogram-low-complexity.png'], 'png');

%%
% FB/DR tests

% target excess
er_list = [0 10];

algo = 'dr';
algo = 'fb';

% FB param
options.niter = 2000;
options.tau = 1.8/norm(Phi)^2;
options.report = @(x)sum(abs(x)>1e-6);
options.verb = 0;
lambda = .28; % for s=10
sigma = .1;
% DR param
options.gamma = 1;
if strcmp(algo, 'dr') %DR converges very fast
    options.niter = 600;
end

% compute FB/DR solutions
SupportFB = {};
for i=1:length(er_list)
    er = er_list(i);
    % those signal with specific excess
    I = find( SupportExcess(rlist==r,:)==er );
    SupportFB{i} = [];
    for j=I
        x0 = X0(:,rlist==r,j);
        w = 1e-3*randn(P,1);
        y = Phi*x0+w;
        % perform FB/DR
        switch algo
            case 'fb'
                [xfb,Elist,G] = perform_l1_reg_fb(y,Phi,lambda, options);
            case 'dr'
                [xfb,Elist,G] = perform_l1_reg_dr(y,Phi,lambda, options);
        end
        SupportFB{i}(:,end+1) = G(:);
    end
end
% setup legend
alpha = .15; % transparency for the plot
clf; hold on;
lgd = {};
for i=1:length(er_list)
    t = (i-1)/(length(er_list)-1);
    col = [t 0 1-t];
    plot(1,1,'color', col, 'LineWidth', 2);
    lgd{end+1} = ['\delta=' num2str(er_list(i))];
end
plot([1 options.niter], [1 1]*r, 'k:', 'LineWidth', 2);
% display
for i=1:length(er_list)
    er = er_list(i);
    t = (i-1)/(length(er_list)-1);
    col = [t 0 1-t];
    plot_semi_transparent(1:size(SupportFB{i},1),SupportFB{i},col,alpha);
    plot(1:size(SupportFB{i},1), mean(SupportFB{i},2), 'color', col, 'LineWidth', 2);
end
legend(lgd);
axis([1 options.niter 0 N]);
SetAR(2/3);
box on; drawnow;
set(gca, 'FontSize', 20);
addstr = '';
if strcmp(algo, 'dr')
    addstr = ['-gamma' num2str(round(1000*options.gamma))];
end
saveas(gcf, [rep algo '-iterations' addstr '.eps'], 'epsc' );
saveas(gcf, [rep algo '-iterations' addstr '.png'], 'png' );


%%
% Singes out samples.

lgd = {};
clf; hold on;
for i=1:length(er_list)
    t = (i-1)/(length(er_list)-1);
    col = [t 0 1-t];
    k = 1+floor(rand*size(SupportFB{i},2));
    plot(1:size(SupportFB{i},1),SupportFB{i}(:,k),'color',col, 'LineWidth', 2);
    lgd{end+1} = ['\delta=' num2str(er_list(i))];
end
plot([1 options.niter], [1 1]*r, 'k:', 'LineWidth', 2);
legend(lgd);
axis([1 options.niter r-1 N]);
SetAR(2/3);
box on; drawnow;
set(gca, 'FontSize', 20);
saveas(gcf, [rep algo '-iterations-single.eps'], 'epsc' );
saveas(gcf, [rep algo '-iterations-single.png'], 'png' );
