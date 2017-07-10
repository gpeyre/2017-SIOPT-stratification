%%
% test for extended support of nuclear norm regularization.

addpath('toolbox/');

%%
% Setting.

n = 20;

rep = ['results/nuclear/'];
if not(exist(rep))
    mkdir(rep);
end

N = n*n;

% number of measurements
P = round(3/4*N);
% measurement matrix
Phi = randn(P,N);

%%
% helpers

mynorm = @(x)norm(x(:));
normnucl = @(x)sum(svd(x));
LowRank = @(r)randn(n,r)*randn(r,n);
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);
setfigname = @(name)set(gcf, 'Name', name, 'NumberTitle','off');

%%
% Compare CVX and FB.

% generate random low rank matrix
r = 1;
randn('state', 123);
x0 = LowRank(r);
y = Phi*x0(:);

lambda=0.1;

% solving with CVX
[xcvx,p,eta] = perform_nucl_reg_cvx(y,Phi,lambda);
% solving with BP
options.niter = 5000;
options.tau = 1.8/norm(Phi)^2;
options.report = @(x)0;
[xfb,Elist] = perform_nucl_reg_fb(y,Phi,lambda, options);
fprintf('Should be 0: %.3e\n', mynorm(xcvx-xfb)/mynorm(xcvx));

% norm(Phi*xfb(:)-y)/norm(y)
clf;
plot([svd(xcvx) svd(xfb)], '.-');
legend('CVX', 'BP');

%%
% Test for minimal norm certificate (aka limitting p when lambda->0)

r = 2;
x0 = LowRank(r);
y = Phi*x0(:);
[x,p,eta] = perform_nucl_reg_cvx(y,Phi,0);
% should be 0 for exact recovery
fprintf('Should be 0: %.3e\n', mynorm(x-x0)/mynorm(x0) );
[p,eta] = compute_certificate_nucl(x,Phi);

% compare support with extended support
clf;
plot([svd(x)/norm(x) svd(eta)], '.-');
axis tight;

%%
% Phase transition plots.

rlist = 1:8; % 8; % round(P/5);
rlist = 4;
nrep = 1e3; % number of replications
nrep = 50;

ExactRecovery = []; % recovery 0/1
RankExcess = []; % exceeding rank
RankFB = [];
X0 = [];
thresh_exactrecov = 1e-4;
thresh_rankexcess = 1e-4;
for i=1:length(rlist)
    progressbar(i,length(rlist));
    r = rlist(i);
    for j=1:nrep
        % generate random instance
        x0 = LowRank(r);
        y = Phi*x0(:);
        X0(:,:,i,j) = x0;
        % solve
        [x,p,eta] = perform_nucl_reg_cvx(y,Phi,0);
        ExactRecovery(i,j) = mynorm(x-x0)/mynorm(x0)<thresh_exactrecov;
        % compute exceeding rank
        [p,eta] = compute_certificate_nucl(x,Phi);
        RankExcess(i,j) = sum( abs(svd(eta)-1)<=thresh_rankexcess ) - r;
    end
end

save([rep 'data-svg'], 'X0', 'RankExcess', 'ExactRecovery');

ExactRecov = sum(ExactRecovery,2)/nrep;
SupportStab = sum(RankExcess==0,2)/nrep;
ermax = max(RankExcess(ExactRecovery==1));

%%
% Plot phase transition.

lw = 2;
ms = 15;
clf; hold on;
% legend
Ier = round(linspace(0,ermax,3));
lgd = {};
for er=Ier
    col = [er/ermax,0,1-er/ermax];
    plot(1,1,'color',col, 'LineWidth', lw);
    lgd{end+1} = ['\delta=' num2str(er)];
end
% plot
for er=0:ermax
    a = sum( (ExactRecovery==1) & (RankExcess<=er) , 2 )/nrep;
    col = [er/ermax,0,1-er/ermax];
    plot(rlist, a, '.-', 'LineWidth', lw, 'MarkerSize', ms, 'color', col);
end
axis tight; box on; drawnow;
legend(lgd)
SetAR(2/3);
set(gca, 'FontSize', 20);
saveas(gcf, [rep 'phase-transition-progress.eps'], 'epsc');


%%
% Histogram of low-complexity index.

% target low-complexity index
r = 4;

% compute histograms of identified manifold dimensions
S = RankExcess(rlist==r,:);
s = 0:max(S(:));
h = hist(S(:), s);
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

algo = 'fb';
algo = 'dr';

% FB param
options.niter = 2000;
options.tau = 1.8/norm(Phi)^2;
options.report = @(x)sum(svd(x)>1e-6);
lambda = 10; % for r=4
sigma = .1;
% DR param
options.gamma = 1;

% target excess
er_list = [0 3];
% compute FB solutions
SupportFB = {};
for i=1:length(er_list)
    er = er_list(i);
    % those signal with specific excess
    I = find( RankExcess(rlist==r,:)==er );
    SupportFB{i} = [];
    for j=I
        x0 = X0(:,:,rlist==r,j);
        w = sigma*randn(P,1);
        y = Phi*x0(:)+w;
        % perform FB/DR
        switch algo
            case 'fb'
                [xfb,Elist,G] = perform_nucl_reg_fb(y,Phi,lambda, options);
            case 'dr'
                [xfb,Elist,G] = perform_nucl_reg_dr(y,Phi,lambda, options);
        end

        SupportFB{i}(:,end+1) = G(:);
    end
    SupportFB{i} = max(SupportFB{i}, r+er);
end
% setup legend
clf; hold on;
lgd = {};
for i=1:length(er_list)
    t = (i-1)/(length(er_list)-1);
    col = [t 0 1-t];
    plot(1,1,'color', col, 'LineWidth', 2);
    lgd{end+1} = ['\delta=' num2str(er_list(i))];
end
% display
alpha = .08; % transparency for the plot
for i=1:length(er_list)
    t = (i-1)/(length(er_list)-1);
    col = [t 0 1-t];
    plot_semi_transparent(1:size(SupportFB{i},1),SupportFB{i},col,alpha);
    plot(1:size(SupportFB{i},1), mean(SupportFB{i},2), 'color', col, 'LineWidth', 2);
end
plot([1 options.niter], [1 1]*r, 'k:', 'LineWidth', 2);
legend(lgd);
axis([1 options.niter r-1 n]);
SetAR(2/3);
box on; drawnow;
set(gca, 'FontSize', 20);
saveas(gcf, [rep 'fb-iterations.eps'], 'epsc' );
saveas(gcf, [rep 'fb-iterations.png'], 'png' );

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
axis([1 options.niter r-1 n]);
SetAR(2/3);
box on; drawnow;
set(gca, 'FontSize', 20);
saveas(gcf, [rep 'fb-iterations-single.eps'], 'epsc' );
saveas(gcf, [rep 'fb-iterations-single.png'], 'png' );
