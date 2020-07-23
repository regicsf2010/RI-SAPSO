% clear, clc, close
function [err, sr, time, i] = main(idf, T, D, STOPC, GC)

%% Parameters
% function
% DIM = 10;
% FNAME = 'levi';
% RANGE = [-10, 10];
faux = str2func('cec17_func');
f = @(z) faux(z', idf);
g = load('g_opt.mat'); g_opt = g.g_opt;
RANGE = [-100, 100];

% PSO parameters
N = 20;
% MAXITER = 5000;
% STOPC = 1e-8;
k = 1;
CMAX = 3;
DT = [1e-2 .25];
dir = 1; % attractive as default
alpha = 3;

params = struct('IW', .7298, 'CC', 0, 'SC', 1.4962, 'GC', GC, 'd', D, 'a', alpha * pi / 180, 'i', eye(D));

%% Main iteration
L = norm(ones(1, D) * (RANGE(2) - RANGE(1))); % diagonal length of the space
VMAX = k * (RANGE(2) - RANGE(1)) / 2;

p = init(N, D, RANGE, VMAX, f);
G = getBestGlobal(p);

sr = 0;

tic
for i = 1 : T
    for j = 1 : N
        if(p.I(j) == 0)
            p.G(j, :) = getGradient(p.X(j, :), f); 
            p.G(j, :) = truncGrad(p.G(j, :), VMAX); % trunc gradient
        end
        
        p.V(j, :) = getVelocity(p.X(j, :), p.V(j, :), p.P(j, :), p.G(j, :), G, p.I(j), params, dir);
        p.V(j, :) = truncVel(p.V(j, :), VMAX); % trunc velocity
        
        p.X(j, :) = p.X(j, :) + p.V(j, :);
        [p.X(j, :), p.I(j), p.C(i)] = truncSpace(p.X(j, :), p.I(j), p.C(j), RANGE); % trunc space
        
        p.XFIT(j) = f(p.X(j, :));
        [p.P(j, :), p.PFIT(j), G] = updateBest(p.X(j, :), p.XFIT(j), p.P(j, :), p.PFIT(j), G);        
    end
    
    [p.I, p.C] = updateImportance(p.X, p.I, p.XFIT, p.OLDXFIT, p.C, p.G, G.X, CMAX, N);
    p.OLDXFIT = p.XFIT;
    
    diversity = getDiversity(p.X, L, N);
    [dir, p.I] = updateDir(dir, diversity, p.I, DT, N);
    
    if((G.XFIT - g_opt(idf)) <= STOPC)
        sr = 1;
        break
    end
end
time = toc;
err = G.XFIT - g_opt(idf);
end