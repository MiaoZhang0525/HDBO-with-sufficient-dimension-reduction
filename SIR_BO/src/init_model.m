function [A,model] = init_model(high_dim,dim, bounds, init_py,init_f, hyp, noise, ...
    cov_model)

A = KSIR(init_py, init_f, dim+1);



init_pt=init_py*A;

model.initsize =size(init_py,1);

model.use_direct = 0;
model.use_CMAES = 1;
model.use_fminsearch = 0;
model.use_gradient = 0;

model.noise = noise;
%%%%%%%%%%%%%%%%%x is d-dimension
if strcmp(cov_model,'se')
    model.cov_model = @(hyp, x, z, records)covSEiso(hyp, x, z);
elseif strcmp(cov_model,'ard')
    model.cov_model = @(hyp, x, z, i)covSEard(hyp, x, z);
end

model.kernel_type = cov_model;

model.hyp = hyp;
model.bounds = bounds;
model.d = dim;

model.cur_hyp = hyp(1);
model.exploit_count = 0;
model.hyper_bound = log([0.01, 50]);

model.prior_mean = 0;

model.max_expoit_count = 5;


model.records = 0;
model.high_dim = 0;

model.L = (model.cov_model(model.hyp, init_pt, init_pt) + ...
    model.noise*eye(size(init_pt,1)));
model.L = chol(model.L, 'lower');

model.X = zeros(3000, high_dim);
model.XT = zeros(3000, dim);
model.X(1:size(init_py,1), :) = init_py;
model.XT(1:size(init_py,1), :) = init_pt;


model.f = init_f;

model.m = size(init_py,1);
model.n = size(init_py,1);


[a,b]=max(init_f);

model.max_val = a;
model.max_x = model.X(b, :);
model.max_xT = model.XT(b, :);
model.display = 1;


model.copts.MaxFunEvals = 2000; model.copts.TolFun = 1e-8; 
model.copts.LBounds = model.bounds(:, 1);
model.copts.UBounds = model.bounds(:, 2); 
model.copts.SaveVariables = 0; model.copts.LogModulo = 0;
model.copts.DispFinal = 0; model.copts.DispModulo = 0;


model.last_display = 0;