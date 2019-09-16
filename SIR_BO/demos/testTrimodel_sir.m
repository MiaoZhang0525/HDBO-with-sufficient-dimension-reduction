%function testAckley(varargin)
clear;
warning off;
store_data=[];
num_expe=10;


for ne=1:num_expe

    total_iter = 500;
    high_dim = 20000;
    dim = 20;

    init_size=100;
    %% Random embedding with true intrisic dimensions
    %  True maximum is ensured to fall in to the bounds

    model = rembo(total_iter-init_size, dim, high_dim,init_size);
    %ditance_log_plot(total_iter, model.f);
    fvalue=model.f;
    fvalue=fvalue(1:500,:);
    store_data=[store_data,fvalue];
    clear model;
end
save SIR_BO-robust-20000-20p-Trimodel.mat store_data

%% Run rembo.
function model = rembo(total_iter, dim, high_dim,init_size)
    % total_iter: total number of iterations.
    % dim: embedding dimension.
    % high_dim: ambient dimension.

   %     A = randn(high_dim, dim);           
 %       A = eye(high_dim, dim);

%        scale = max(1.5*log(dim), 1);
    bounds = stardardBounds(high_dim);                 % Initialize bounds.
    obj_fct = @(x) multi_model(x);     % Initialize the objective function.

    init_py = rand(init_size, high_dim).*2-1;     % Initial point.

    for i=1:size(init_py,1)
        init_f(i,1) = obj_fct(init_py(i,:));  
    end% Evaluate initial point.

    hyp = [ones(dim, 1)*0.1 ; 1];          % Setup initial hyper-parameters.
    hyp = log(hyp);
    % Initialize model.
    [A,model] = init_model(high_dim,dim, bounds, init_py,init_f, hyp, 1e-10, 'ard');  
    % Do optimization.
    model = sparse_opt(obj_fct, total_iter-1, model,A,init_size,high_dim);
end

