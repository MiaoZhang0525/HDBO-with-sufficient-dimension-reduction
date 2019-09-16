%function testAckley(varargin)
clc;
clear;
warning off;
store_data=[];
num_expe=10;
high_dim = 200;
dim = 10;
init_size=250;
total_iter = 500;
total_iter = total_iter -init_size;
for ne=1:num_expe
    model = rembo(total_iter, dim, high_dim,init_size);
    %ditance_log_plot(total_iter, model.f);
    fvalue=model.f;   
    store_data=[store_data,fvalue];
    clear model;
end
save KISIR_BO-200-50p-Trimodel.mat store_data

%% Run rembo.
function model = rembo(iter, dim, high_dim,init_size)

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
    model = sparse_opt(obj_fct, iter, model,A,init_size,high_dim);
end
