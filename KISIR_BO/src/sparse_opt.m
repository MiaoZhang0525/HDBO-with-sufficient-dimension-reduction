function model = sparse_opt(objective_func, num_iter, model,A,init_size,high_dim)

dopt.maxevals = 500;
dopt.maxits = 200;
dopt.showits = 0;
reinit=0;
r=0;
for i = 1:num_iter
    if reinit==20&&(i+init_size+r*init_size)<num_iter
        for j=1:init_size
            final_xatmin=rand(1, high_dim).*2-1;
            f_t = objective_func(final_xatmin);
            [A,model] = update_model(model, f_t, final_xatmin,A,init_size);
        end
        reinit=0;  
        r=r+1;
    else
        final_xatmin = maximize_acq(model,A, dopt, 'ucb');

        if model.high_dim > 0
            [f_t, record] = objective_func(final_xatmin');
            model.records(i+1, :) = record;
        else
            f_t = objective_func(final_xatmin');
            record = 0;
        end
        if f_t<max(model.f)
            reinit=reinit+1;
        end

        [A,model] = update_model(model, f_t, final_xatmin',A,init_size);
    end
    if num_iter-i-r*init_size==0
        break
    end
    
end


end