function [acq, grad] = acq(model, xx,A, type)

si = 1;

if nargin < 4
	type = 'ucb';
end

%%%%%%%%%%%%%x translate based on A original input x'
x=(xx*A);
[mu, var] = mean_var(model, x);
sigma = sqrt(var);

% type
if strcmp('ucb', type)
    coeff = si*sqrt(2*log(model.n^(model.d/2+2) *pi^2 /(3*0.1)));
	acq = -(mu+coeff*sigma);
	
elseif strcmp('ei', type)
	% diff = (mu - model.max_val-0.01);
	% Z = diff./sigma;
 	% [npdf, ncdf] = stand_normal_stats(Z);
	% acq = -(ncdf.*diff + npdf.*sigma);

	% if compute_gard
	% 	Z_grad = (mu_grad*sigma - sigma_grad*diff)/var;
	
	% 	grad = ncdf*mu_grad + diff*npdf*Z_grad;
	% 	grad = grad + npdf*sigma_grad + sigma*npdf*(-Z)*Z_grad;
	% 	grad = -grad/model.noise;
	% end

	acq = -log_exp_imp(-model.max_val, -(mu-0.0001), sigma);
end


function [npdf, ncdf] = stand_normal_stats(Z)
	npdf = (2*pi)^(-0.5)*exp(-(Z.^2)/2);
	ncdf = 0.5*(1+erf(Z/sqrt(2)));
end

end