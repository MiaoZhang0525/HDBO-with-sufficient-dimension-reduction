function [mean, var] = mean_var(model, x)


k_tt = model.cov_model(model.hyp, x, x);
k_x = model.cov_model(model.hyp, model.XT(1:model.n,:), x);

intermediate = model.L'\(model.L\k_x);
mean = model.prior_mean+intermediate'*(model.f-model.prior_mean);
var = diag(k_tt - k_x'*intermediate);