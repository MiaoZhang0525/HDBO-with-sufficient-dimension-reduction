function y = Trimodel(x)
  %  x=x(1:10);
    probs = [0.1; 0.8; 0.1];
    numDims=2;
    gaussVar = 0.01 * numDims^0.1;
    covarDiag = gaussVar * ones(1, numDims);
    centres12 = [0.62 -0.38; 0.18 0.58; -0.58 -0.56];
    centresRest = [-0.66 -0.19 0.62]';
    centres = [centres12, repmat(centresRest, 1, numDims -2)];
    y=log(probs(1) * mvnpdf(x(1:2), centres(1, :), covarDiag ) +...
            probs(2) * mvnpdf(x(1:2), centres(2, :), covarDiag ) +...
            probs(3) * mvnpdf(x(1:2), centres(3, :), covarDiag ) );
end