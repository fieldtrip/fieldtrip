function test_predictive_gamma_symm

X = create_mixed_source(8);

params.algorithm = 'symm';
params.gammaf.h = @gamma_predictive_symm;

state = denss(X, params);
