function test_basic_defl

X = create_mixed_source(4);

params.algorithm = 'defl';

state = denss(X, params);
