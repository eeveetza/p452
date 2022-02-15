success = 0;
fail = 0;
[s,f] = test_deltaBullington();
success = s;
fail = f;
%
[s,f] = test_path_1();
success = success + s;
fail = fail + f;
%
[s,f] = test_path_2();
success = success + s;
fail = fail + f;
%
[s,f] = test_ClutterLoss();
success = success + s;
fail = fail + f;

%%
fprintf(1,'*** %d (out of %d) tests succeeded.\n', success, success+fail);
fprintf(1,'*** %d (out of %d) tests failed.\n', fail, success + fail);
