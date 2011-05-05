addpath(pwd);
tests = dir('tests');
cd('tests');

num_tests = 0;
num_failed = 0;
for i=1:length(tests)
  test = tests(i).name;
  % detect tests b/c directory contains other stuff (e.g., surdirs and
  % helper files)
  if ( (~tests(i).isdir) & strncmp(test, 'test', 4) )
    tic
    f = str2func(test(1:end-2));
    num_tests = num_tests + 1;
    disp(['** Running test(s) in: ' test ]);
    [pass,str] = f();
    if pass
      disp(['** PASS: ' str]);
    else
      disp(['** FAIL: ' str]);
      num_failed = num_failed + 1;
    end
    toc
  end
end

if (num_failed > 0)
  disp(['***** WARNING: ' num2str(num_failed) ...
        ' of ' num2str(num_tests) ' tests failed *****']);
else
  disp(['***** All ' num2str(num_tests) ' tests passed *****']);
end
cd('..')