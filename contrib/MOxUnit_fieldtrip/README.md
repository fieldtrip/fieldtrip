`MOxUnit_fieldtrip` adds unit test functionality to [FieldTrip] using [MOxUnit].
Contributor: Nikolaas N. Oosterhof <nikolaas.oosterhof@unitn.it>

### Features
- Run tests using Matlab or GNU Octave.
- Based on [MOxUnit].
- Tests can be filtered using 'maxwalltime', 'maxmem', or dependencies.
- Results can be uploaded to fieldtrip's dashboard.
- Tests can be run on open-source continuous integration platforms with Travis-ci and shippable.com using GNU Octave.

### Implementation notes
`MOxUnit_fieldtrip` adapts the testing functionality from MOxUnit by subclassing the following classes:

- `MOxUnitFieldTripTestCase` implements a test case by subclassing `MOxUnitFunctionHandleTestCase`. For a FieldTrip test, it defines a function handle that, when evaluated, `cd`s into the directory of the test and then runs it. After running the test it tries to close any figure windows that were openened by the test.
  When this test is instantiated, it reads the test .m file and looks for lines that define 'WALLTIME', 'MEM' or 'TEST' dependencies. Any of those values are stored internally and can be retrieved using the 'get' method.
- `MOxUnitFieldTripTestSuite` implements a test suite by subclassing MOxUnitTestSuite. It overrides addFromFile method by using MOxUnitFieldTripTestCase's constructor. The filter method can be applied to exclude tests that exceed a particular 'WALLTIME', 'MEM' or 'TEST' dependency requirement; excluded tests are replaced by a test case that raises a skipped-test exception that is caught by the test framework.
- `MOxUnitFieldTripTestReport` manages test reports by subclassing `MOxUnitTestReport`. It implements a `sendToFieldtripDashboard` function that sends test results to the Fieldtrip dashboard.

`moxunit_fieldtrip_runtests` is the main function to run tests. It supports the same arguments as `ft_test run`, including `maxwalltime`, `maxmem` `upload` and `dependency`. To run tests through `MoxUnit_fieldtrip` one can also use `ft_moxunit_test run`.

### Continuous integration (CI) testing
Currently CI testing is supported through the travis-ci and shippable.com services. The configuration is set in `${FIELDTRIP_ROOT}/.travis.yml`.

### How to write tests
- Test files should be in the `${FIELDTRIP_ROOT}/test` directory and should start with 'test' or 'failed' in order to be included in the test suite.
- Test files should consist of a function that, when called, runs the test. If a test cannot be run it is considered to have errored. If it can be run but raised an exception, the test runner considers it as failed or skipped (depending on the exception; see below). A test that runs without throwing an exception is considered as passed.
- Tests can be skipped by calling `moxunit_throw_test_skipped_exception`. This throws an exception which is caught by the test runner and reported afterwards as a skipped test.
- The comment section can contain lines indicating maximume execution time (e.g. WALLTIME 00:10:00 for 10 minutes) or maximum memory (e.g. MEM 1gb for one gigabyte of memory).

### License and authors
For `MOxUnit_fieldtrip`'s copyright information and license terms, see the COPYING file distributed with `MOxUnit_fieldtrip`.


[FieldTrip]: http://www.fieldtriptoolbox.org
[MOxUnit]: https://github.com/MOxUnit/MOxUnit
[Travis-ci]: https://travis-ci.org
[Shippable]: https://www.shippable.com



