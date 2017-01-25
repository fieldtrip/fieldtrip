`MOxUnit_fieldtrip` adds unit test functionality to [FieldTrip] using [MOxUnit].
Contributor: Nikolaas N. Oosterhof <nikolaas.oosterhof@unitn.it>

### Features
- Run tests using Matlab or GNU Octave.
- Tests can be filtered using 'maxwalltime', 'maxmem', or dependencies.
- Results can be uploaded to fieldtrip's dashboard
- Tests can be run using Travis-ci

### Implementation notes
`MOxUnit_fieldtrip` adapts the testing functionality from MOxUnit by subclassing the following classes:
- MOxUnitFieldTripTestCase implements a test case by subclassing MOxUnitFunctionHandleTestCase. For a FieldTrip test, it defines a function handle that, when evaluated, `cd`s into the directory of the test and then runs it. After running the test it tries to close any figure windows that were openened by the test. 
  When this test is instantiated, it reads the test .m file and looks for lines that define 'WALLTIME', 'MEM' or 'TEST' dependencies. Any of those values are stored internally and can be retrieved using the 'get' method.
- MOxUnitFieldTripTestSuite implements a test suite by subclassing MOxUnitTestSuite. It overrides addFromFile method by using MOxUnitFieldTripTestCase's constructor. The filter method can be applied to exclude tests that exceed a particular 'WALLTIME', 'MEM' or 'TEST' dependency requirement; excluded tests are replaced by a test case that raises a skipped-test exception that is caught by the test framework. 



[FieldTrip]: http://www.fieldtriptoolbox.org
[MOxUnit]: https://github.com/MOxUnit/MOxUnit
[Travis-ci]: https://travis-ci.org



