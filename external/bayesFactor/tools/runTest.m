% Script to run the unit tests that test Bayes Factor Toolbox
% functionality.
import matlab.unittest.TestCase
import matlab.unittest.TestSuite
import matlab.unittest.constraints.IsEqualTo
import matlab.unittest.constraints.AbsoluteTolerance
import matlab.unittest.constraints.RelativeTolerance


suiteClass = TestSuite.fromClass(?bfUnitTest);
result = run(suiteClass);
table(result)