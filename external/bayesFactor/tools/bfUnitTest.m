classdef bfUnitTest < matlab.unittest.TestCase
    % Class to test the Bayes Factor toolbox using the Matlab Unit Test
    % framework.
    % BK - Dec 2019
    
    methods (TestClassSetup)
        function importPackages(testCase) %#ok<MANU>
            import matlab.unittest.*
            import matlab.unittest.constraints.*            
        end
    end
    
    methods (Test)
        function testFromF(testCase)
            actual = bf.bfFromF(10,1,3,5);
            expected = 17.4812;
            testCase.verifyThat(actual,matlab.unittest.constraints.IsEqualTo(expected,'Within',matlab.unittest.constraints.AbsoluteTolerance(0.0001)))
        end
        
        function testFromT(testCase)
            actual = bf.bfFromT(2,5);
            expected = 2.9574;
            testCase.verifyThat(actual,matlab.unittest.constraints.IsEqualTo(expected,'Within',matlab.unittest.constraints.AbsoluteTolerance(0.0001)))
        end

    end
end