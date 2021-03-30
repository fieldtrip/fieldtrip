classdef duneuro_function < handle
    properties (Hidden = true)
        cpp_handle;
        driver;
    end
    methods
        % Constructor
        function this = duneuro_function(driver)
            this.driver = driver;
            this.cpp_handle = duneuro_matlab('make_domain_function', driver.cpp_handle);
        end
        % Destructor
        function delete(this)
            duneuro_matlab('delete_function', this.cpp_handle);
        end
    end
end
