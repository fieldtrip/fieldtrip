function sendToFieldtripDashboard(obj)
% send test results to FieldTrip dashboard
%
% Notes:
%   - unlike ft_run_test, this function sends all results after all the
%     tests have been run

    url = 'http://dashboard.fieldtriptoolbox.org/api/';

    n_outcomes=countTestOutcomes(obj);

    verbosity=getVerbosity(obj);
    output_stream=getStream(obj);

    % helper function to report progress
    report=@(msg)report_with_verbosity(output_stream,verbosity,msg);
    report(sprintf('Sending %d test results to %s\n',n_outcomes,url));

    n_sent=0;
    for k=1:n_outcomes
        test_outcome=getTestOutcome(obj,k);
        test_case=getTest(test_outcome);

        passed=isNonFailure(test_outcome); % passed or skipped
        skipped=passed && ~isSuccess(test_outcome); % skipped only

        functionname=getName(test_case);
        location=getLocation(test_case);

        if skipped
            % TODO decide how to deal with skipped tests.
            % For now the skipped tests are not reported at all to the
            % dashboard.
            if verbosity>=2
                msg=sprintf('Skipped %03d / %03d: %s (%s)\n',...
                            k,n_outcomes,functionname,location);
            elseif verbosity>=1
                msg='s';
            else
                msg='';
            end

            report(msg);

            continue
        end

        runtime=getDuration(test_outcome);

        revision=ft_version();

        result=struct();
        result.matlabversion    = version('-release');
        result.fieldtripversion = revision;
        result.branch           = ft_version('branch');
        result.arch             = computer('arch');
        result.hostname         = run_private_fieldtrip('gethostname');
        result.user             = run_private_fieldtrip('getusername');
        result.passed           = passed;
        result.runtime          = runtime;
        result.functionname     = functionname;

        % send results
        moxunit_fieldtrip_util_send_json(url, result);

        outcome_str=getOutcomeStr(test_outcome,verbosity);
        if verbosity>=2
            msg=sprintf('Sent %03d / %03d [%s]: %s (%s)\n',...
                            k,n_outcomes,...
                            outcome_str,functionname,location);
        elseif verbosity>=1
            msg=outcome_str;
        else
            msg='';
        end
        report(msg);
        n_sent=n_sent+1;
    end

    report(sprintf('Completed sending %d / %d test results to %s\n',...
                                    n_sent,n_outcomes,url));

function report_with_verbosity(stream,verbosity,msg)
    if ~isempty(msg)
        if verbosity>0
            fprintf(stream,'%s',msg);
        end
    end


function result=run_private_fieldtrip(func_name)
% run function in ${FIELDTRIP_ROOT}/utilities/private
%
% somewhat ugly, we have to 'cd' to the directory in order to access the
% private function.
    orig_pwd=pwd();
    pwd_cleaner=onCleanup(@()cd(orig_pwd));

    fieldtrip_root_dir=fileparts(which('ft_defaults'));
    fieldtrip_util_private_dir=fullfile(fieldtrip_root_dir,...
                                        'utilities','private');

    cd(fieldtrip_util_private_dir);
    func=str2func(func_name);
    result=func();
