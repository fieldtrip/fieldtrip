function test_example_glm_trials

% MEM 4gb
% WALLTIME 00:10:00

%
%% Using General Linear Modeling over trials
%
% Most MEG/EEG experimental paradigms consists of a repeated task or stimulus presentation, during which data is acquired. Subsequently the data is segmented (based on the events), pre-processed, and for example averaged to construct an ERP. These ERPs can be computed separately for different conditions, and contrasted over conditions to see when and where the cortcial processing differs between experimental conditions. In these trial-based analysis it is common to treat ever channel and timepoint asround the stimulus separate from all others, and do the averaging or statistical analysis for each channel-time point. Cluster-based permutation tests can be used to deal with the multiple comparisons, but the clusters themselves do not form an explicit model of cortical activity.
%
% An alternative to modeling all channel-time points over trials independently from each other is to make explicit models for the precise temporal structure in the data around each event (i.e. within the trial). This is common for event-related fMRI and fNIRS, where the cognitive activity related to events is very short compared to the low temporal resolution of the data. As the precise millisecond timing of the cortical activity is not represented in the data, the stimuli can be modeled as stick-functions and convolved with the BOLD haemodynamic response function. In MEG and EEG it is also possible to make explicit models including peri-stimulus time and use GLMs to estimate the event-related activity (REFS).
%
% The example page on [Using General Linear Modeling on time series data](/example/glm_timeseries) shows how to use GLMs to model the event-related time series. This specific page will show how to use GLMs to model the trial-based structure in the data.
%
%% # GLMs for two conditions
%
% When considering the contrast between two experimental conditions, we can describe this as an ANOVA with one factor and two levels. In this simple experimental design, we do not have to use ANOVAs and F-tests, but we can resort to t-tests.
%
% For example, in the [eeg-language](/tutorial/eeg_language) dataset, there are trials with a visual stimulus presentation and trials with an auditory stimulus presentation. The procedure for the analysis looks like this:
%
%*   preprocess the data in both conditions
%*   average per condition and compute condition difference
%*   compute contrast with t-test
%
%% # GLMs for three conditions
%
% Considering an experimental manipulation with three conditions, we can describe this as an ANOVA with one factor and three levels. We could use pair-wise t-tests, which puts us back in the previous situation. We need an F-test for testing an overall effect of the three conditions (or levels).
%
%*   preprocess the data in all three conditions
%*   compute an F-test
%
%% # GLMs for a multi-factorial design
%
% Considering an experiment from the general perspective of an ANOVA, we can have multiple factors and/or multiple levels per factor.
%
% In the [eeg-language](/tutorial/eeg_language) dataset there is one factor for "stimulus modality" with three levels for pictures, visual, and auditory. Another factor is the "stimulus category" with two levels for animals and tools. This can be represented in a 3-by-2 table.
%
% Note that the data in the [eeg-language](/tutorial/eeg_language) dataset has additional characteristics that are experimentally controlled, such as different stimulus items in each category, targets versus non-targets, and also responses on the targets with varying response latencies. We are not further considering these here, but you could extend the table with these factors.
%
% The procedure for the analysis looks like this:
%
%*   preprocess the data for all 3x2 conditions combined
%*   average per condition and compute pairwise condition differences
%*   estimate GLMs
%*   test different contrasts (t, F)
%
%% # Within-subject versus group analysis
%
% As demonstrated above, GLMs can be used for within-subject analysis where the variability over trials is used in combination with the model estimates to determine significance. For group analysis it is common to make 1-st level model estimates for each subject, and do a 2-nd level analysis on the group level.
%
% The procedure for the analysis looks like this:
%
%*   preprocess the data in all 3x2 conditions for each subject
%*   1st level: estimate GLMs for each subject
%*   2nd level: do statistics across subjects for different contrasts (t, F)
