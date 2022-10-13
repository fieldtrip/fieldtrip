# Parallel computing branch of FieldTrip

This is a experimental branch of the FieldTrip code base that resulted from a collaborative project with Vijay Iyer at [MathWorks](https://www.mathworks.com) and Aljen Uitbeijerse at [VORtech](https://www.vortech.nl). The goal of the project was to explore whether "under the hood" parallelization _inside_ the FieldTrip functions, just above the level of the actual algorithms, would reduce the execution time and the time that researchers need for their analyses.

The project design is documented [here](https://www.fieldtriptoolbox.org/development/project/parallel/) on the FieldTrip website. The results of the project are organized and discussed on the following places:

-   GitHub project [5](https://github.com/fieldtrip/fieldtrip/projects/5)
-   GitHub issue [1851](https://github.com/fieldtrip/fieldtrip/issues/1851), [1852](https://github.com/fieldtrip/fieldtrip/issues/1852), [1853](https://github.com/fieldtrip/fieldtrip/issues/1853), [2068](https://github.com/fieldtrip/fieldtrip/pull/2068), [2069](https://github.com/fieldtrip/fieldtrip/pull/2069)
-   the repositories <https://github.com/apjanke/fieldtrip-parallel-support> and <https://github.com/AljenU/fieldtrip-parallel-support>

The project outcome was that - although there were some speedups - we decided not to merge this in the release branch and not continue with a more wide-spread implementation of this "under the hood" parallelization. The impact on the flow of the code is too significant, and interferes with efficient memory handling.

## What was changed

Parallel constructs were added in ft_timelockanalysis, ft_networkanalysis and ft_freqanalysis, as recorded in [this branch](https://github.com/fieldtrip/fieldtrip/tree/parallel). An earlier effort had identified ft_timelockanalysis, ft_networkanalysis and ft_freqanalysis as functions that might benefit from adding parallel constructs. Thus, to evaluate application of parallel constructs in high-level Fieldtrip analysis functions, these three functions were adapted. In particular, these functions have a for loop as part of their analysis, and changing those to a parfor loop should easily integrate parallelization in these functions.

## How it was changed

In each of the three functions, the for loop was replaced with a parfor loop. Using a parfor loop introduces constraints on how data can be used. The consequences will be detailed here, with reference to an example code change. When referring to examples, this is to the implementation on the 'parallel' branch on the fieldtrip github repository.

Using a parfor loop requires that the interpreter can determine that each loop iteration can be executed independent of the other iterations. Thus there are restrictions on how data is used inside the loop. These are well documented in the Matlab documentation. It mainly affects variables that are both written to inside the loop, and also used outside the loop.

1.  If data is assigned to a subset of a variable, it must be done in exactly the same way in each iteration of the loop. In current fieldtrip implementations, often a switch-case or if-else construct is used to select different data assignments, also within loops. However, depending on the exact use-case, it is possible to refactor the fieldtrip implementations to have similar assignments within each branch within the loop, and still have a different data layout and usage for further processing. The changes with regard to the variable `output` in ft_networkanalysis show an example of how this can be done, without increasing the memory footprint of the variable. In this case by adding extra dimensions of size 1 to the data, and using a (cheap) reshape after the loop to remove those extra dimensions.

2.  A variable cannot be both a reduction variable and a subset-assigned variable within the loop. This is a variant of the restriction that assignments need to be similar in each loop iteration. This situation is present in ft_timelockanalysis for the variable covsig, and in this case there were two ways to deal with it. The first way to deal with it is to use only the subset-way of assigning to the variable, and do the reduction after the loop. However, this would increase the memory footprint, since all intermediate results would be kept until the loop is finished. And the reduction appears to be present exactly to reduce the memory footprint. The second way to deal with the non-similar assignment, is to reverse the order of entering-the-loop and switching-for-different-data-assignment. Then it is possible to keep the memory footprint of the variable the same. The implementation in ft_timelockanalysis shows an example of this. There, each branch of the `if keeptrials` switch has its own parfor loop (instead of a single loop with an if statement inside it).

3.  The limits of for loops nested in the parfor-loop must be specific. "You must define the range of a for-loop nested in a parfor-loop by constant numbers or broadcast variables." The change in ft_networkanalysis for the upper limit of the inner loop shows an example of how this can be done.

For more efficient memory use, there are additional recommendations, also available through the Matlab documentation. These are on variables that are used, but not changed, inside the loop. If only a subset of the data in the variable is needed in an iteration, the interpreter should be enabled to determine that indeed only a subset is needed on the parallel worker.

1.  Use only the needed subset of struct type variables. The full struct will be sent to each parallel worker when, within the loop, any fields are referenced through dot notation within the struct variable. If only a limited number of fields is needed within the loop, it is possible to limit how much of a struct type variable is sent to a worker. This is done by adding a variable that only holds the value of the needed field, and use that new variable inside the loop. The use of the variable datacov_trial in ft_timelockanalysis is an example. Note that, due to Matlab's copy-on-write, the value will NOT be duplicated in the main thread.

2.  Use the recommended, as per the Matlab documentation, way of indexing into variables where only a subset of the data is used in a particular iteration. If the parfor loop index is used for the subset identification, only the subset used in the iteration will be sent to the parallel worker. The use of `data_trial` in ft_freqanalysis is an example.

To give the user control over the amount of parallelization, including not doing the loops in parallel, an extra input option was added to the three functions. A cfg option was introduced to control the amount of workers used, `cfg.parformaxworkers = ft_getopt(cfg, 'parformaxworkers' , Inf)`. Setting cfg.parformaxworkers to 0, the parfor loops will not use parallel computing. The default setting of Inf will use all available workers in the running parallel pool.

## Different amounts of adaptation were needed

In the previous section, the examples on what was changed have referred mainly to ft_timelockanalysis and ft_networkanalysis. This was done, because the amount of changes between the parallel branch and the parent are very limited for these functions. Thus, it is more easy to see the individual restrictions and recommendations being applied.

In ft_freqanalysis, many more changes were made. The function had a for loop that assigns to multiple variables, and does the assignment differently for each of the variables, using the common fieldtrip pattern of having logical branching within the loop. The git commit history of ft_freqanalysis shows the various steps that were taken, to finally conform with the restrictions for using parfor. Together they show the effort needed to limit the memory footprint of variables, while implementing parallel constructs inside a fieldtrip function that uses multiple non-uniform data assignments.

A pattern that needed to be replaced specifically in ft_freqanalysis, was the allocation of memory for loop-output variables within the first iteration of the loop. To be able to use a parfor loop, the memory for the loop-output variables needs to be assigned before the loop is started. Only then can the first restriction mentioned earlier be adhered to. In order to do so without extensive code duplication, the memory-allocation part of the loop, and the before-the-memory-allocation part of the loop, were each put into their own subfunction. That way, the code duplication is limited to a single line: the call of the new before-the-memory-allocation subfunction, which is called both one time before the loop, and once in each iteration of the loop.

## Measured effects on runtime

The effects of the changes on runtime are different for the three functions, and (as expected) also depend on the specific input to a function. The effects were measured on a computer with an Intel Core i7-3770 and 32 Gb memory. The detailed timing results are available [here](https://github.com/fieldtrip/fieldtrip/issues/1853#issuecomment-1261971975). The results compare three cases: the original code, the parfored code version when run in parallel, and the parfored code version when run with parallelization disabled.

-   For ft_timelockanalysis changing the main loop from for to parfor has only small positive effects, if any, with the tested inputs. The small positive effects (a runtime decrease up to about 15%) do happen mainly for the longer runs, and when the loop is indeed run in parallel. This suggests that the change might be beneficial in actual workflows, where even more / larger data is processed.

-   For ft_networkanalysis changing the main loop from for to parfor has only negative effects, if any, with the tested inputs. The negative effects occur when running in parallel, where the runtime can be double the original. When parallelization is disabled, the parfored code version has the same runtime as the original.

-   For ft_freqanalysis changing the main loop from for to parfor has generally positive effects. When run in parallel, the runtime is almost halved in quite some cases, and decreased in nearly all. The only cases where the runtime does not decrease, are the ones that have a very short runtime to start with. When parallelization is disabled, the parfored code version has approximately the same runtime as the original, though some slight increases in runtime are observed.

## Side-effects

Running an analysis in parallel will increase the total amount of memory needed. Even though the implementations did take into account the per-variable memory footprint, dividing the calculations over multiple workers will increase the total memory footprint. Because functionality and data that is needed in each iteration, also from Matlab itself, will be duplicated. Thus, if an analysis is already memory-limited when run sequentially, running it in parallel on the same computer can cause memory-overflow, and actually greatly increase the runtime of the calculation.

Using a parfor loop where each iteration does very little work can increase the runtime. The overhead of sending data to and receiving data from the parallel workers, will in that case take more time than what is saved by dividing the total number of calculations over multiple workers.

# Conclusions

Changing for loops to parfor loops can reduce the runtime of certain high-level analyses, see for example the measured effects on ft_freqanalysis. In particular when an analysis takes more than a few seconds to run, and extra computer memory is available, enabling parallelization can speed up the workflow.

Whether a specific user will see a speedup for his workflow, when parfors are introduced, will at least depend on the specifications of the computer the analysis is run on. It might improve significantly, it could be similar, but it can also be much worse due to memory overflow. Thus, if parfors are introduced, it might be prudent to make the default number of workers equal to 0, and require the user to explicitly enable parallelization throught the cfg option.

Whether a specific user will see a speedup for his workflow, when parfors are introduced, will at least depend on the analyses that are part of the workflow. Certain analyses are more likely to benefit from parallelization than others. Therefore, deciding which analysis should have parallelization introduced should take into account the typical workflow, including typical data sizes, that include such an analysis.

Looking at the way many high-level analysis functions are now structured, introducing parallelization requires significant code changes. High-level analysis functions currently can have many lines of code per function, and extensive logical branches within the main calculation part of the funcion, including within loops. To avoid code duplication and / or to adhere to the parfor restrictions, such functions need significant changes to their structure. Most of these structure changes are actually also recommendable from a general software development perspective, but will need to be weighed against other aims the current structures may have.
