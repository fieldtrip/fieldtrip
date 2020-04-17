Rename after an issue number/pull request has been set.

Test script will be updated here soon.
@mcpiastra and I will prepare a test script.


The current version of ft_sourcenanalysis recomputes the leadfields (lines 373 to 403). 
It ignores the fact that the input may already contain pre-computed leadfields. So for any leadfield computed in a non-standard way (for example FEM computed externally, but the issue is not unique to FEM).
Suggested changes would best be done in "prepare_headmodel", which could have everything needed set, before it goes on with the inverse solution (suggestion by @robertoostenveld).
