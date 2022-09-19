# Contributing to FieldTrip

There are several ways in which you can contribute to the ongoing development of FieldTrip:

## Bug reporting

1. Please search on both the [FieldTrip discussion list](http://www.fieldtriptoolbox.org/discussion_list)
   and the [GitHub issue list](https://github.com/fieldtrip/fieldtrip/issues)
   to see if anybody else has lodged a similar observation.

2. How confident are you that the behaviour you have observed is in fact a
   genuine bug, and not a misunderstanding?

   -  *Confident*: Please [open a new GitHub issue](https://github.com/fieldtrip/fieldtrip/issues/new);
      select the "bug report" issue template to get started.

   -  *Not so confident*: That's fine! Consider instead creating a new topic
      on the [FieldTrip discussion list](http://www.fieldtriptoolbox.org/discussion_list);
      others can then comment on your observation and determine the
      appropriate level of escalation.

## Requesting a new feature

Please search the [GitHub issue list](https://github.com/fieldtrip/fieldtrip/issues)
to see if anybody else has made a comparable request:

   -  If a corresponding issue already exists, please add a comment to that
      issue to escalate the request. Additionally, describe any
      aspect of that feature not yet described in the existing issue.

   -  If no such listing exists, then you are welcome to create a [new
      issue](https://github.com/fieldtrip/fieldtrip/issues/new) outlining the
      request. Be sure to select the "feature request" option to get started
      with writing the issue.

## Asking questions

General questions regarding FieldTrip download, usage, or any other
aspect that is not specific to the FieldTrip *code*, should be directed to
the [discussion list](http://www.fieldtriptoolbox.org/discussion_list). Also check
the [online documentation](https://fieldtriptoolbox.org/), and specifically the [FAQ](https://fieldtriptoolbox.org/faq), in case your issue has already been described there.

## Making direct contributions

Thanks for your interest in making direct contributions to FieldTrip!
We are excited to expand the breadth of researchers involved in improving
and expanding this software, and to ensure that all who make such
contributions receive appropriate acknowledgement through Git.

The instructions below give an short overview of how to go about generating a
proposed change to FieldTrip. A more detailed tutorial on using Git and contributing
to the code (or website) is presented as [online tutorial](http://www.fieldtriptoolbox.org/development/git/)
on the FieldTrip website.

1. You will need to create a *fork* of the [FieldTrip repository](https://github.com/fieldtrip/fieldtrip)
   into your GitHub account, where unlike the main FieldTrip repository,
   you will have full write access to make the requisite changes.

2. Create a Git branch that is named appropriately according to the
   modifications that are being made. The existing code branch on which
   this new derived branch should be based depends on the nature of the
   proposed change (described later below).

3. Generate one or more Git commits that apply your proposed changes to
   the repository:

   -  Individual commits should ideally have a clear singular purpose,
      and not incorporate multiple unrelated changes. If your proposed
      changes involve multiple disparate components, consider breaking
      those changes up into individual commits.

      Conversely, if multiple code changes are logically grouped with /
      linked to one another, these should ideally be integrated into a
      single commit.

   -  Commits should contain an appropriate message that adequately
      describes the change encapsulated within.

      If the change demands a longer description, then the commit message
      should be broken into a synopsis (less than 80 characters) and message
      body, separated by two newline characters (as this enables GitHub to
      parse them appropriately).

      This can be achieved at the command-line as follows:

      `$ git commit -m $'Commit synopsis\n\nHere is a much longer and wordier description of my proposed changes that doesn\'t fit into 80 characters.\nI can even spread the message body across multiple lines.'`

      (Note also the escape character "`\`" necessary for including an
      apostrophe in the message text)

   -  Where relevant, commit messages can also contain references to
      GitHub issues or pull requests (type the "`#`" character followed
      by the issue / PR number), and/or other individual commits (copy
      and paste the first 8-10 characters of the commit hash).

   -  If multiple persons have contributed to the proposed changes, it is
      possible to modify individual Git commits to have [multiple
      authors](https://help.github.com/en/articles/creating-a-commit-with-multiple-authors),
      to ensure that all contributors receive appropriate acknowledgement.

   As a general rule: Git commits and commit messages should be constructed
   in such a way that, at some time in the future, when one is navigating
   through the contribution history, the evolution of the code is as clear
   as possible.

4. Identify the appropriate classification of the change that you propose
   to make, and read the relevant instructions there:

   -  "[**Fix**](#fix)": If the current code behaviour is
      *clearly incorrect*.

   -  "[**Enhancement**](#enhancement)": If the proposed change improves the *performance* or extends the *functionality* of a particular command or process.

   -  "[**Documentation**](#documentation)": If you made some changes to the function description in the help section of the function.

5. Check that your modified code does not prevent FieldTrip from
   passing existing tests, wherever possible (all test files are in the FieldTrip test directory).

6. For code contributions, if possible, a unit test or reproducibility
   test should be added. This not only demonstrates the behaviour of the
   proposed code, but will also preclude future regression of the behaviour
   of that code. For tests that require some custom data, please
   see [here](http://www.fieldtriptoolbox.org/faq/how_should_i_send_example_data_to_the_developers/)
   and/or coordinate with one of the FieldTrip developers.

1. Once completed, a Pull Request should be generated that merges the
   corresponding branch in your forked version of the FieldTrip repository
   into the appropriate target branch of the main FieldTrip repository
   itself:

   -  If your intended changes are complete, and you consider them ready
      to be reviewed by an FieldTrip team member and merged imminently,
      then create a standard Pull Request.

   -  If your changes are ongoing, and you are seeking feedback from the
      FieldTrip developers before completing them, then create a
      [draft Pull Request](https://github.blog/2019-02-14-introducing-draft-pull-requests/).

#### Fix

1. If there does not already exist a [GitHub issue](https://github.com/fieldtrip/fieldtrip/issues)
   describing the bug, consider reporting the bug as a standalone issue
   prior to progressing further; that way developers can confirm the issue,
   and possibly provide guidance if you intend to resolve the issue yourself.
   Later, the Pull Request incorporating the necessary changes should then
   reference the listed issue (simply add somewhere in the description
   section the "`#`" character followed by the issue number).

2. Bug fixes are merged directly to `master`; as such, modifications to the
   code should be made in a branch that is derived from `master`, and the
   corresponding Pull Request should select `master` as the target branch
   for code merging.

3. A unit test or reproducibility test should ideally be added: such a
   test should fail when executed using the current `master` code, but pass
   when executed with the proposed changes.

#### Enhancement

1. New features, as well as any code changes that extend the functionality of
   FieldTrip, are merged to the `master` branch, which contains
   all resolved changes since the most recent tag update. As such, any
   proposed changes that fall under this classification should be made
   in a branch that is based off of the `master` branch, and the corresponding
   Pull Request should select `master` as the target branch for code merging.

#### Documentation

If you want to contribute to the documentation on the Fieldtrip website, please refer to the website's [Github repository](https://github.com/fieldtrip/website).

#### Coding conventions

Please follow the [code guidelines](http://www.fieldtriptoolbox.org/development/guideline/code/) and the [FieldTrip website guidelines pages](http://www.fieldtriptoolbox.org/tag/guidelines/) for guidelines on contributing to the FieldTrip code base.

Non-exhaustive summary of coding guidelines:

* 2 space indents; indent using spaces and not tabs
* No spaces between function name and opening parenthesis of argument list
* One space after the comma that separates function arguments
* Vertically align code with spaces in case it improves readability

#### References

This document is based on the excellent CONTRIBUTING.md document from the [MRTRIX repository](https://github.com/MRtrix3/mrtrix3/), and adjusted accordingly. 
