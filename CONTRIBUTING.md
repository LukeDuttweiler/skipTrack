# Contributing to skipTrack

While the skipTrack package is small and specific, we welcome contributions from the community! In particular we appreciate:

1. Bug reports or feature requests filed as issues.
1. Suggested changes filed as pull requests.

Questions may be addressed to @LukeDuttweiler through GitHub or by email at lduttweiler@hsph.harvard.edu.

## Issues

Issues may be filed on the [skipTrack GitHub Issue page](https://github.com/LukeDuttweiler/skipTrack/issues).

When filing a feature request, please include a short summary of the functionality desired and, if possible, an example.

When filing a bug report, please include a minimal reproducible example. This example should include:

1. **Packages**. These should be loaded at the top of the script.

1. **Data**. The easiest way to get data for the skipTrack package is through the `skipTrack.simulate()` function. If filing a bug report with data from `skipTrack.simulate()`
please ensure that you use `set.seed()` to ensure that the bug is reproducible.

1. **Code**. Please ensure that your code is well-documented and easy for others to read.

Before submitting the bug report, please open a fresh R session and run your example to ensure that it causes the bug as expected.

## Pull requests

Community suggested changes to skipTrack **must** be addressing a reported issue directly. If you are adding a new feature or fixing a bug, please create an issue as described above
prior to submitting a pull request or submit a change related to an already reported issue.

Pull requests may be created as outlined in [this guide](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request).

In order for a change to be considered for skipTrack the pull request must include:

<!--
* [ ] One paragraph motivating the change in NEWS.
      In parentheses, please reference your github user name and the issue you are addressing:
      `(@LukeDuttweiler, #1234)`
* [ ] Only changes relevant to the addressed issue.
* [ ] Updated documentation with roxygen2 re-run.
* [ ] A test for the addressed bug (if addressing a bug report).
* [ ] A minimal example in documentation.
--->
