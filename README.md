# SeedbankTree

This repository contains the SeedbankTree BEAST 2 package, as well as additional testing and data files pertaining to "Bayesian phylodynamic inference of population dynamics with dormancy".

### Usage

As of now, SeedbankTree isn't officially included in the BEAST 2 package manager. To use SeedbankTree, clone this repository and then run `ant build` from inside the directory to build a zip file containing the entire package located in the `dist` folder. Then, create a new "SeedbankTree" subdirectory in a BEAST 2 package directory, the location of which depends on your OS. If you are on Linux and using version 2.7 of BEAST 2, the default location to put the new "SeedbankTree" subdirectory is `~/.beast/2.7/SeedbankTree`. Copy the zip file into this directory, unzip it, and you should be good to go!

Check out the "Install by hand" section under [this](https://www.beast2.org/managing-packages/) BEAST 2 article for more detailed information.

### Files

#### src

`src` contains all the java code that makes up the SeedbankTree package.

- `SeedbankTreeDensity.java` \
  Probability distribution of the seedbank coalescent tree.

- `SeedbankClockModel.java` \
  Branch rate model for computing branch-specific rate multipliers according to the percentage of dormancy on a branch in between two sample/coalescent events.

- `SeedbankTree.java`, `SeedbankNode.java` \
  The seedbank tree data structure.

- `SeedbankTreeFromNewick`, `SeedbankTreeInitializer` \
  Initializers for the first seedbank tree in a MCMC run. `SeedbankTreeFromNewick` builds a `SeedbankTree` corresponding to an annotated newick string, whereas `SeedbankTreeInitializer` simulates a random `SeedbankTree` from the number of samples.

- `TransitionModel.java` \
  Model for transitions between the active and dormant states in the seedbank coalescent.

- `operators\` \
  Folder containing all the MCMC operators. `SeedbankTreeOperator` and `UniformizationRetypeOperator` are parent classes that contain general topology helper functions and retype helper functions respectively.

- `util\` \
  Folder containing some utility loggers.

#### validation_data

- `analytic_validation` \
  This folder contains the log files and example XML files for the analytic validation conducted in the paper. The four subdirectories correspond to the four sets of analytic validation: 2-sample trees using only `NodeShiftRetype`, 2-sample trees using `SeedbankTreeScale` and `RecolorBranch`, and 10-sample and 100-sample trees using all operators. `TABLE_VALUES.txt` compiles all the analytic and MCMC values used in the tables for all the statistics apart from active and dormant length percentage. `length_percentages.ipynb` calculates the active and dormant length percentages.

- `synthetic_validation` \
  This folder contains the log files and utility scripts for the synthetic validation conducted in the paper. The run log data compiled into the 95% HPD statistics are contained in the four `model_X_logs` directories. Model 1 is c=K=1, alpha=0.1; model 2 is c=K=1, alpha=0.5; model 3 is c=1, K=10, alpha=0.1; model 4 is c=1, K=10, alpha=0.5. Note that more log files are included than are actually used in the statistics—we generate batches of XML files and run them in parallel, but look at runs individually and sequentially. \
  Although the XML files are not included, the `make_batch_runs` subdirectory contains the utility scripts used to automate the pipeline of making synthetic data and creating XML files. The modified `treetime` package used to mutate nucleotide sequences is included in `make_batch_runs`. Feel free to use the functions in `simulate_tree.py` and `simulate_mutation.py` files directly to generate synthetic data of your own.
