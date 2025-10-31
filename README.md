# hmetad_stan

An implementation of the (hierarchical) meta-d' model of confidence
ratings in [Stan](https://mc-stan.org).

This software is adapted from the [HMeta-d
Toolbox](https://github.com/metacoglab/HMeta-d), which was a previous
implementation in JAGS. Beyond the change in estimation software
(which should allow for faster/more reliable model fitting and access
to tools surrounding `cmdstanr`), updates include:

  - Aggregation of confidence data into counts is done automatically
    (previously, this was a mandatory pre-processing step).
  - Automatic generation of estimated pseudo-type 1 ROC and Type 2 ROC
    curves
  - Automatic calculation of type-1 and type-2 response probabilities
  - Provided functions for simulating data from the meta-d' model
  - Minor changes in the structure of the statistical models

Currently, the most significant change is in the specification of the
type-2 confidence thresholds. Previously, these thresholds were given
a normal prior and the ordering constraint was enforced by sorting the
thresholds. In this version, the type-2 confidence thresholds are
ordered using the [ordered
transform](https://mc-stan.org/docs/reference-manual/transforms.html#ordered-vector).
So, instead of enforcing a prior directly on the type-2 confidence
thresholds, this implementation uses a log-normal prior on the
difference between subsequent thresholds.

As of now, there are four implementations, all in the `stan`
directory:
 
 1. `meta_d.stan`: the basic meta-d' model
 2. `meta_d_condition.stan`: the basic meta-d' model, estimated for
    multiple between-participant conditions
 3. `hmeta_d.stan`: the hierarchical meta-d' model for multiple
    participants
 4. `hmeta_d_condition.stan`: the hierarchical meta-d' model for
    multiple participants and multiple within-participant conditions.

Requests for alternative experimental designs can be sent to
<kevin.o'neill@ucl.ac.uk>.
