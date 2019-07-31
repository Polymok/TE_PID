# Transfer Entropy Partial Information Decomposition

Compute partial information decomposition (PID) on transfer entropy for an input matrix of time-series for a single trial, or an input cell containing multiple matrices corresponding to multiple trials. Transfer entropy is normalized by the entropy of the target time-series. The redundancy partial information term is given by the minimum information function described by Timme et al., 2016.

## Prerequisites

* MATLAB: all functions found here are `.m` files.

## Usage

To calculate transfer entropy PID for all possible neuron triplets, call `TE_PID.m` with two required arguments: a matrix or cell, and a positive integer time-delay. For an input cell containing multiple matrices for multiple trials, the input cell must be 1-dimensional. Each matrix or cell column should contain the entire time-series of a single neuron. A time-delay is necessary to calculate transfer entropy. Output matrices have 7 columns indicating in increasing order: `target_index`, `source1_index`, `source2_index`, `synergy`, `redundancy`, `unique1`, `unique2`. PID output terms are given in units of bits.

### Examples

Call `TE_PID.m` with test case matrices in `test_logicgates.mat` and time-delay `1`.

```MATLAB
>> cd ~/TE_PID  
>> load test_logicgates.mat  
>> TE_PID(identity, 1)  
...  
ans =  
  1 2 3 0 0 1 0  
  ...  
>> TE_PID(xor, 1)  
...  
ans =  
  1 2 3 1 0 0 0  
  ...  
>> TE_PID(and, 1)  
...  
ans =  
  1 2 3 0.5439  0.3113  0.0724  0.0724  
  ...
```

## Functions

### `TE_PID.m`

`TE_PID.m`(*time-series*, *time-delay*)  
Given input matrices and a scalar time-delay, calculate PID terms for transfer entropy.

For *N* neurons, the number of possible neuron triplets is *N\*(N-1)\*(N-2)/2*. Consequently for large *N*, `TE_PID.m` is computationally intensive.

#### Parameters:

*time-series*: MATLAB matrix or 1-dimensional cell. Columns should contain entire time-series of a given neuron.  
*time-delay*: positive integer scalar.

#### Returns:

*PID_matrix*: MATLAB matrix or 1-dimensional cell. Output matrices have 7 columns indicating in increasing order: `target_index`, `source1_index`, `source2_index`, `synergy`, `redundancy`, `unique1`, `unique2`. PID output terms are given in units of bits.

### Call structure

`TE_PID.m`            calls `I_min_TE.m`, `TE.m`, `TE_2dim.m`  
`I_min_TE.m`          calls `I_spec.m`  
`TE_tripletfinder.m`  calls `TE_timelag.m`  
`TE_timelag.m`        calls `TE.m`  
`TE.m`                calls `cond_MI.m`  
`TE_2dim.m`           calls `cond_MI_2dim.m`

### Unfinished and exploratory

* `ndim_cond_MI.m` is not necessary to call `TE_PID.m`. It is, however, an exploratory work-in-progress to extend conditional mutual information to accept *n*-dimensional input time-series. Such an extension from 2-dimensional inputs has proven non-trivial.

* `TE_tripletfinder.m` finds functional neuron triplets `{i,j,k}` that satisfy `TE(j->i)` and `TE(k->i)` greater than some significance value. This function is unfinished and has yet to be incorporated into `TE_PID.m`.

* `MI.m` computes the mutual information between two input time-series. This function is not required to call `TE_PID.m`.

## Bugs

* `I_spec.m`, `cond_MI.m` and `cond_MI_2dim.m` may condition on events with zero probability. When zero probability events are encountered, a notification is printed to the MATLAB console. For now, all cases involving zero probability are discarded. Such a decision may require theoretical justification. Possible alternative solutions include:
  * Time-binning. For Yuqing's model, all AdEx neurons have some time constant `tau`. Might make sense to bin at time resolution equal to `tau`.
  * Choosing an arbitrarily small value in place of probability zero, e.g. `eps` in MATLAB.

* `TE.m` and `TE_2dim.m` calculate the transfer entropy normalized by the entropy of the target. If target entropy is zero, a notification is printed to the MATLAB console. For now, all cases involving zero target entropy use instead the unnormalized transfer entropy.

* Non-zero values (both positive and negative) with small magnitude may obtain when value zero is expected. Possible reasons include:  
  * Rounding error may occur. Information measures in `I_spec.m`, `cond_MI.m`, and `cond_MI_2dim.m` involve both ratios of probabilities as well as the absolute value of probabilities. Ratio calculations use the number of occurrences unnormalized by the total length of the time-series, whereas the absolute value probabilities are divided by the total length.
  * Transfer entropy is normalized by the entropy of the target_future time-series in `TE.m` and `TE_2dim.m` to calculate transfer entropy terms in `TE_PID.m`. This may lead to small discrepancies between the redundancy and the transfer entropy between the individual sources and the target in some cases.

## References

Williams, Paul L. and Randall D. Beer. "Nonnegative Decomposition of Multivariate Information." *CoRR* abs/1004.2515 (2010). url: http://arxiv.org/abs/1004.2515v1.

Williams and Beer. "Generalized Measures of Information Transfer." abs/1102.1507 (2011). url: https://arxiv.org/abs/1102.1507v1.

Timme, Nicholas, Wesley Alford, Benjamin Flecker, and John M. Beggs. "Synergy, Redundancy, and Multivariate Information Measures." *J Comput Neurosci* 36(2014): 119-140. doi: 10.1007/s10827-013-0458-4.

Timme et al. "High-Degree Neurons Feed Cortical Computations." *PLoS Comput Biol* 12(5):e1004858 (2016). doi: 10.1371/journal.pcbi.1004858.

## Authors

* Mofei Wu.  
mofei@uchicago.edu.  
mwumofei@gmail.com.  
