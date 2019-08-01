# Transfer Entropy Partial Information Decomposition

Compute partial information decomposition (PID) on transfer entropy for an input matrix of time-series for a single trial, or an input cell containing multiple matrices corresponding to multiple trials. Transfer entropy is normalized by the entropy of the target time-series. The redundancy partial information term is given by the minimum information function described by Timme et al., 2016.

## Prerequisites

* MATLAB: all functions found here are `.m` files.

## Usage

We identify a neuron with its time-series. To calculate transfer entropy PID for all possible neuron triplets, call `TE_PID.m` with two required arguments: a matrix or cell, and a positive integer time-delay. For an input cell containing multiple matrices for multiple trials, the input cell must be 1-dimensional. Each matrix or cell column should contain the entire time-series of a single neuron. A time-delay is necessary to calculate transfer entropy. Output matrices have 7 columns indicating in increasing order: *target_index*, *source1_index*, *source2_index*, *synergy*, *redundancy*, *unique1*, *unique2*. PID output terms are given in units of bits. Each row of the output matrix contains terms for a specific triplet.

### Examples

Call `TE_PID.m` with test case matrix variables `{identity, xor, and}` stored in `test_logicgates.mat` and time-delay `1`.

```
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

In all three test variables, the first column contains the target time-series of interest. Therefore, we check the first row of the output matrix where *target_index = 1* to verify that our code runs as expected. Note that since transfer entropy splits the target time-series into a future time-series and a past time-series, test cases have their first column shifted by one.

As expected, all transfer entropy information is found in *unique1* for the `identity` test case, *synergy* for the `xor` relation.

In the case of `and`, small values obtain for *unique1* and *unique2*, contrary to the expected zero value, and the *synergy* is slightly larger than the expected *0.5* bits. This is due to a discrepancy between the unnormalized and normalized transfer entropies. The entropy of the target time-series is *H(target) = -0.25\*log(0.25)-0.75\*log(0.75) = 0.8113*. The unnormalized transfer entropies are *TE(source1->target) = TE(source2->target) = 0.3113 = redundancy*, and *TE({source1,source2}->target) = 0.8113*. Therefore, the normalized transfer entropies are *normed_TE(source1->target) = normed_TE(source2->target) = 0.3837 > redundancy* and *normed_TE({source1,source2}->target) = 1*.

## Functions

### `TE_PID.m`

`TE_PID.m`(*time-series*, *time-delay*)  
Given input matrices and a scalar time-delay, calculate PID terms for transfer entropy.

For *N* neurons, the number of possible neuron triplets is *N\*(N-1)\*(N-2)/2*. Consequently for large *N*, `TE_PID.m` is computationally intensive.

#### Parameters:

*time-series*: MATLAB matrix or 1-dimensional cell. Columns should contain entire time-series of a given neuron.  
*time-delay*: positive integer scalar.

#### Returns:

*PID_matrix*: MATLAB matrix or 1-dimensional cell. Output matrices have 7 columns indicating in increasing order: `target_index`, `source1_index`, `source2_index`, `synergy`, `redundancy`, `unique1`, `unique2`. Each row of the output matrix corresponds to a specific triplet. PID terms are given in units of bits.

### Call structure

`TE_PID.m`            calls `I_min_TE.m`, `TE.m`, `TE_2dim.m`  
`I_min_TE.m`          calls `I_spec.m`  
`TE_tripletfinder.m`  calls `TE_timelag.m`  
`TE_timelag.m`        calls `TE.m`  
`TE.m`                calls `cond_MI.m`  
`TE_2dim.m`           calls `cond_MI_2dim.m`

### Unfinished

* `TE_tripletfinder.m` finds functional neuron triplets *{i,j,k}* that satisfy *TE(j->i)* and *TE(k->i)* greater than some significance value. This function is unfinished and has yet to be incorporated into `TE_PID.m`.

### Exploratory

* `MI.m` computes the mutual information between two scalar-valued, discrete time-series. This function is not required to call `TE_PID.m`.

* `ndim_cond_MI.m` computes the conditional mutual information for vector-valued, discrete time-series. This function is not required to call `TE_PID.m`. It is possible in the future to replace both `cond_MI.m` and `cond_MI_2dim.m` by a function that allows *n*-dimensional inputs, namely `ndim_cond_MI.m`. Similarly, `TE.m` and `TE_2dim.m` may be subsumed under a single function.

## Bugs

* `I_spec.m`, `cond_MI.m` and `cond_MI_2dim.m` may condition on events with zero probability. When zero probability events are encountered, a notification is printed to the MATLAB console. For now, all cases involving zero probability are discarded. Such a decision may require theoretical justification. Possible alternative solutions include:
  * Time-binning. For Yuqing's model, all AdEx neurons have some time constant `tau`. Might make sense to bin at time resolution equal to `tau`.
  * Choosing an arbitrarily small value in place of probability zero, e.g. `eps` in MATLAB.

* `TE.m` and `TE_2dim.m` calculate the transfer entropy normalized by the entropy of the target. If target entropy is zero, a notification is printed to the MATLAB console. For now, all cases involving zero target entropy use instead the unnormalized transfer entropy.

* Non-zero values (both positive and negative) with small magnitude may obtain when value zero is expected. Possible reasons include:  
  * Rounding error may occur. Information measures in `I_spec.m`, `cond_MI.m`, and `cond_MI_2dim.m` involve both ratios of probabilities as well as the absolute value of probabilities. Ratio calculations use the number of occurrences unnormalized by the total length of the time-series, whereas the absolute value probabilities are divided by the total length.
  * Transfer entropy is normalized by the entropy of the target_future time-series in `TE.m` and `TE_2dim.m` to calculate transfer entropy terms in `TE_PID.m`. This may lead to small discrepancies between the redundancy and the transfer entropy between the individual sources and the target in some cases. We can consider using the unnormalized transfer entropy for all calculations.

## References

Williams, Paul L. and Randall D. Beer. "Nonnegative Decomposition of Multivariate Information." *CoRR* abs/1004.2515 (2010). url: http://arxiv.org/abs/1004.2515v1.

Williams and Beer. "Generalized Measures of Information Transfer." abs/1102.1507 (2011). url: https://arxiv.org/abs/1102.1507v1.

Timme, Nicholas, Wesley Alford, Benjamin Flecker, and John M. Beggs. "Synergy, Redundancy, and Multivariate Information Measures." *J Comput Neurosci* 36(2014): 119-140. doi: 10.1007/s10827-013-0458-4.

Timme et al. "High-Degree Neurons Feed Cortical Computations." *PLoS Comput Biol* 12(5):e1004858 (2016). doi: 10.1371/journal.pcbi.1004858.

## Authors

* Mofei Wu.  
mofei@uchicago.edu.  
mwumofei@gmail.com.  
