# Transfer Entropy Partial Information Decomposition

Compute partial information decomposition (PID) on transfer entropy for an input matrix of time-series for a single trial, or an input cell containing multiple matrices corresponding to multiple trials. Transfer entropy is normalized by the entropy of the target time-series. The redundancy partial information term is given by the minimum information function described by Timme et al., 2016.

## Prerequisites

* MATLAB: all functions found here are `.m` files.

## Usage

We identify a neuron with its time-series. To calculate transfer entropy PID for all possible neuron triplets, call `TE_PID.m` with two required arguments: a matrix or cell, and a positive integer time-delay. Optionally, supply a list of neuron triplet indices for which to calculate PID. For an input cell containing multiple matrices for multiple trials, the input cell must be 1-dimensional. Each matrix or cell column should contain the entire time-series of a single neuron, i.e. columns should represent neurons while rows represent observations at incremental times. A time-delay is necessary to calculate transfer entropy. Output matrices have 7 columns indicating in increasing order: *target_index*, *source1_index*, *source2_index*, *synergy*, *redundancy*, *unique1*, *unique2*. PID output terms are given in units of bits. Each row of the output matrix contains terms for a specific triplet.

### Examples

Call `TE_PID.m` with test case matrix variables `{identity, xor, and}` stored in `test_logicgates.mat` and time-delay `1`.

```
>> cd ~/TE_PID  
>> load test_logicgates.mat  
>> TE_PID(identity, 1);  
Time bin input time-series? y/n: n  
Enter output file name: identity_results.csv  
...  
>> type identity_results.csv  

Target, Source1, Source2, Synergy, Redundancy, Unique1, Unique2  
1, 2, 3, 1, 1, 0, -1  
3, 1, 2, 0.235339, 0.0689627, 0.000788227, 0.000788227  
2, 3, 1, 0.311278, 0, 0, 0.188722  
>> TE_PID(xor, 1);  
Time bin input time-series? y/n: n  
Enter output file name: xor_results.csv  
...  
>> type xor_results.csv  

Target, Source1, Source2, Synergy, Redundancy, Unique1, Unique2  
1, 2, 3, 1, 0, 0, 0  
3, 1, 2, 0.128529, 0.062128, 0.000175672, 0.0314032  
2, 3, 1, 0.106115, 0.0865255, 0.000244657, -0.0242218  
>> TE_PID(and, 1);  
Time bin input time-series? y/n: n  
Enter output file name: and_results.csv  
...  
>> type and_results.csv  

Target, Source1, Source2, Synergy, Redundancy, Unique1, Unique2  
1, 2, 3, 0.543901, 0.311278, 0.0724104, 0.0724104  
3, 1, 2, 0.264169, 0.257856, 0.00294723, -0.23318  
2, 3, 1, 0.0612781, 0.0487949, -2.77556e-17, 0.155639
```

In all three test variables, the first column contains the target time-series of interest. Therefore, we check the first row of the output matrix where *target_index = 1* to verify that our code runs as expected. Note that since transfer entropy splits the target time-series into a future time-series and a past time-series, test cases have their first column shifted by one.

As expected, all transfer entropy information is found in *unique1* for the `identity` test case, *synergy* for the `xor` relation.

In the case of `and`, small values obtain for *unique1* and *unique2*, contrary to the expected zero value, and the *synergy* is slightly larger than the expected *0.5* bits. This is due to a discrepancy between the unnormalized and normalized transfer entropies. The entropy of the target time-series is  

*H(target) = -0.25\*log(0.25)-0.75\*log(0.75) = 0.811278*.

The unnormalized transfer entropies are  

*TE(source1->target) = TE(source2->target) = 0.311278 = redundancy*  
*TE({source1,source2}->target) = 0.811278*.

Therefore, the normalized transfer entropies are

*normed_TE(source1->target) = normed_TE(source2->target) = 0.383688 > redundancy*  
*normed_TE({source1,source2}->target) = 1*.

## Functions

### `TE_PID.m`

`TE_PID.m`(*time-series*, *time-delay*, *triplet_list*)  
Given input matrices and a scalar time-delay, calculate PID terms for transfer entropy.

For *N* neurons, the number of possible neuron triplets is *N\*(N-1)\*(N-2)/2*. Consequently for large *N*, `TE_PID.m` is computationally intensive if *triplet_list* is not given.

#### Parameters:

*time-series*: MATLAB matrix or 1-dimensional cell. Columns should contain entire time-series of a given neuron.  
*time-delay*: positive integer scalar.  
*triplet_list*: optional *nx3* matrix of neuron indices. The first column should represent the target neuron index. If *time-series* is a cell, then *triplet_list* should either be a cell with equal dimensions—each cell element containing a triplet list—or a matrix—in which case PID calculations for all trials are restricted to triplets contained in the single matrix. If not given, PID is calculated for all possible triplets.

#### Returns:

*PID_matrix*: MATLAB matrix or 1-dimensional cell. Output matrices have 7 columns indicating in increasing order: `target_index`, `source1_index`, `source2_index`, `synergy`, `redundancy`, `unique1`, `unique2`. Each row of the output matrix corresponds to a specific triplet. PID terms are given in units of bits.

### Call structure

| Parent function      | Child function [optional]         |
|----------------------|-----------------------------------|
| `TE_PID.m`           | `I_min_TE.m` `TE.m` [`timebin.m`] |
| `I_min_TE.m`         | `I_spec.m`                        |
| `TE_tripletfinder.m` | `TE_timelag.m`                    |
| `TE_timelag.m`       | `TE.m`                            |
| `TE.m`               | `cond_MI.m`                       |

## Bugs

* `I_spec.m` and `cond_MI.m` may condition on events with zero probability. When zero probability events are encountered, a notification is printed to the MATLAB console. For now, all cases involving zero probability are discarded. Such a decision may require theoretical justification. Possible alternative solutions include:
  * Time-binning. For Yuqing's model, all AdEx neurons have some time constant *tau*. Might make sense to bin at time resolution equal to *tau*.
  * Choosing an arbitrarily small value in place of probability zero, e.g. `eps` in MATLAB.

* `TE.m` calculates the transfer entropy normalized by the entropy of the target. If target entropy is zero, a notification is printed to the MATLAB console. For now, all cases involving zero target entropy use instead the unnormalized transfer entropy.

* Non-zero values (both positive and negative) with small magnitude may obtain when value zero is expected. Possible reasons include:  
  * Rounding error may occur. Information measures in `I_spec.m` and `cond_MI.m` require both ratios of probabilities as well as the absolute value of probabilities. Ratio calculations use the number of occurrences unnormalized by the total length of the time-series, whereas the absolute value probabilities are divided by the total length.
  * Transfer entropy is normalized by the entropy of the target_future time-series in `TE.m` to calculate transfer entropy terms in `TE_PID.m`. This may lead to small discrepancies between the redundancy and the transfer entropy between the individual sources and the target in some cases. We can consider using the unnormalized transfer entropy for all calculations.

## References

Williams, Paul L. and Randall D. Beer. "Nonnegative Decomposition of Multivariate Information." *CoRR* abs/1004.2515 (2010). url: http://arxiv.org/abs/1004.2515v1.

Williams and Beer. "Generalized Measures of Information Transfer." abs/1102.1507 (2011). url: https://arxiv.org/abs/1102.1507v1.

Timme, Nicholas, Wesley Alford, Benjamin Flecker, and John M. Beggs. "Synergy, Redundancy, and Multivariate Information Measures." *J Comput Neurosci* 36(2014): 119-140. doi: 10.1007/s10827-013-0458-4.

Timme et al. "High-Degree Neurons Feed Cortical Computations." *PLoS Comput Biol* 12(5):e1004858 (2016). doi: 10.1371/journal.pcbi.1004858.

## Authors

* Mofei Wu.  
mofei@uchicago.edu.  
mwumofei@gmail.com.  
