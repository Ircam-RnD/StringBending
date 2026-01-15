# TODO List String Bending project

## 1-Real-time implementation

- [x] Fix initialization of audio processor.
- [ ] Introduce a preset system in the Max patch to remember good sounding configurations. 
- [ ] Find a way to share some info between MAX/MSP objects (Nmodes).
- [ ] Decide wether to keep the architecture as is (multiple objects, each capable of audio operations), or instead use multiple objects for symbolic representation and centralize audio computations in one processor (could be useful for bidirectional interactions: e.g. bowing). Or find a way to keep both.
- [x] Set the nonlinearity type as an object member.

## 2-Numerical method

- [ ] Convergence curves in multiple cases.
- [ ] Perceptive metrics? Maybe at some point organize a simple listening test.
- [ ] A python implementation in probably needed for all that...
- [ ] Compensate for warping on linear part. (sr / pi -> sr / 2 at stability condition).



## 3-Theory
- [ ] Fix modal amplitudes, or provide a control parameter mapping to allow for easy setting of modal amplitudes in dB.
- [ ] Do we want to listen to displacement or velocity or other? Maybe set to displacement, to easily handle nonlinearity levels defined from it.
- [ ] Cleanly write the nonlinearity corresponding to the string (KC and geometric).
- [ ] Re-write the scheme cleanly, with the stability condition.  
- [ ] Design a method for modulation of linear parameters.
