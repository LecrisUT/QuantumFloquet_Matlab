# Matlab library for Quantum Floquet systems

The main goal of this library is to provide a simple implementation of the
various operations frequently used in quantum Floquet calculations. This
library is not well optimized as other back-ends, but it is an easy to
debug and thinker medium.

This library provides calculations of frequently used physical operations,
post-processing of basic input data, and plotting functions of commonly
plotted variables.

Implemented functionalities:
- Eigenstate operations
  - Calculate average energy and quasi-energy eigenstates
  - Adiabatically continue the Floquet eigenstates
  - Calculate the eigenstates variationally
- Other operations
  - Calculate the energy spectra and the overlap measure between energy
    spectra
  - Calculate the Lindblad master equations and their steady states
- Various QOL functions

### (Driven) Physical systems implemented

- Two-level system
- 1D Harmonic oscillator
- 1D Particle in a box

All systems are represented in their discrete static eigenstate
representations. For continuous systems we include projections to
continuous position representations at fixed points.

-------------

#### Contributions

Contributions are always welcome via Github pull requests or discussions in
the issues section.

Some welcome help are in:
1. Documentation: The main documentations should be in HTML format, but
   equivalent short documentations should be included in the code as well.
2. Performance: This must not come at the cost of readablity and ease of
   modification. Instead simply improving the matlab code according to
   matlab best practices, e.g. use `pagemtimes` instead of for loops
3. New physical systems
4. Additional functionalities for common Floquet operations
5. Export to HDF5 and other efficient well structured formats
6. Add convenience plot and post-processing functions
7. Link with other back-ends

#### TODO for v1.0

- [ ] Package as a toolbox
- [ ] Publish HTML documentation

#### How to cite

Cite either:
- PhD thesis: 
- This Github pade and the papers that developed the theory behind these methods:
  - [PRA.105.052213](https://doi.org/10.1103/PhysRevA.105.052213)

If more formal publications are made, this section will be updated to reflect that.