Markus new contraction code 2014

- modularized building the propagator and Gamma structure into
  ../modules/Basicoperator.cpp
  new Initialization of class new class ReadWrite (see below)
  introduced init_operator() to set up D_u^-1
  introduced get_operator() to implement gamma_5 trick and gamma structure

in Basicoperator.cpp:
- changed layout of perambulator and basicoperator to reflect blockdiagonal
  structure in Dirac space. 
- changed memory layout of basicoperator to array of 4 blocks (rest is 0)
- introduced contraction[rnd_i][blocknr] array holding the 16 blocks of D_u^-1
- changed gamma structure from multiplication to reordering for efficiency
- deleted create_operator(), mul_l_gamma() and mul_r_gamma()
- introduced a struct for lookup tables
- changed create_gamma() into a lookup table
- re-implemented all 16 Dirac structures (not tested yet)
- implemented char as switch between interlace and block time dilution

in ReadWrite.cpp:
- introduced new class ReadWrite to handle the file input
- shifted read_eigenvectors_from_file(), read_perambulators_from_file() and 
  read_rnd_vec_from_file() into ReadWrite.cpp


TODO: implement 4-pt functions in init_operator() and LapHs.cpp
TODO: implement momenta and check dispersion relation

TODO: implement disconnected diagramms
TODO: IMPLEMENT RHO MESON



