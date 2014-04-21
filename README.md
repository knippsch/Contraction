Markus new contraction code 2014

- modularized building the propagator and Gamma structure into
  ../modules/Basicoperator.cpp
  introduced init_operator() to set up D_u^-1
  introduced get_operator() to implement gamma_5 trick and gamma structure

in Basicoperator.cpp:
- changed layout of perambulator and basicoperator to reflect blockdiagonal
  structure in Dirac space. 
- changed memory layout of basicoperator to array of 4 blocks (rest is 0)
- introduced contraction[rnd_i][col] array holding the columns of D_u^-1
- changed gamma structure from multiplication to reordering for efficiency
- deleted create_operator(), mul_l_gamma() and mul_r_gamma()

TODO: write I/O class to modularize read_eigenvectors_from_file(), 
      read_perambulators_from_file() and read_rnd_vec_from_file() 
      Check.
      make it subclass of basicoperator
TODO: EXPLAIN DEVIATION 
      Check.
TODO: implement better solution for init_operator() regarding 4-pt functions
TODO: implement lookup-table for get_operator
TODO: implement more flexible solution for dilution scheme
TODO: implement momenta

TODO: implement disconnected diagramms
TODO: IMPLEMENT RHO MESON
