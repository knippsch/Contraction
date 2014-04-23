Markus new contraction code 2014

in LapHs.cpp
- modularized building the propagator and Gamma structure into
  ../modules/Basicoperator.cpp
  new Initialization of new class ReadWrite (see below)
  introduced init_operator() to set up D_u^-1
  introduced get_operator() to implement gamma_5 trick and gamma structure
- implemented C2 and C4_1

in Basicoperator.cpp:
- changed layout of perambulator and basicoperator to reflect blockdiagonal
  structure in Dirac space. 
- changed memory layout of basicoperator to array of 4 blocks (rest is 0)
- introduced contraction[rnd_i][blocknr] array holding the 16 blocks of D_u^-1
- changed gamma structure from multiplication to reordering for efficiency
- deleted create_operator(), mul_l_gamma() and mul_r_gamma()
- introduced a struct for lookup tables
- changed create_gamma() into a lookup table
- re-implemented all 16 Dirac structures
- implemented char as switch between interlace and block time dilution

in ReadWrite.cpp:
- introduced new class ReadWrite to handle the file input
- shifted read_eigenvectors_from_file(), read_perambulators_from_file() and 
  read_rnd_vec_from_file() into ReadWrite.cpp
- changed file layout to Bastians 2014 perambulator code output structure

TODO:
- implement 4-pt functions in init_operator() and LapHs.cpp
- implement switchable solution for interlace and block time dilution
  in ReadWrite.cpp -> choose dilution scheme in input file?
- implement momenta and check dispersion relation
- check if change in get_operator in BasicOperator.cpp results in 
  speedup
- implement Dirac structure for Four-Point functions (second instance of 
  basic?)

- implement disconnected diagramms
- IMPLEMENT RHO MESON

known bugs:
- imaginary part always has additional minus sign compared to Bastian
- two point function switch gamma_10 & gamma_11



