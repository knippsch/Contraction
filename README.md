Markus new contraction code 2014

in LapHs.cpp
- modularized building the propagator and Gamma structure into
  ../modules/Basicoperator.cpp
  new Initialization of new class ReadWrite (see below)
  introduced init_operator() to set up D_u^-1
  introduced get_operator_charged() for u quark and get_operator_g5() to 
  implement gamma_5 trick and gamma structure
- implemented C2 and C4
- transferred read_eigenvectors_from_file() to build_source_matrix()
- implemented get_operator_uncharged for neutral particles
- introduced Corr to hold traces for all momenta, dirac structures, times
  and randomvectors. Build C2, C4_1, C4_2 from Corr to save ca. 25% time
- changed all memory allocation to boost
- reworked output format. Now creates a new folder for every GEVP entry
  (different p^2, dirac, displacement) and writes correlation function
  configurationwise

in Basicoperator.cpp:
- changed layout of perambulator and basicoperator to reflect blockdiagonal
  structure in Dirac space. 
- changed memory layout of basicoperator to array of 4 blocks (rest is 0)
- introduced contraction[rnd_i][blocknr] array holding the 16 blocks of D_u^-1
- changed gamma structure from multiplication to reordering for efficiency
- deleted create_operator(), mul_l_gamma() and mul_r_gamma()
- introduced a struct for lookup tables
- transferred build_source_matrix() to ReadWrite.cpp
- changed create_gamma() into a lookup table
- re-implemented all 16 Dirac structures
- implemented char as switch between interlace and block time dilution
- implemented displacement (not working yet!)
- introduced s matrix and changed from dilution in build_source_matrix to
  "live" dilution in get_operator()
- changed all memory allocation to boost

in ReadWrite.cpp:
- introduced new class ReadWrite to handle the file input
- shifted read_eigenvectors_from_file(), read_perambulators_from_file() and 
  read_rnd_vec_from_file() into ReadWrite.cpp
- changed file layout to Bastians 2014 perambulator code output structure
- implemented momenta in build_source_matrix()
- implemented displacement in build_source_matrix()
- transferred dilution from build_source_matrix() to get_operator() in
  BasicOperator.cpp
- changed read_eigenvectors_from_file() to timeslice-wise reading of 
  vectors
- changed all memory allocation to boost

in config_utils.cpp
- introduced from Christopher for displacement

in GlobalData.cpp
- added dirac, displacement, lattice_name (A40, B55, etc.) and outpath
  to infile

- introduced create_runs.sh and start_runs.sh to create infiles and
  submit jobs to cluster

TODO:
- implement switchable solution for interlace and block time dilution
  in ReadWrite.cpp -> choose dilution scheme in input file?
  Change from TI2 to TB24 layout
- check if create_runs.sh can be writen better in perl (or python)
- check if change in get_operator in BasicOperator.cpp results in 
  speedup
- find out why only the code has only 200% CPU-usage
- write new class Contractions.cpp with code for 2- and 4-pt functions

- implement disconnected diagramms? not necassary?
- IMPLEMENT RHO MESON

known bugs:
- imaginary part always has additional minus sign compared to Bastian
- two point function switch gamma_10 & gamma_11 (solved?)
- C4_2 doesnt cover timevalue t_sink = Lt-1, t_sink_1 = 0 (done for
  consistency with Bastians Code)



