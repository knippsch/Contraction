~/ExternLibraries/lime/examples/lime_extract_record C2_pi+-_conf0714.dat 2 3 new.dat
od -t fD new.dat > new_ASCII.dat
diff new_ASCII.dat old_ASCII.dat
