#this is correct
puts [lindex $::smiles_parser::atom_id $atom]

#This is WRONG
puts [lindex 1 $atom]
