#Parse SMILES according to http://www.opensmiles.org/opensmiles.html
#Copyright (C) 2013 Laurent Fasnacht
#
#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; version 2 of the License.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along
#with this program; if not, write to the Free Software Foundation, Inc.,
#51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

lappend auto_path [file dirname [file dirname [file normalize [info script]]]]
package require smiles_parser

set data [::smiles_parser::smiles_parse "C\[C@@\](C)(O1)C\[C@@H\](O)\[C@@\]1(O2)\[C@@H\](C)\[C@@H\]3CC=C4\[C@\]3(C2)C(=O)C\[C@H\]5\[C@H\]4CC\[C@@H\](C6)\[C@\]5(C)Cc(n7)c6nc(C\[C@@\]89(C))c7C\[C@@H\]8CC\[C@@H\]%10\[C@@H\]9C\[C@@H\](O)\[C@@\]%11(C)C%10=C\[C@H\](O%12)\[C@\]%11(O)\[C@H\](C)\[C@\]%12(O%13)\[C@H\](O)C\[C@@\]%13(C)CO"]
#set data [::smiles_parser::smiles_parse "C1CC1C"]
puts [::smiles_parser::flatten $data]
exit

set fp [open "test_smiles/smiles.txt" r]
set file_data [read $fp]
close $fp

foreach line [split $file_data "\n"] {
    puts $line
    set parsed [::smiles_parser::smiles_parse $line]
    puts $parsed
    set parsed_flat [::smiles_parser::flatten $parsed]
    puts $parsed_flat
}
exit


exit 
puts $data
set atomcount [lindex $data $::smiles_parser::chain_atomcount]
set bondcount [lindex $data $::smiles_parser::chain_bondcount]
set chaincount [lindex $data $::smiles_parser::chain_chaincount]

#Check that we can get all atoms, and that it's correct
for {set i 0} {$i<$atomcount} {incr i} {
    set atom [::smiles_parser::smiles_get_atom $data $i]
    if {[lindex $atom $::smiles_parser::atom_id] != $i} {
        puts "Error with i=$i"
        puts $atom
    }
}

for {set i 0} {$i<$bondcount} {incr i} {
    set bond [::smiles_parser::smiles_get_bond $data $i]
    if {[lindex $bond $::smiles_parser::bond_id] != $i} {
        puts "Error with i=$i"
        puts $bond
    }
}

for {set i 0} {$i<$chaincount} {incr i} {
    puts $i
    set chain [::smiles_parser::smiles_get_chain $data $i]
    if {[lindex $chain $::smiles_parser::chain_id] != $i} {
        puts "Error with i=$i"
        puts $chain
    }
}


puts [::smiles_parser::smiles_get_atom $data 0]
puts [::smiles_parser::smiles_get_bond $data 0]

exit
#puts [::smiles_parser::smiles_parse_chain "CNC"]
puts [::smiles_parser::smiles_parse_atom "C"]
puts [::smiles_parser::smiles_parse_atom "c"]
puts [::smiles_parser::smiles_parse_atom "*"]
puts [::smiles_parser::smiles_parse_atom "\[15N\]"]
puts [::smiles_parser::smiles_parse_atom "\[C\]"]
puts [::smiles_parser::smiles_parse_atom "\[CH0+4\]"]
puts [::smiles_parser::smiles_parse_atom "\[CH+3\]"]
puts [::smiles_parser::smiles_parse_atom "\[CH2+2\]"]
puts [::smiles_parser::smiles_parse_atom "\[CH3+\]"]
puts [::smiles_parser::smiles_parse_atom "\[CH4\]"]
puts [::smiles_parser::smiles_parse_atom "\[OH3+\]"]
puts [::smiles_parser::smiles_parse_atom "\[CH4++\]"]
puts [::smiles_parser::smiles_parse_atom "\[OH-\]"]
puts [::smiles_parser::smiles_parse_atom "\[O--\]"]
puts [::smiles_parser::smiles_parse_atom "\[O-2\]"]
puts [::smiles_parser::smiles_parse_atom "\[C@H\]"]
puts [::smiles_parser::smiles_parse_atom "\[C@H:0\]"]
puts [::smiles_parser::smiles_parse_atom "\[C@H:3\]"]
puts [::smiles_parser::smiles_parse_atom "\[C@H:13\]"]

set fp [open "test_smiles/smiles.txt" r]
set file_data [read $fp]
close $fp

foreach line [split $file_data "\n"] {
    if {[string length $line] > 15} {
        break
    }
    set parsed [::smiles_parser::smiles_parse $line]
    puts "$line $parsed"
    if {0} {
        foreach bracket_atom [regexp -all -inline "\\\[\[^\\\[\]+\\\]" $line] {
            if {[string length $bracket_atom] > 4} {
                set parsed_bracket_atom [::smiles_parser::smiles_parse_atom $bracket_atom]
                puts "$bracket_atom $parsed_bracket_atom"
            }
        }
    }
  # do some line processing here
}
