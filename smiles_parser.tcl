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


#Structures:
#['atom'...]
# symbol: one of element_symbols
# aromatic: 1/0
# isotope: number, -1 if not specified
# chiral: '', or one of chiral
# hcount: number, -1 if not specified
# charge: number (0 by default)
# class: number (0 by default)

#['bond'...]
# count: number (1-4)
# aromatic: 1/0/-1 (-1 = unknown)
# direction: '/' '\' '-'

#['chain'...]

#['ringbond'...]
#ringbond: "number|ref" (we can use multiple time the same number, buf ref will increase)

#return of a parser func: [length data]

package provide smiles_parser 1.0



namespace eval ::smiles_parser:: {
    #Object definition for atom
    set atom_default [list "atom" -1 -1 -1 -1 "*" 0 -1 "" 0 0 0 [list] [list] [list]]
    
    set atom_id 1
    set atom_atomcount 2
    set atom_bondcount 3
    set atom_chaincount 4
    set atom_symbol 5
    set atom_aromatic 6
    set atom_isotope 7
    set atom_chiral 8
    set atom_hcount 9
    set atom_charge 10
    set atom_class 11
    set atom_ringbonds 12
    set atom_branches 13
    set atom_linkedby 14
    
    #Object definition for bond
    set bond_default [list "bond" -1 0 -1 "-" [list]]
    
    set bond_id 1
    set bond_count 2
    set bond_aromatic 3
    set bond_direction 4
    set bond_linking 5
    
    #Object definition for chain
    set chain_default [list "chain" -1 -1 -1 -1 [list]]
    
    set chain_id 1
    set chain_atomcount 2
    set chain_bondcount 3
    set chain_chaincount 4
    set chain_list 5
    
    #Object definition for ringbond
    set ringbond_default [list "ringbond" -1 -1 -1 $bond_default]
    
    set ringbond_to_atom 1
    set ringbond_number 2
    set ringbond_ref 3
    set ringbond_bond 4
    
    
    #Temporary variables
    set chain_refs [list]
    set cur_id_atom 0
    set cur_id_bond 0
    set cur_id_chain 0
    
    #for numbering
    set cur_ringbond_links [dict create]
    
    set cur_flatten_atom [list]
    set cur_flatten_bond [list]
    set cur_flatten_chain [list]
    
    
    #Regular expression parts (from grammar)
    set re_aliphatic_organic "B|C|N|O|S|P|F|Cl|Br|I|\\*"
    set re_aromatic_organic "b|c|n|o|s|p"
    set re_isotope "\[0-9\]*"
    set re_element_symbols "H|He|Li|Be|B|C|N|O|F|Ne|Na|Mg|Al|Si|P|S|Cl|Ar|K|Ca|Sc|Ti|V|Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Ge|As|Se|Br|Kr|Rb|Sr|Y|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Sb|Te|I|Xe|Cs|Ba|Hf|Ta|W|Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|At|Rn|Fr|Ra|Rf|Db|Sg|Bh|Hs|Mt|Ds|Rg|Cn|Fl|Lv|La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Ac|Th|Pa|U|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr"
    set re_aromatic_symbols "b|c|n|o|p|s|se|as"
    set re_chiral "|@|@@|@TH1|@TH2|@AL1|@AL2|@SP1|@SP2|@SP3|@TB1|@TB2|@TB3|@TB4|@TB5|@TB6|@TB7|@TB8|@TB9|@TB10|@TB11|@TB12|@TB13|@TB14|@TB15|@TB16|@TB17|@TB18|@TB19|@TB20|@OH1|@OH2|@OH3|@OH4|@OH5|@OH6|@OH7|@OH8|@OH9|@OH10|@OH11|@OH12|@OH13|@OH14|@OH15|@OH16|@OH17|@OH18|@OH19|@OH20|@OH21|@OH22|@OH23|@OH24|@OH25|@OH26|@OH27|@OH28|@OH29|@OH30"
    set re_hcount "H?|H\[0-9\]+"
    set re_charge "-?|\\+?|-\[0-9\]|\\+\[0-9\]+|--|\\+\\+"
    set re_class "|:\[0-9\]+"
    set re_bond "|-|=|#|\$|:|/|\\\\"
    set re_ringbond_num "\[0-9\]|%\[0-9\]\[0-9\]"
}

#proc ::smiles_parser::smiles_parse_thing {SMILES} {
#    set l 0
#    set ret [list "thing"]
#    
#    return [list $l $ret]
#}

proc ::smiles_parser::smiles_parse_bond {SMILES} {
    #length parsed
    set length 1
    #['bond' count aromatic direction]
    set ret $::smiles_parser::bond_default
    
    set bond_desc [string range $SMILES 0 0]
    
    if {$bond_desc == "-"} {
        lset ret $::smiles_parser::bond_count 1
        lset ret $::smiles_parser::bond_aromatic 0
    } elseif {$bond_desc == "/"} {
        lset ret $::smiles_parser::bond_count 1
        lset ret $::smiles_parser::bond_aromatic 0
        lset ret $::smiles_parser::bond_direction "/"
    } elseif {$bond_desc == "\\"} {
        lset ret $::smiles_parser::bond_count 1
        lset ret $::smiles_parser::bond_aromatic 0
        lset ret $::smiles_parser::bond_direction "\\"
    } elseif {$bond_desc == "="} {
        lset ret $::smiles_parser::bond_count 2
        lset ret $::smiles_parser::bond_aromatic 0
    } elseif {$bond_desc == "#"} {
        lset ret $::smiles_parser::bond_count 3
        lset ret $::smiles_parser::bond_aromatic 0
    } elseif {$bond_desc == "\$"} {
        lset ret $::smiles_parser::bond_count 4
        lset ret $::smiles_parser::bond_aromatic 0
    } elseif {$bond_desc == ":"} {
        lset ret $::smiles_parser::bond_count 1
        lset ret $::smiles_parser::bond_aromatic 1
    } else {
        lset ret $::smiles_parser::bond_count 1
        set length 0
    }
    return [list $length $ret]
}

proc ::smiles_parser::smiles_parse_atom {SMILES} {
    #length parsed
    set length 0
    #['atom' symbol aromatic isotope chiral hcount charge class [list of ringbonds] [list of branches (chains)]]
    set ret $::smiles_parser::atom_default
    #match aliphatic_organic
    if { [regexp "^(${::smiles_parser::re_aliphatic_organic})" $SMILES match] } {
        lset ret $::smiles_parser::atom_symbol $match
        set length [string length $match]
        set SMILES [string range $SMILES $length end]
    #match aromatic_organic
    } elseif { [regexp "^(${::smiles_parser::re_aromatic_organic})" $SMILES match] } {
        lset ret $::smiles_parser::atom_symbol [string toupper $match]
        lset ret $::smiles_parser::atom_aromatic 1
        
        set length [string length $match]
        set SMILES [string range $SMILES $length end]
    #match bracket_atom
    } elseif { [regexp "^\\\[(${::smiles_parser::re_isotope})(${::smiles_parser::re_element_symbols}|${::smiles_parser::re_aromatic_symbols}|\\*)(${::smiles_parser::re_chiral})(${::smiles_parser::re_hcount})(${::smiles_parser::re_charge})(${::smiles_parser::re_class})\\\]" $SMILES match m_isotope m_element m_chiral m_hcount m_charge m_class] } {
        if { $m_isotope != "" } {
            lset ret $::smiles_parser::atom_isotope $m_isotope
        }
        lset ret $::smiles_parser::atom_symbol $m_element
        lset ret $::smiles_parser::atom_chiral $m_chiral
        
        #aromatic_symbols?
        if { [regexp "^(b|c|n|o|p|s|se|as)" $m_element]} {
            if { $match == "se"} {
                lset ret $::smiles_parser::atom_symbol "Se"
            } elseif { $match == "as" } {
                lset ret $::smiles_parser::atom_symbol "As"
            } else {
                lset ret $::smiles_parser::atom_symbol [string toupper $match]
            }
            lset ret $::smiles_parser::atom_aromatic 1
        }
        
        if { $m_hcount != "" } {
            if {$m_hcount == "H"} {
                lset ret $::smiles_parser::atom_hcount 1
            } else {
                lset ret $::smiles_parser::atom_hcount [string range $m_hcount 1 end]
            }
        }
        if { $m_charge != "" } {
            if {$m_charge == "-"} {
                lset ret $::smiles_parser::atom_charge -1
            } elseif {$m_charge == "--"} {
                lset ret $::smiles_parser::atom_charge -2
            } elseif {$m_charge == "+"} {
                lset ret $::smiles_parser::atom_charge 1
            } elseif {$m_charge == "++"} {
                lset ret $::smiles_parser::atom_charge 2
            } else {
                lset ret $::smiles_parser::atom_charge [expr $m_charge]
            }
        }
        
        if { $m_class != "" } {
            lset ret $::smiles_parser::atom_class [string range $m_class 1 end]
        }
        
        set length [string length $match]
        set SMILES [string range $SMILES $length end]
    } else {
        #Didn't match anything!
        return [list 0 [list]]
    }
    
    set ringbonds [list]
    
    while {[regexp "^(${::smiles_parser::re_bond})(${::smiles_parser::re_ringbond_num})" $SMILES match m_bond m_number] } {
        #Remove % if needed
        if { [string range $m_number 0 0] == "%" } {
            set m_number [string range $m_number 1 end]
        }
        #get reference counter
        set match_index_rc [lindex $::smiles_parser::chain_refs $m_number]
        #increase it
        lset ::smiles_parser::chain_refs $m_number [expr $match_index_rc + 1]
        
        
        set ringbond $::smiles_parser::ringbond_default
        lset ringbond $::smiles_parser::ringbond_number $m_number
        lset ringbond $::smiles_parser::ringbond_ref [expr $match_index_rc / 2]
        #Bond is the second term of the return value
        lset ringbond $::smiles_parser::ringbond_bond [lindex [::smiles_parser::smiles_parse_bond $m_bond] 1]
        
        lappend ringbonds $ringbond
        
        set length [expr $length + [string length $match]]
        set SMILES [string range $SMILES [string length $match] end]
    }
    
    lset ret $::smiles_parser::atom_ringbonds $ringbonds
    
    set branches [list]
    #Create the subchain for all branches
    while {[string range $SMILES 0 0] == "(" } {
        set parenthese_count 0
        for {set x 1} {$x<[string length $SMILES]} {incr x} {
            if {[string range $SMILES $x $x] == "("} {
                incr parenthese_count
            } elseif {[string range $SMILES $x $x] == ")"} {
                if {$parenthese_count == 0} {
                    break
                } else {
                    incr parenthese_count -1
                }
            }
        }
        set subchain_end [expr $x - 1]
        
        set subchain [string range $SMILES 1 $subchain_end]
        
        set parse_chain_result [::smiles_parser::smiles_parse_chain $subchain 1]
        set chain_length [lindex $parse_chain_result 0]
        set chain_result [lindex $parse_chain_result 1]
        
        if {[string length $subchain] != $chain_length} {
            puts "Could not fully parse subchain $subchain, abort!"
            break
        }
        
        lappend branches $chain_result
        
        set length [expr $length + $subchain_end + 2]
        set SMILES [string range $SMILES [expr $subchain_end + 2] end]
    }
    
    lset ret $::smiles_parser::atom_branches $branches
    
    return [list $length $ret]
}

proc ::smiles_parser::smiles_parse_chain {SMILES start_by_bond} {
    set length 0
    set chain [list]
    set ret $::smiles_parser::chain_default
    
    while {[string length $SMILES] > 0} {
        if {$start_by_bond} {
            set start_by_bond 0
        } else {
            set parse_atom_result [::smiles_parser::smiles_parse_atom $SMILES]
            set atom_length [lindex $parse_atom_result 0]
            set atom_result [lindex $parse_atom_result 1]
            
            #if we get don't get an atom it's a failure!
            if {$atom_length == 0} {
                puts "Parse atom failed (parsed 0 bytes). Remaining SMILES: $SMILES"
                break
            }
            if {[lindex $atom_result 0] != "atom"} {
                puts "Parse atom failed (didn't return an atom). Remaining SMILES: $SMILES"
                break
            }
            lappend chain $atom_result
            
            #Remove atom from unparsed SMILES
            set length [expr $length + $atom_length]
            set SMILES [string range $SMILES $atom_length end]
        }
        
        #If something is left in the SMILES, we have a bond (which may be implicit)
        if {[string length $SMILES] > 0} {
            #We skip the dots (dots are NOT bonds)
            if {[string range $SMILES 0 0] == "."} {
                set length [expr $length + 1]
                set SMILES [string range $SMILES 1 end]
                continue
            }
            set parse_bond_result [::smiles_parser::smiles_parse_bond $SMILES]
            set bond_length [lindex $parse_bond_result 0]
            set bond_result [lindex $parse_bond_result 1]
            
            if {[lindex $bond_result 0] != "bond"} {
                puts "Parse bond failed (didn't return an bond). Remaining SMILES: $SMILES"
                break
            }
            
            lappend chain $bond_result
            #Remove atom from unparsed SMILES
            set length [expr $length + $bond_length]
            set SMILES [string range $SMILES $bond_length end]
        }
    }
    
    lset ret $::smiles_parser::chain_list $chain
    
    return [list $length $ret]
}

proc ::smiles_parser::_chain_renumber {obj} {
    if {[lindex $obj 0] == "chain" } {
        set ret $::smiles_parser::chain_default
        
        #set chain id
        lset ret $::smiles_parser::chain_id $::smiles_parser::cur_id_chain
        incr ::smiles_parser::cur_id_chain
        
        set startatom $::smiles_parser::cur_id_atom
        set startbond $::smiles_parser::cur_id_bond
        set startchain $::smiles_parser::cur_id_chain
        
        set elts [list]
        foreach subelt [lindex $obj $::smiles_parser::chain_list] {
            lappend elts [::smiles_parser::_chain_renumber $subelt]
        }
        
        lset ret $::smiles_parser::chain_list $elts
        
        lset ret $::smiles_parser::chain_atomcount [expr $::smiles_parser::cur_id_atom - $startatom]
        lset ret $::smiles_parser::chain_bondcount [expr $::smiles_parser::cur_id_bond - $startbond]
        lset ret $::smiles_parser::chain_chaincount [expr $::smiles_parser::cur_id_chain - $startchain]
        return $ret
    } elseif {[lindex $obj 0] == "bond" } {
        set ret $obj
        lset ret $::smiles_parser::bond_id $::smiles_parser::cur_id_bond
        incr ::smiles_parser::cur_id_bond
        return $ret
    } elseif {[lindex $obj 0] == "atom" } {
        set ret $obj
        lset ret $::smiles_parser::atom_id $::smiles_parser::cur_id_atom
        incr ::smiles_parser::cur_id_atom
        
        set startatom $::smiles_parser::cur_id_atom
        set startbond $::smiles_parser::cur_id_bond
        set startchain $::smiles_parser::cur_id_chain
        
        set elts [list]
        foreach subelt [lindex $obj $::smiles_parser::atom_branches] {
            lappend elts [::smiles_parser::_chain_renumber $subelt]
        }
        
        lset ret $::smiles_parser::atom_branches $elts
        
        set rbs [list]
        #For ringbonds, we only append RESOLVED ringbonds
        foreach rb [lindex $obj $::smiles_parser::atom_ringbonds] {
            set rb_n [lindex $rb $::smiles_parser::ringbond_number]
            set rb_i [lindex $rb $::smiles_parser::ringbond_ref]
            set rb_id "${rb_n}_${rb_i}"
            
            if {[dict exists $::smiles_parser::cur_ringbond_links $rb_id]} {
                lset rb $::smiles_parser::ringbond_to_atom [dict get $::smiles_parser::cur_ringbond_links $rb_id]
                lappend rbs [::smiles_parser::_chain_renumber $rb]
            } else {
                dict set ::smiles_parser::cur_ringbond_links $rb_id [lindex $ret $::smiles_parser::atom_id]
            }
        }
        
        lset ret $::smiles_parser::atom_ringbonds $rbs
        
        lset ret $::smiles_parser::atom_atomcount [expr $::smiles_parser::cur_id_atom - $startatom]
        lset ret $::smiles_parser::atom_bondcount [expr $::smiles_parser::cur_id_bond - $startbond]
        lset ret $::smiles_parser::atom_chaincount [expr $::smiles_parser::cur_id_chain - $startchain]
        
        return $ret
    } elseif {[lindex $obj 0] == "ringbond" } {
        set ret $obj
        lset ret $::smiles_parser::ringbond_bond [::smiles_parser::_chain_renumber [lindex $obj $::smiles_parser::ringbond_bond]]
        return $ret
    } else {
        set datatype [lindex $obj 0] 
        error "Unsupported data type $datatype";
    }
}

proc ::smiles_parser::_chain_find_neighbors_add_link {atom1 bond atom2} {
    set _bond [lindex $::smiles_parser::cur_flatten_bond $bond]
    set _atom1 [lindex $::smiles_parser::cur_flatten_atom $atom1]
    set _atom2 [lindex $::smiles_parser::cur_flatten_atom $atom2]
    
    #If bond is unspecified, then set aromatic if both ends are aromatic atoms
    if {[lindex $_bond $::smiles_parser::bond_aromatic] == -1} {
        if {[lindex $_atom1 $::smiles_parser::atom_aromatic] == 1 && [lindex $_atom2 $::smiles_parser::atom_aromatic] == 1} {
            lset _bond $::smiles_parser::bond_aromatic 1
        } else {
            lset _bond $::smiles_parser::bond_aromatic 0
        }
    }
    
    #Sanity check: verify that an aromatic bond is between two aromatic atoms
    if {[lindex $_bond $::smiles_parser::bond_aromatic] == 1} {
        if {[lindex $_atom1 $::smiles_parser::atom_aromatic] == 0 || [lindex $_atom2 $::smiles_parser::atom_aromatic] == 0} {
            error "Bond $bond should be between two aromatic atoms ($atom1 $atom2)"
        }
    }
    
    lset _bond $::smiles_parser::bond_linking [list $atom1 $atom2]
    
    set atom1_linkedby [lindex $_atom1 $::smiles_parser::atom_linkedby]
    lappend atom1_linkedby $bond
    lset _atom1 $::smiles_parser::atom_linkedby $atom1_linkedby
    
    set atom2_linkedby [lindex $_atom2 $::smiles_parser::atom_linkedby]
    lappend atom2_linkedby $bond
    lset _atom2 $::smiles_parser::atom_linkedby $atom2_linkedby
    
    lset ::smiles_parser::cur_flatten_bond $bond $_bond
    lset ::smiles_parser::cur_flatten_atom $atom1 $_atom1
    lset ::smiles_parser::cur_flatten_atom $atom2 $_atom2
}

proc ::smiles_parser::_chain_find_neighbors {chain {previous_atom -1}} {
    if {[lindex $chain 0] != "chain" } { error "_chain_find_neighbors needs a chain" }
    
    set previous_bond -1
    
    set new_list [list]
    
    foreach e [lindex $chain $::smiles_parser::chain_list] {
        if {[lindex $e 0] == "bond" } {
            if {$previous_bond != -1} {error "BUG: two bonds following each other?!"}
            if {$previous_atom == -1} {error "BUG: no atom before this bond?!"}
            
            set previous_bond [lindex $e $::smiles_parser::bond_id]
        } elseif {[lindex $e 0] == "atom" } {
            set current_atom [lindex $e $::smiles_parser::atom_id]
            if {$previous_bond != -1} {
                #We link something
                ::smiles_parser::_chain_find_neighbors_add_link $previous_atom $previous_bond $current_atom
                set previous_bond -1
            }
            set previous_atom $current_atom
            
            foreach branch [lindex $e $::smiles_parser::atom_branches] {
                ::smiles_parser::_chain_find_neighbors $branch $current_atom
            }
            foreach ringbond [lindex $e $::smiles_parser::atom_ringbonds] {
                set rb_bond [lindex $ringbond $::smiles_parser::ringbond_bond]
                ::smiles_parser::_chain_find_neighbors_add_link $current_atom [lindex $rb_bond $::smiles_parser::bond_id] [lindex $ringbond $::smiles_parser::ringbond_to_atom]
            }
        }
    }
}

proc ::smiles_parser::flatten {chain} {
    set ::smiles_parser::cur_flatten_atom [list]
    set ::smiles_parser::cur_flatten_bond [list]
    set ::smiles_parser::cur_flatten_chain [list]
    if {[lindex $chain $::smiles_parser::chain_atomcount] > 0} {
        set ::smiles_parser::cur_flatten_atom [lrepeat [lindex $chain $::smiles_parser::chain_atomcount] 0]
    }
    if {[lindex $chain $::smiles_parser::chain_bondcount] > 0} {
        set ::smiles_parser::cur_flatten_bond [lrepeat [lindex $chain $::smiles_parser::chain_bondcount] 0]
    }
    set ::smiles_parser::cur_flatten_chain [lrepeat [expr 1+[lindex $chain $::smiles_parser::chain_chaincount]] 0]
    
    return [::smiles_parser::_flatten $chain]
}
    
proc ::smiles_parser::_flatten {chain} {
    if {[lindex $chain 0] != "chain" } { error "flatten_atoms has to be run on a chain" }
    
    lset ::smiles_parser::cur_flatten_chain [lindex $chain $::smiles_parser::chain_id] $chain
    
    foreach e [lindex $chain $::smiles_parser::chain_list] {
        if {[lindex $e 0] == "atom"} {
            lset ::smiles_parser::cur_flatten_atom [lindex $e $::smiles_parser::atom_id] $e
            foreach branch [lindex $e $::smiles_parser::atom_branches] {
                ::smiles_parser::_flatten $branch
            }
            foreach ringbond [lindex $e $::smiles_parser::atom_ringbonds] {
                set rb_bond [lindex $ringbond $::smiles_parser::ringbond_bond]
                lset ::smiles_parser::cur_flatten_bond [lindex $rb_bond $::smiles_parser::bond_id] $rb_bond
            }
        } elseif {[lindex $e 0] == "bond"} {
            lset ::smiles_parser::cur_flatten_bond [lindex $e $::smiles_parser::bond_id] $e
        } else {
            error "Unexpected object"
        }
    }
    return [list $::smiles_parser::cur_flatten_atom $::smiles_parser::cur_flatten_bond $::smiles_parser::cur_flatten_chain]
}

proc ::smiles_parser::_unflatten { obj } {
    #This gets data from flattened cache to reconstruct tree
    if {[lindex $obj 0] == "chain" } {
        set chain_id [lindex $obj $::smiles_parser::chain_id]
        set chain [lindex $smiles_parser::cur_flatten_chain $chain_id]
        
        set new_chain_list [list]
        foreach e [lindex $chain $::smiles_parser::chain_list] {
            lappend new_chain_list [::smiles_parser::_unflatten $e]
        }
        lset chain $::smiles_parser::chain_list $new_chain_list
        return $chain
    } elseif {[lindex $obj 0] == "bond" } {
        set bond_id [lindex $obj $::smiles_parser::bond_id]
        set bond [lindex $smiles_parser::cur_flatten_bond $bond_id]
        return $bond
    } elseif {[lindex $obj 0] == "atom" } {
        set atom_id [lindex $obj $::smiles_parser::atom_id]
        set atom [lindex $smiles_parser::cur_flatten_atom $atom_id]
        
        set new_atom_ringbonds [list]
        foreach e [lindex $atom $::smiles_parser::atom_ringbonds] {
            lappend new_atom_ringbonds [::smiles_parser::_unflatten $e]
        }
        lset atom $::smiles_parser::atom_ringbonds $new_atom_ringbonds
        
        set new_atom_branches [list]
        foreach e [lindex $atom $::smiles_parser::atom_branches] {
            lappend new_atom_branches [::smiles_parser::_unflatten $e]
        }
        lset atom $::smiles_parser::atom_branches $new_atom_branches
        
        return $atom
    } elseif {[lindex $obj 0] == "ringbond" } {
        set ringbond $obj
        set bond [lindex $obj $smiles_parser::ringbond_bond]
        lset ringbond $smiles_parser::ringbond_bond [::smiles_parser::_unflatten $bond]
        return $ringbond
    }   else {
        error "Unexpected object"
    }
}


proc ::smiles_parser::smiles_parse {SMILES} {
    #100 times 0
    set ::smiles_parser::chain_refs [lrepeat 100 0]
    
    set smiles_chain [::smiles_parser::smiles_parse_chain $SMILES 0]

    if {[string length $SMILES] != [lindex $smiles_chain 0]} { error "Didn't parse full SMILES string" }
    
    #Remove length (not useful)
    set smiles_chain [lindex $smiles_chain 1]
    
    set ::smiles_parser::cur_id_atom 0
    set ::smiles_parser::cur_id_bond 0
    set ::smiles_parser::cur_id_chain 0
    set ::smiles_parser::cur_ringbond_links [dict create]
    set smiles_chain [::smiles_parser::_chain_renumber $smiles_chain]
    
    #flatten (to initialize flat notation)
    ::smiles_parser::flatten $smiles_chain
    
    #Find neighbors
    ::smiles_parser::_chain_find_neighbors $smiles_chain
    
    #Unflatten to update the tree
    set smiles_chain [::smiles_parser::_unflatten $smiles_chain]
    
    return $smiles_chain
}
