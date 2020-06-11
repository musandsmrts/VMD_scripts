proc oxygenate {mol0 name} {
  # Randomly selects water molecules and renames them so that PSFGEN mutates
  # them to dioxygen molecules.
  # Also creates a new system with the to-be-mutated molecules removed from the
  # system.
  # Both PSF and PDB of system are required.
    
    # Calculate volume
    set minmax0 [measure minmax [atomselect $mol0 all]]
    set dims0 [vecsub [lindex $minmax0 1] [lindex $minmax0 0]]
    puts [format "Dimensions | x: %.03f y: %.03f z: %.03f" [lindex $dims0 0] [lindex $dims0 1] [lindex $dims0 2]]
    set vol0 [expr {[lindex $dims0 0]*[lindex $dims0 1]*[lindex $dims0 2]}]
    puts [format "Volume | %.03f Angstrom^3" $vol0]

  # Suggest number of oxygen molecules and the resulting molar concentration
  # # 100 mM => 100*6.022*10^23*10^-3 molecules / liter
  # #        => 6.022*10^22 molecules / liter
  # #        => 6.022*10^22 molecules / 10^3 cm^3
  # #        => 6.022*10^22 molecules / 10^6 mm^3
  # #        => 6.022*10^22 molecules / 10^15 um^3
  # #        => 6.022*10^22 molecules / 10^24 nm^3
  # #        => 6.022*10^22 molecules / 10^27 Ang^3
  # #        => 6.022*10^-5 molecules / Ang^3
    
    set rec0 [expr {$vol0*6.022*10.**-5}]
    puts "How many molecules would you like to add?"
    puts [format "%.02f molecules would be ~100mM" $rec0]
    gets stdin addNum

    # Select waters
    # Make list of water residue numbers
    set wat0 [atomselect $mol0 "water and name OH2"]
    set wat0L [$wat0 get residue]

    # Randomly select some water residues
    set selL ""
    while {[llength $selL]<$addNum} {
        #puts $selL
        puts [llength $selL]
        #puts $wat0L
        set picker [expr {int(floor(rand()*[llength $wat0L]))}]
        lappend selL [lindex $wat0L $picker]
        set wat0L [lreplace $wat0L $picker $picker]
    }
    puts $selL
    #puts $wat0L

    # Save selected waters to PDB
    [atomselect $mol0 all] set beta 0
    set saver0 [atomselect $mol0 [concat "residue " $selL]]
    $saver0 set beta 1
    set saver0 [atomselect $mol0 [concat "name OH2 H1 and residue " $selL]]

    # We reassign the resID because redundancies cause PSFGEN to skip them when
    # we mutate them to O2.
    # This approach to reassigning might fail if the number of residues exceeds
    # 9999
    set resIDLst ""
    for {set i 0} {$i < $addNum} {incr i} {
        lappend resIDLst [expr {$i+1}] [expr {$i+1}]
        #                 ^             ^
	# Each of the two O2 atoms must have their resID reassigned
    }
    #puts $resIDLst

    $saver0 set resid $resIDLst

    $saver0 writepdb 01-oxySeed.pdb

    # Save remaining atoms from system
    set saver1 [atomselect $mol0 "beta 0"]
    #mol fromsels $saver1
    #mol reanalyze top
    #set saver1 [atomselect top all]
    $saver1 writepsf "$name.psf"
    $saver1 writepdb "$name.pdb"
}
