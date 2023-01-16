set nH2O 128
set lO {}
set lH1 {}
set lH2 {}
for {set x 0} {$x <$nH2O} {incr x} {lappend lO [ expr 3*$x+0]}
for {set x 0} {$x <$nH2O} {incr x} {lappend lH1 [ expr 3*$x+1]}
for {set x 0} {$x <$nH2O} {incr x} {lappend lH2 [ expr 3*$x+2]}
puts $lO
puts $lH1
puts $lH2
set selO [atomselect top "index $lO"]
set selH1 [atomselect top "index $lH1"]
set selH2 [atomselect top "index $lH2"]
for {set x 0} {$x <$nH2O} {incr x} {
set id0 [lindex $lO $x];
set id1 [lindex $lH1 $x];
set id2 [lindex $lH2 $x];
topo addbond $id0 $id1;
topo addbond $id0 $id2;
topo addangle $id1 $id0 $id2;
}
$selO set charge -0.8476
$selH1 set charge 0.4238
$selH2 set charge 0.4238
topo writelammpsdata water.data
