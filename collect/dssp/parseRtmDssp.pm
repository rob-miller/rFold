package parseRtmDssp;
require(Exporter);
@ISA = qw(Exporter);

@EXPORT = qw( 
	      parseRtmDssp
);

my $debug=0;

sub parseRtmDssp {
  my $FH = shift;
  my $reading=0;
  my %rsltH=();
  my @rsltA=();
  print "\n" if ($debug);
  while (<$FH>) {
    #print "$_";
    chomp;
    if (/  \#  RESIDUE AA STRUCTURE/) {
      $reading=1;
      print $_ ."\n" if ($debug);
    } elsif ($reading) {
      #  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO   KAPPA  ALPHA   PHI    PSI    OMG    X-CA   Y-CA   Z-CA   X-CB   Y-CB   Z-CB 
      #    3   21 A F  T >> S+     0   0  167      1,-0.2     4,-1.4     2,-0.1     3,-1.1   0.838   95.6   67.3  -70.3  -34.4 -175.6  -29.8  -21.1   56.4  -31.3  -20.7   55.9
      #0123456789012345678901234567890123456789

      my $inline = $_;
      my $posStr = substr($inline,0,14);
      my $strucStr = substr($inline,16,9);
      my $bpStr1 = substr($inline,25,4);
      my $bpStr2 = substr($inline,29,5);
      my $valStr = substr($inline,34);
      if ($debug) {
	print "$inline\n";
	print "posStr=\.$posStr\.\nstrucStr=\.$strucStr\.\nbpStr1=\.$bpStr1\.\nbpStr2=\.$bpStr2\.\nvalStr=\.$valStr\.\n";
      }

      my %dsspRec=();
      if ($posStr =~ /^\s*(\d+)\s+(\-?\d+\w?)\s*(\w)\s+(\w)\s*$/) {
	$dsspRec{'resNdx'} = $1;
	$dsspRec{'resNum'} = $2;
	$dsspRec{'chain'}  = $3;
	$dsspRec{'residue'}= $4;

	print "ndx: $1  num: $2 chain: $3 res: $4\n" if ($debug);
      } elsif ($posStr =~ /^\s*\d+\s+\!\*?\s*$/) {  # chain break
	next;
      } else {
	die "fail to parse posStr: $posStr => _" . (join '.',(split //,$posStr)) ."_\n";
      }

      $dsspRec{'structure'} = $strucStr;
      
      if ($bpStr1 =~/\s*(\w+)\s*$/) {
	$dsspRec{'bp1'} = $1;
      }
      if ($bpStr2 =~/\s*(\w+)\s*$/) {
	$dsspRec{'bp2'} = $1;
      }

      if ($valStr =~ /^\s*(\d+)\s+(\S+)\,\s*(\S+)\s+(\S+)\,\s*(\S+)\s+(\S+)\,\s*(\S+)\s+(\S+)\,\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$/) {
	
	$dsspRec{'acc'} = $1;

	$dsspRec{'nho1offset'} = $2;
	$dsspRec{'nho1energy'} = $3;
	$dsspRec{'ohn1offset'} = $4;
	$dsspRec{'ohn1energy'} = $5;

	$dsspRec{'nho2offset'} = $6;
	$dsspRec{'nho2energy'} = $7;
	$dsspRec{'ohn2offset'} = $8;
	$dsspRec{'ohn2energy'} = $9;

	$dsspRec{'tco'}   = $10;
	$dsspRec{'kappa'} = $11;
	$dsspRec{'alpha'} = $12;

	$dsspRec{'phi'} = $13;
	$dsspRec{'psi'} = $14;
	$dsspRec{'omg'} = $15;

	$dsspRec{'cax'} = $16;
	$dsspRec{'cay'} = $17;
	$dsspRec{'caz'} = $18;

	$dsspRec{'cbx'} = $19;
	$dsspRec{'cby'} = $20;
	$dsspRec{'cbz'} = $21;

	print "cbz= $21\n" if ($debug);
      } else {
	die "fail to parse valStr: _" . (join '.',(split //,$valStr)) ."_  \nposStr=$posStr\nstrucStr=$strucStr\n$valStr\n";
      }

      #print $dsspRec{'resNdx'} .'       '. $dsspRec{'cbz'} ."\n";
      push @rsltA,\%dsspRec;
    } elsif (/^HEADER /) {
      if (/\s(\S+)\s*\.$/) {
	#print "pdbid= $1\n";
	$rsltH{'pdbid'}= $1;
      }
    }
  }
  
  $rsltH{'resArr'} = \@rsltA;
  return \%rsltH;
}

1;




