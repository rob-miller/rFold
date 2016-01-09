#!/usr/bin/perl -w
use parseRtmDssp;

my $cmd="/Users/rob/proj/fold/dssp/rtm-dssp-2.2.1/mkdssp /Users/rob/Desktop/4mqs.pdb";
print "$cmd\n";

open DSSP,"$cmd |" or die "unable to open dssp pipe: $!\n";
my %dsspH = %{ parseRtmDssp(\*DSSP) };

if (exists $dsspH{'pdbid'}) {
  print 'pdbID: ' . $dsspH{'pdbid'} ."\n";
}

my @resArr = @{ $dsspH{'resArr'} };
print "resarr has ". $#resArr ." cardinality\n";

foreach my $r (@resArr) {
  my $dsspRec = %{$r};
  printf "%4d   %4.1f %4.1f %4.1f    %4.1f %4.1f %4.1f \n", $r->{'resNdx'},$r->{'cax'},$r->{'cay'},$r->{'caz'},$r->{'cbx'},$r->{'cby'},$r->{'cbz'};
}

print "done\n";
