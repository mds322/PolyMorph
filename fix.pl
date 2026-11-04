#!/bin/perl
if (/^CONECT/) { print; next; }
if (/^(ATOM  |HETATM)/) {
  my @f = split " ", $_;

  my $serial = 0 + $f[1];
  my $aname  = $f[2];
  my $resid  = 0 + $f[4];
  my $x      = 0 + $f[5];
  my $y      = 0 + $f[6];
  my $z      = 0 + $f[7];
  my $occ    = 0 + ($f[8]  // 1.00);
  my $bfac   = 0 + ($f[9]  // 0.00);
  my $elem   =       ($f[10] // "");

  # atom name cols 13-16 (4 chars): if starts with letter and len<4, leading space
  $aname = " $aname" if length($aname) < 4 && $aname =~ /^[A-Za-z]/;
  $aname = sprintf("%-4.4s", $aname);

  $elem =~ s/[^A-Za-z0-9]//g;
  $elem = sprintf("%2.2s", $elem);

  printf("%-6s%5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s\n",
         "HETATM",$serial,$aname,"BGLC","G",$resid,
         $x,$y,$z,$occ,$bfac,"G1",$elem);
  next;
}
print;
