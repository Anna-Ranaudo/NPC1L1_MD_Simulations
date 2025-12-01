#!/usr/bin/perl
use strict;
use warnings;
use Carp;

$| = 1;

use lib $ENV{HOME};
use Data::Dumper;

INIT: {
    my $splash = <<END
********************************************************************************
PDB4AMBER
release 15.4.a

Copyright (c) 2014, Dario CORRADA <dario.corrada\@gmail.com>

This work is licensed under a Creative Commons
Attribution-NonCommercial-ShareAlike 3.0 Unported License.
********************************************************************************
END
    ;
    print $splash;
    my $usage = <<END

This script convert a generic PDB file in a format compliant to AMBER software.
Only the lines with flag "ATOM" will be parsed by this script

SYNOPSIS

    \$ PDB4AMBER.pl <protein.pdb>
END
    ;
    unless ($ARGV[0]) {
        print $usage;
        goto FINE;
    }
}

CORE: {

    my $pdb_content = [ ];
    my $hetatm_content = [ ];
    my $header;

    open (PDBIN, '<' . $ARGV[0]) or croak "E- unable to read [$ARGV[0]]\n\t";
    while (my $newline = <PDBIN>) {
        chomp $newline;
        if ($newline =~ /^(TITLE|MODEL|EXPDTA|REMARK)/) {
            $header .= "$newline\n";
        } elsif ($newline =~ /^ATOM/) {
            my @splitted = unpack('Z6Z5Z1Z4Z1Z3Z1Z1Z4Z1Z3Z8Z8Z8Z6Z6Z4Z2Z2', $newline);
            push(@{$pdb_content}, [ @splitted ]);
        } elsif ($newline =~ /^HETATM/) {
            push(@{$hetatm_content}, $newline);
        } else {
            next;
        }
    }
    close PDBIN;

    print "\nSearching terminal residues";
    my @Nterms; my @Cterms;
    my $chain = 'null';
    for (my $i = 0; $i < scalar @{$pdb_content}; $i++) {
        if ($pdb_content->[$i][7] eq $chain) {
            next;
        } else {
            push(@Nterms, $pdb_content->[$i][7] . $pdb_content->[$i][8]);
            push(@Cterms, $pdb_content->[$i-1][7] . $pdb_content->[$i-1][8]) if $i > 0;
            $chain = $pdb_content->[$i][7];
        }
    }

    print "\n\nFixing histidines";
    my $is_protonated;
    foreach my $newatm (@{$pdb_content}) {
        if ($newatm->[3] =~ /^\sH/) {
            $is_protonated = 1;
            last;
        }
    }
    if ($is_protonated) {
        my %histidine_list;
        for (my $i = 0; $i < scalar @{$pdb_content}; $i++) {
            if ($pdb_content->[$i][5] =~ /^HIS$/) {
                my $id_tag = $pdb_content->[$i][7] . $pdb_content->[$i][8];
                if ($pdb_content->[$i][3] =~ /HD1/) {
                    $histidine_list{$id_tag} .= 'HD1';
                } elsif ($pdb_content->[$i][3] =~ /HE2/) {
                    $histidine_list{$id_tag} .= 'HE2';
                }
            }
        }
        my $newtag = 'null';
        for (my $i = 0; $i < scalar @{$pdb_content}; $i++) {
            if ($pdb_content->[$i][5] =~ /^HIS$/) {
                my $id_tag = $pdb_content->[$i][7] . $pdb_content->[$i][8];
                if (exists $histidine_list{$id_tag}) {
                    if ($histidine_list{$id_tag} =~ /HD1HE2/) {
                        $pdb_content->[$i][5] = 'HIP';
                        printf("\n\tHIS %s -> HIP", $id_tag) unless ($id_tag eq $newtag);
                    } elsif ($histidine_list{$id_tag} =~ /HE2/) {
                        $pdb_content->[$i][5] = 'HIE';
                        printf("\n\tHIS %s -> HIE", $id_tag) unless ($id_tag eq $newtag);
                    } else {
                        $pdb_content->[$i][5] = 'HID';
                        printf("\n\tHIS %s -> HID", $id_tag) unless ($id_tag eq $newtag);
                    }
                    $newtag = $id_tag;
                }
            }
        }
    } else {
        print "\n\tW- all histidines will be treated as HIE\n";
        for (my $i = 0; $i < scalar @{$pdb_content}; $i++) {
            if ($pdb_content->[$i][5] =~ /^HIS$/) {
                $pdb_content->[$i][5] = 'HIE';
            }
        }
    }

    print "\n\nFixing isoleucines";
    for (my $i = 0; $i < scalar @{$pdb_content}; $i++) {
        if ($pdb_content->[$i][5] =~ /^ILE$/ && $pdb_content->[$i][3] =~ /^ CD $/) {
            $pdb_content->[$i][3] = ' CD1';
        }
    }

    print "\n\nChecking disulfides";
    my $SG = [ ];
    for (my $i = 0; $i < scalar @{$pdb_content}; $i++) {
        if ($pdb_content->[$i][5] =~ /^CYS$/ && $pdb_content->[$i][3] =~ /^ SG $/) {
            push(@{$SG}, $pdb_content->[$i]);
        }
    }
    my %bridge;
    for (my $i = 0; $i < scalar @{$SG}; $i++) {
        for (my $j = $i + 1; $j < scalar @{$SG}; $j++) {
            my $xa = $SG->[$i][11]; $xa =~ s/\s//g; my $xb = $SG->[$j][11]; $xb =~ s/\s//g;
            my $ya = $SG->[$i][12]; $ya =~ s/\s//g; my $yb = $SG->[$j][12]; $yb =~ s/\s//g;
            my $za = $SG->[$i][13]; $za =~ s/\s//g; my $zb = $SG->[$j][13]; $zb =~ s/\s//g;
            my $distance = ($xa - $xb)**2 + ($ya - $yb)**2 + ($za - $zb)**2;
            $distance = sqrt $distance;
            if ($distance < 2.1) {
                my $key1 = $SG->[$i][7].$SG->[$i][8];
                my $key2 = $SG->[$j][7].$SG->[$j][8];
                $bridge{$key1} = 1 unless (exists $bridge{$key1});
                $bridge{$key2} = 1 unless (exists $bridge{$key2});
            }
        }
    }
    my $newbridge = 'null';
    for (my $i = 0; $i < scalar @{$pdb_content}; $i++) {
        my $pattern = $pdb_content->[$i][7] . $pdb_content->[$i][8];
        if (exists $bridge{$pattern}) {
            $pdb_content->[$i][5] = 'CYX';
            printf("\n\tCYS %s -> CYX", $pattern) unless ($pattern eq $newbridge);
        }
        $newbridge = $pattern;
    }

    print "\n\nC-term fix";
    for (my $i = 0; $i < scalar @{$pdb_content}; $i++) {
        my $pattern = $pdb_content->[$i][7] . $pdb_content->[$i][8];
        foreach my $ref (@Cterms) {
            if ($pattern =~ /$ref/ && $pdb_content->[$i][3] =~ /^ (O1 |OC1)$/) {
                $pdb_content->[$i][3] = ' O  ';
            } elsif ($pattern =~ /$ref/ && $pdb_content->[$i][3] =~ /^ (O2 |OC2)$/) {
                $pdb_content->[$i][3] = ' OXT';
            }
        }
    }

    print "\n\nRemoving hydrogens";
    my $pdb_noH = [ ];
    while (my $newline = shift @{$pdb_content}) {
        if ($newline->[3] =~ /^H/ || $newline->[3] =~ /^.H/) {
            next;
        } else {
            push(@{$pdb_noH}, $newline);
        }
    }

    print "\n\nWriting output";
    my $filename = $ARGV[0];
    $filename =~ s/\.pdb$/.amber.pdb/;
    open (PDBOUT, ">", $filename);
    print PDBOUT "REMARK 888 COMPLIANT WITH AMBER FORMAT\n";
    print PDBOUT "REMARK 888 WRITTEN BY PDB4AMBER.pl\n";

    foreach my $newline (@{$pdb_noH}) {
        print PDBOUT join('', @$newline), "\n";
    }
    foreach my $newline (@{$hetatm_content}) {
        print PDBOUT $newline, "\n";
    }
    print PDBOUT "ENDMDL\n";
    close PDBOUT;

    print "\n\nOutput PDB written to <$filename>";
}

FINE: {
    print "\n\n*** REBMA4BDP ***\n";
    exit;
}

