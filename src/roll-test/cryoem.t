#!/usr/bin/perl -w
# cryoem roll installation test.  Usage:
# cryoem.t [nodetype]
#   where nodetype is one of "Compute", "Dbnode", "Frontend" or "Login"
#   if not specified, the test assumes either Compute or Frontend

use Test::More qw(no_plan);

my $appliance = $#ARGV >= 0 ? $ARGV[0] :
                -d '/export/rocks/install' ? 'Frontend' : 'Compute';
my $installedOnAppliancesPattern = '.';
my @packages = (
  'relion','frealign','eman2'
);
my @relion_packages = (
'relion_autopick', 'relion_autopick_mpi', 'relion_display', 'relion_find_tiltpairs', 'relion_image_handler', 'relion_manualpick', 'relion_mask_create', 'relion_particle_polish', 'relion_particle_polish_mpi', 'relion_particle_sort', 'relion_particle_sort_mpi','relion_postprocess','relion_preprocess','relion_preprocess_mpi','relion_project','relion_reconstruct','relion_refine','relion_refine_mpi','relion_run_ctffind','relion_run_ctffind_mpi','relion_stack_create','relion_star_compare','relion_tiltpair_plot');

my @frealign_packages = (
'frealign_v9_mp.exe'
);

my @eman2_packages = (
'fftspeed.py'
);

if($appliance =~ /$installedOnAppliancesPattern/) {
  foreach my $package(@packages) {
    ok(-d "/opt/$package", "$package installed");
  }
} else {
  ok(! $isInstalled, 'eman2 not installed');
}

$packageHome = '/opt/relion';
SKIP: {
  skip 'relion not installed', 1 if ! -d $packageHome;
  foreach my $package(@relion_packages) {
     $output = `module load relion; $package  2>&1`;
     $output = lc($output);
     like($output, qr/options/, "$package works");
  }
}

$packageHome = '/opt/frealign';
SKIP: {
  skip 'frealign not installed', 1 if ! -d $packageHome;
  foreach my $package(@frealign_packages) {
     $output = `module load frealign; $package </dev/null  2>&1`;
     like($output, qr/Use is subject to Janelia Farm Research Campus Software Copyright 1.1/, "$package works");
  }
}

$packageHome = '/opt/eman2';
SKIP: {
  skip 'eman2 not installed', 1 if ! -d $packageHome;
  foreach my $package(@eman2_packages) {
     $output = `module load eman2; python $packageHome/examples/$package </dev/null  2>&1`;
     like($output, qr/ALL TESTS RAN/, "eman2 test $package works");
  }
}


foreach my $package(@packages) {
      `/bin/ls /opt/modulefiles/applications/$package/[0-9]* 2>&1`;
      ok($? == 0, "$package module installed");
     `/bin/ls /opt/modulefiles/applications/$package/.version.[0-9]* 2>&1`;
       ok($? == 0, "$package version module installed");
      ok(-l "/opt/modulefiles/applications/$package/.version",
      "$package  version module link created");
}
