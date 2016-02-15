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
  'relion'
);
my @relion_packages = (
'relion_autopick', 'relion_autopick_mpi', 'relion_display', 'relion_find_tiltpairs', 'relion_image_handler', 'relion_manualpick', 'relion_mask_create', 'relion_particle_polish', 'relion_particle_polish_mpi', 'relion_particle_sort', 'relion_particle_sort_mpi','relion_postprocess','relion_preprocess','relion_preprocess_mpi','relion_project','relion_reconstruct','relion_refine','relion_refine_mpi','relion_run_ctffind','relion_run_ctffind_mpi','relion_stack_create','relion_star_compare','relion_tiltpair_plot');

$packageHome = '/opt/relion';
SKIP: {
  skip 'relion not installed', 1 if ! -d $packageHome;
  foreach my $package(@relion_packages) {
     $output = `module load relion; $package  2>&1`;
     $output = lc($output);
     like($output, qr/options/, "$package works");
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
