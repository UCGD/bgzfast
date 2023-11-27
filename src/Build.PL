#!/usr/bin/perl
eval 'exec /usr/bin/env perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

BEGIN {
    if($0 ne 'Build.PL'){
	eval 'require FindBin';
	chdir("$FindBin::Bin");
	exec($^X, 'Build.PL')
    }

    our ($cwd) = split(/\n/, `pwd`);
}

use lib("$cwd/inc/lib");
use lib("$cwd/lib");
use lib("$cwd/../perl/lib");

use strict;
use warnings;
use App::Build;

#configure Build
my $build = App::Build->new(
    module_name => 'BGZFast',
    dist_author => 'Carson Holt <carson.holt@genetics.utah.edu>',
    dist_version_from  => 'lib/BGZFast/bgzf_manipulator.pm',
    dist_abstract => 'Fast manipulation of BGZF formatted files',
    license  => 'artistic_2',
    auto_configure_requires => 1, #Module::Build adds self to config_requires
    include_dirs => ['lib/Bio/DB/', 'lib/BGZFast/'], #for included c headers
    requires => {
	'perl'                  => '5.8.0',
	'Bio::Root::Version'    => '1.006923',
	'Inline::C'             => '0.53',
	'forks'                 => '0.34',
	'forks::shared'         => '0.34',
	'Perl::Unsafe::Signals' => '0.02',
	'Perl::Destruct::Level' => '0.02',
	'Storable'              => '0.607',
	'IPC::Semaphore'        => '0',
	'Scalar::Util'          => '0',
	#'Other::Module' => '>= 1.2, != 1.5, < 2.0',
    },
    build_requires => {
    },
    config_requires => {
    },
    recommends => {
    },
    inline_modules => { #additional modules to compile
	'Bio::DB::FAI' => '0.01',
	'BGZFast::tbx_extra' => '0.01',
	#'Other::Module' => '0.01',
    }
);

#ask for MPI
$build->config_hts();

$build->create_build_script;
$build->log_info("\n\nThe script 'Build' has been created for you to finish installing BGZFast.\n");
$build->dispatch('status');
