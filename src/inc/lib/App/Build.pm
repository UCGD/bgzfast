#------------------------------------------------------------------------
#----                           App::Build                           ----
#------------------------------------------------------------------------
package App::Build;

use strict;
use warnings;
use Cwd qw();

BEGIN {
    my $Bundled_MB = 0.4205;  #version included in my distribution

    # Find out what version of Module::Build is installed or fail quietly.
    # This should be cross-platform.
    my $Installed_MB = `$^X -e "eval q{require Module::Build; print Module::Build->VERSION} or exit 1"`;
    chomp $Installed_MB;
    $Installed_MB = 0 if $?;

    # Use the bundled copy of Module::Build if it's newer than the installed.
    if ($Bundled_MB > $Installed_MB){
	my $base;
	eval 'require FindBin; $base = $FindBin::Bin';
	$base ||= Cwd::cwd;
	my $I = "$base/inc/bundle/Module-Build/lib";
	unshift(@INC, $I) if($INC[0] ne $I);
	my @PERL5LIB = split(/:/, $ENV{PERL5LIB}) if($ENV{PERL5LIB});
	unshift(@PERL5LIB, $I);
	$ENV{PERL5LIB} = join(':', @PERL5LIB);
    }
    require Module::Build;
}

use POSIX;
use Config;
use DynaLoader;
use File::Copy;
use File::Path;
use File::Basename;
use IPC::Open3;

use base qw(Inline::Build);
__PACKAGE__->add_property( 'exe_requires' );
__PACKAGE__->add_property( 'lib_requires' );

#------------------------------------------------------------------------
#--------------------------------- METHODS ------------------------------
#------------------------------------------------------------------------
sub new {
    my $class = shift @_;
    my $self = $class->SUPER::new(my %args = @_);

    #override install location
    $self->install_base($self->base_dir.'/..');
    $self->install_base_relpaths('exe'  => 'exe');
    $self->install_base_relpaths('data' => 'data');
    $self->install_base_relpaths('lib' => 'lib');
    $self->install_base_relpaths('arch' => 'lib');
    $self->install_base_relpaths('libdoc' => undef); #no libdocs for apps

    #set INC properly to find temporary cpan installs
    unshift(@INC, $self->base_dir().'/inc/perl/lib');

    #performs a check for eternal algorithm dependencies
    $self->exe_requires({}) if(!$self->exe_requires);
    $self->lib_requires({}) if(!$self->lib_requires);
    $self->check_exes;
    $self->check_libs;

    return $self;
}

sub resume {
    my $class = shift @_;
    my $self = $class->SUPER::resume(@_);

    #override install location (needed because of local::lib)
    $self->install_base($self->base_dir.'/..');

    return $self;
}

#returns hts lib
sub config_hts {
    my $self = shift;

    my $base = $self->base_dir;
    my $ebase = $self->install_destination('exe');
    my @exes = grep {/(^|[\/])tabix$/} (<$base/../exe/*/*>, <$base/../exe/*/bin/*>);
    my $tabix = "$ebase/htslib/bin/tabix" if(-f "$ebase/htslib/bin/tabix");
    ($tabix) = grep {/(^|[\/])tabix$/} (<$base/../exe/*/*>, <$base/../exe/*/bin/*>) if(!$tabix);
    $tabix = $self->config('tabix') if(! $tabix && $self->config('tabix') && $self->config('tabix') =~ /(^|[\/])tabix$/);
    ($tabix) = App::Which::where('tabix') if(!$tabix || ! -f $tabix);

    my $htsdir = $tabix || '';
    $htsdir =~ s/\/[^\/]+$//;

    #directories to search for mpi.h
    my @includes = (<$htsdir/include>,
		    <$htsdir/../include>,
		    <$ebase/*/include>,
		    </usr/include>,
		    </usr/hts*/include>,
		    </usr/local/include>,
		    </usr/local/hts*/include>,
		    </usr/lib/include>,
		    </usr/lib/hts*/include>,
		    </usr/local/lib/include>,
		    </usr/local/lib/hts*/include>);

    my ($HTS_INC) = grep {-f "$_/htslib/hts.h"} @includes;
    $HTS_INC //= '';
    $HTS_INC = $self->prompt("\nPlease specify the path to the directory containing 'htslib/hts.h':", $HTS_INC);
    while(! -f "$HTS_INC/htslib/hts.h"){
	$HTS_INC = $self->prompt("\nCannot find 'htslib/hts.h'\n.".
				"Please specify the containing directory path (leave blank to cancel):", '');
	return if(! $HTS_INC);
    }

    #directories to search for mpi.h
    my @libs = (<$HTS_INC/../lib>,
       	        <$htsdir/lib>,
       	        <$htsdir/../lib>,
		<$ebase/*/lib>,
		</usr/lib>,
		</usr/lib/hts*>,
		</usr/hts*/lib>,
		</usr/local/lib>,
		</usr/local/lib/hts*>,
		</usr/local/hts*/lib>);

    my ($HTS_LIB) = grep {-f "$_/libhts.so"} @libs;
    $HTS_LIB //= '';
    $HTS_LIB = $self->prompt("\nPlease specify the path to the directory containing 'libhts.so':", $HTS_LIB);
    while(! -f "$HTS_LIB/libhts.so"){
	$HTS_LIB = $self->prompt("\nCannot find 'libhts.so'\n.".
				"Please specify the containing directory path (leave blank to cancel):", '');
	return if(! $HTS_LIB);
    }
    
    $self->config_data(HTS_INC => $HTS_INC);
    $self->config_data(HTS_LIB => $HTS_LIB);
    $self->add_lib_requires(hts => "$HTS_LIB/libhts.so");
}

sub config_exe_loc {
    my $self = shift;
    my $tag = shift;

    my $base = $self->base_dir;
    my $ebase = $self->install_destination('exe');
    my ($exe) = grep {/(^|[\/])$tag$/ && -f $_} (<$base/../exe/*/*>, <$base/../exe/*/bin/*>);
    ($exe) = App::Which::where($tag) if(!$exe || ! -f $exe);
    $exe = $self->prompt("\nPlease specify the path to '$tag' on your system:", $exe);
    if(!$exe){
	print "Skipping '$tag'...\n";
	return $exe;
    }
    while($exe !~ /(^|[\/])$tag$/ || ! -f $exe){
        $exe = $self->prompt("\nCannot find '$tag'.\n".
			     "Please specify the path (leave blank to skip):", '');
        return if(! $exe);
    }

    $self->config_data($tag => $exe);
    return $exe;
}



#add a Module to the requires list
sub add_requires {
    my $self = shift;
    return unless(@_);

    if (@_ > 1 && @_ % 2 == 0){
	for (my $i = 0; $i < @_; $i += 2){
	    $self->{properties}{requires}{$_[$i]} = $_[$i+1]
	}
    }
    elsif(ref $_[0] eq 'HASH'){
	map {$self->{properties}{requires}{$_} = $_[0]->{$_} } keys %{$_[0]}
    }
    else{
	map {$self->{properties}{requires}{$_} = 0} @_;
    }
}

#add a Module to the build_requires list
sub add_build_requires {
    my $self = shift;
    return unless(@_);

    if (@_ > 1 && @_ % 2 == 0){
	for (my $i = 0; $i < @_; $i += 2){
	    $self->{properties}{build_requires}{$_[$i]} = $_[$i+1]
	}
    }
    elsif(ref $_[0] eq 'HASH'){
	map {$self->{properties}{build_requires}{$_} = $_[0]->{$_} } keys %{$_[0]}
    }
    else{
	map {$self->{properties}{build_requires}{$_} = 0} @_;
    }
}

#add a Module to the exe_requires list
sub add_exe_requires {
    my $self = shift;
    return unless(@_);

    if (@_ > 1 && @_ % 2 == 0){
	for (my $i = 0; $i < @_; $i += 2){
	    $self->{properties}{exe_requires}{$_[$i]} = $_[$i+1]
	}
    }
    elsif(ref $_[0] eq 'HASH'){
	map {$self->{properties}{exe_requires}{$_} = $_[0]->{$_} } keys %{$_[0]}
    }
    else{
	map {$self->{properties}{exe_requires}{$_} = 0} @_;
    }
}

#add a Library to the lib_requires list
sub add_lib_requires {
    my $self = shift;
    return unless(@_);

    if (@_ > 1 && @_ % 2 == 0){
	for (my $i = 0; $i < @_; $i += 2){
	    $self->{properties}{lib_requires}{$_[$i]} = $_[$i+1]
	}
    }
    elsif(ref $_[0] eq 'HASH'){
	map {$self->{properties}{lib_requires}{$_} = $_[0]->{$_} } keys %{$_[0]}
    }
    else{
	map {$self->{properties}{lib_requires}{$_} = 0} @_;
    }
}

#replaces Module::Build's config method
sub config {
    my $self = shift;
    
    #hack to get around bug in Module::Build 
    #otherwise it will ignore changes to config 'cc' and 'ld'
    $self->{stash}->{_cbuilder} = undef;

    return $self->SUPER::config(@_);
}

#override default build
sub ACTION_code {
    my $self = shift;

    my @perl = map {keys %{$_->{requires}}} $self->prereq_failures();
    my @exes = map {keys %{$_->{exe_requires}}} $self->exe_failures();
    my @libs = map {keys %{$_->{lib_requires}}} $self->lib_failures();

    if(@perl || @exes || @libs){
	print "\nERROR: Cannot '".$self->invoked_action."', missing prerequisites\n";
	exit;
    }

    #add @INC to PERL5LIB before compiling just in case
    local $ENV{PERL5LIB} = join(':', @INC);

    # All installable stuff gets created in blib/ .
    # Create blib/arch to keep blib.pm happy
    my $blib = $self->blib;
    $self->add_to_cleanup($blib);
    File::Path::mkpath( File::Spec->catdir($blib, 'arch') );
    File::Path::mkpath( File::Spec->catdir($blib, 'script') );
    File::Path::mkpath( File::Spec->catdir($blib, 'lib') );

    #compile MPI Communicator module
    if($self->feature('mpi_support')){
	my $inline_modules = $self->inline_modules;
	my $module = 'Parallel::Architect::Comm::MPI';
	my $version = $inline_modules->{$module};

	my $perl   = $self->perl;
	my $arch   = $self->blib.'/arch';

	my %inline;

	#adaptation of Inline::Build command line
	(my $inl_file = "$module.inl") =~ s/\:\:/-/g;
	if(!-f $inl_file){
	    my $text = '';
	    foreach my $key (sort keys %inline){
		my $value = $inline{$key};
		$text .= qq($key => "$value", );
	    }
	    my $command = join(' ', $perl => ("-Mblib",
					      "-MInline=NOISY,_INSTALL_",
					      "-M$module",
					      '-e' => ('\''.$module.'->bind('.$text.');'.
						       'my %A = (modinlname => "'.$inl_file.'", module => "'.$module.'");'.
						       'my %S = (API => \%A);'.
						       'Inline::satisfy_makefile_dep(\%S);\''),
					      "$version",
					      "$arch"));
	    print STDERR $command."\n";
	    my $stat = $self->do_system($command);
	    
	    die "ERROR: $module failed to compile\n" unless($stat);
	}
	$self->add_to_cleanup($inl_file);
	$self->add_to_cleanup('_Inline');

	#write ConfigData.pm file for MPI module
	my $notes_name = $module.'::ConfigData';
	my $notes_pm = File::Spec->catfile($self->blib, 'lib', split(/::/, "$notes_name.pm"));
	File::Path::mkpath(File::Basename::dirname($notes_pm));

	if(!-f $notes_pm){
	    Module::Build::Notes->write_config_data(file => $notes_pm,
						    module => $module,
						    config_module => $notes_name,
						    config_data => {HTS_INC => $self->config_data('HTS_INC'),
								    HTS_LIB => $self->config_data('HTS_LIB')},
						    auto_features => {});
	}
    }
    $self->notes(code_done => 1);

    #compile the rest
    $self->SUPER::ACTION_code(@_);
}

#replaces Module::Build's ACTION_installdeps
#so prereqs install locally inside of .../perl/lib
sub ACTION_installdeps{
    my $self = shift;

    my $prereq = $self->prereq_failures();
    if($prereq && ($prereq->{requires} || $prereq->{recommends})){
	my @perl = map {keys %{$_->{requires}}} $self->prereq_failures();
	my $access = 1 if(-w $Config{installsitelib} && -w $Config{installarchlib});
	if(!$access){
	    my ($root_id, $grp_id) = (getpwnam('root'))[2,3];
	    $access = 1 if($< == $root_id || $> == $root_id);
	}

	my $local;
	if(! $access){
	    $local = $self->y_n("You do not appear to have write access to globally install\n".
				"missing Modules. However, if you are using local::lib or\n".
				"have manually edited CPAN configuration to do local installs\n".
				"then this message may be in error. If you want, I can try\n".
				"and install all modules locally (i.e. only for this package)\n".
				"in the .../<package>/perl/lib directory, or you can run\n".
				"'./Build installdeps' as root or using sudo and try again.\n".
				"Do want to try and build a local installation?", 'N');

	}	

	if(!$access && !$local){
	    print "\n\nWARNING: You do not appear to have write access to install missing\n".
		  "Modules. Please run './Build installdeps' as root or using sudo.\n\n";

	    my $go = $self->y_n("Do you want to continue anyway?", 'N');

	    exit(0) unless($go);
	}

	my %errors;
	foreach my $m (keys %{$prereq->{build_requires}}){    
	    $self->cpan_install($m, $local);
	}
	foreach my $m (keys %{$prereq->{requires}}){    
	    $self->cpan_install($m, $local);
	}

	print "Checking optional dependencies:\n";
	foreach my $m (keys %{$prereq->{recommends}}){    
	    my $go = $self->y_n("Install $m?", 'Y');
	    if($go){
		$errors{$m}++ unless($self->cpan_install($m, $local));
	    }
	}
	
	print "\nRechecking dependencies to see if installation was successful\n";
	$self->check_prereq;
	
	if($self->prereq_failures() && keys %errors){
	    if(! keys %{$prereq->{requires}}){
		print "\nWARNING: You have all required dependencies, but some optional components\n".
		    "failed to install.\n";
	    }
	    else{
		my ($usr_id) = (getpwnam('root'))[2];
		print "\nWARNING: Installation failed (please review any previous errors).\n";
		print "Try installing the missing packages as 'root' or using sudo.\n" if($< != $usr_id);
		print "You may need to configure and install these packages manually.\n";
		return 0;
	    }
	}
    }

    $self->status;
}

#these install individual external algorithms
sub ACTION_htslib{ shift->_exe_action('htslib'); }


#runs all the algorithm installs that are missing
sub ACTION_installexes{
    my $self = shift;

    my $exe_failures = $self->exe_failures();

    return if(! $exe_failures || ! $exe_failures->{exe_requires});
    foreach my $name (keys %{$exe_failures->{exe_requires}}){
	next if(!$exe_failures->{exe_requires}{$name});
	$name = lc($name); #all dispatch parameters are lower case
	$self->dispatch($name);
	$exe_failures = $self->exe_failures(); #solve recusive install
    }
    
    $self->check_exes;
}

#prints out a simple configuration status message
sub ACTION_status {
    my $self = shift;
    
    $self->status;
}

#syncronize the maker/src/bin and maker/src/inc/bin directories
#to maker/bin because of user edits to scripts
sub ACTION_sync {
    my $self = shift;

    $self->sync_dirs('lib', 'bin');
}

#update current MAKER from the MAKER repository
sub ACTION_update {
    my $self = shift;
    my $base = $self->base_dir;

    return if(!-d "$base/../.svn" && !-d "$base/../../.svn"); #not a subversion repository

    $self->depends_on('sync');
    my ($s_svn) = $self->svn_w_args('info') =~ /Revision\:\s*(\d+)/;
    print $self->svn_w_args('update');
    my ($f_svn) = $self->svn_w_args('info') =~ /Revision\:\s*(\d+)/;

    #there were changes so re-run install
    if($s_svn != $f_svn){
	my $go = $self->y_n("There were changes made to your source repository. ".
			    "Do you wish to reinstall?", 'Y');
	if($go){
	    $self->depends_on('clean');
	    $self->depends_on('install');
	}
    }

    print "\nSVN STATUS:\n";
    print $self->svn_w_args('status');
}

#commit current MAKER to subversion repository
sub ACTION_commit {
    my $self = shift;
    my @args = @{$self->args->{ARGV}};

    die "ERROR: Failure to provide log message for commit. Use -m 'message'\n"
	if(@args < 2 || $args[0] ne '-m');

    print "Updating repository before commit...\n";
    $self->depends_on('update');
    
    print $self->svn_w_args('commit', @args);
}

#create a versioned release
sub ACTION_release {
    my $self = shift;

    App::Which::which('tar') || die "ERROR: Cannot find tar to build the release\n";

    #update current repository
    print "\nUpdating to most current repository...\n";
    print $self->svn_w_args('update');
    my ($s_svn) = $self->svn_w_args('info') =~ /Revision\:\s*(\d+)/;

    #do prerelease commit
    my $status = $self->svn_w_args('status');
    if($status =~ /^[^\?]/m){
	print "\nPre-release commit of any user changes...\n";
	print $self->svn_w_args('commit', '-m', 'pre-release commit');
	$self->svn_w_args('update');
    }

    #clean and check/update version for release
    $self->dispatch('clean');
    my $ver = $self->check_update_version(); #returns new version
    $self->{properties}->{dist_version} = $ver;

    #build tarball for users to download
    my $base = $self->base_dir;
    (my $ibase = $base) =~ s/\/[^\/]+$//;
    my $tgz = $self->dist_name."-$ver.tgz";
    if(! -f "$base/$tgz"){
	my ($name, $dir) = File::Basename::fileparse($ibase);
        my $exclude = $self->svn_w_args('status');
        $exclude = join("\n", map {"$name/$_"} ($exclude =~ /\?\s+([^\n]+)/g))."\n";
        open(OUT, "> $base/.exclude");
        print OUT $exclude;
	print OUT "$name/src/$tgz\n";
	print OUT "$name/src/.exclude\n";
        close(OUT);

        print "\nBuilding tarball for distribution...\n";
        my $command = "tar -C $dir --exclude \"*~\" --exclude \".svn\" -X $base/.exclude -zcf $base/$tgz $name";
        $self->do_system($command) || unlink($tgz);
        unlink("$base/.exclude");
        die "ERROR: tarball creation failed\n" if(! -f $tgz);
    }

    #there were changes so re-run install (updates version info in scripts)
    my ($f_svn) = $self->svn_w_args('info') =~ /Revision\:\s*(\d+)/;
    if($s_svn != $f_svn){
        print "\nNow reinstalling scripts to reflect version changes...\n";
        $self->depends_on('realclean'); #clean up all old files
        $self->create_build_script; #update stored Build script
        $self->depends_on('install');
    }
}


#checks external algorithms to see if they're present. Anologous to check_prereqs
sub check_exes{
    my $self = shift;

    my $exe_failures = $self->exe_failures();
    return if(! $exe_failures);

    print "Checking external program dependencies...\n";
    while(my $cat = each %$exe_failures){
	my $s = ($cat eq 'exe_requires') ? '!' : '*';
	$cat =~ /exe_(.*)/;

	print "  $1:\n";

	while(my $name = each %{$exe_failures->{$cat}}){
	    print "    $s  ".$exe_failures->{$cat}{$name}{message} ."\n";
	}
    }
    print "\nERRORS/WARNINGS FOUND IN PREREQUISITES.  You may wish to install the programs\n".
	"indicated above before proceeding with this installation.\n".
	"Run 'Build installexes' to install missing prerequisites.\n\n";
}

#checks external libraries to see if they're present. Anologous to check_prereqs
sub check_libs{
    my $self = shift;

    my $lib_failures = $self->lib_failures();
    return if(! $lib_failures);

    print "Checking external library dependencies...\n";
    while(my $cat = each %$lib_failures){
	my $s = ($cat eq 'lib_requires') ? '!' : '*';
	$cat =~ /lib_(.*)/;

	print "  $1:\n";

	while(my $name = each %{$lib_failures->{$cat}}){
	    print "    $s  ".$lib_failures->{$cat}{$name}{message} ."\n";
	}
    }
    print "\nERRORS/WARNINGS FOUND IN PREREQUISITES.  You may wish to install the libraries\n".
	"indicated above before proceeding with this installation\n".
	"\nInstaller cannot do this for you.  You will have to do it manually.\n";
}

#returns missing exes, anologous to prereq_failures
{
my $set_PATH = 0;
sub exe_failures {
    my $self = shift;
    my $other = shift || {};
    my %exes = (%{$self->exe_requires}, %{$other});

    if(! $set_PATH && $self->config_data('PATH')){
	$ENV{PATH} = ($ENV{PATH}) ?
	    $ENV{PATH}.":".$self->config_data('PATH') : $self->config_data('PATH');
	$set_PATH = 1;
    }

    my %exe_failures;
    while (my $name = each %exes){
	my @cmds = map {$_ =~ s/\s+$|^\s+//g; $_} split(/,/, $exes{$name});
	my $base = $self->install_destination('exe');
	my $dest = (-d "$base/$name") ? "$base/$name" : "$base/".lc($name);

	my $loc;
	foreach my $cmd (@cmds){
	    last if(($loc) = grep {$_ && -f $_} $self->config_data($cmd));
	    last if(($loc) = grep {-f $_} ("$dest/$cmd", "$dest/bin/$cmd", App::Which::which($cmd)));
	    last if(($loc) = grep {-f $_ && -x $_} ($cmd));
	}

        if(! $loc){
	    $exe_failures{'exe_requires'}{$name}{have} = '<none>';
	    $exe_failures{'exe_requires'}{$name}{message} = "$name is not installed";
	    $exe_failures{'exe_requires'}{$name}{need} = $exes{$name};
	}
    }

    return (keys %exe_failures) ? \%exe_failures : undef;
}
}

#returns missing libs, anologous to prereq_failures
sub lib_failures {
    my $self = shift;
    my $other = shift || {};
    my %libs = (%{$self->lib_requires}, %{$other});

    my $user_extra = $self->config_data('include_extra');
    push(@DynaLoader::dl_library_path, @{$user_extra}) if($user_extra);
    push(@DynaLoader::dl_library_path, @{$self->include_dirs});
    push(@DynaLoader::dl_library_path, '/usr/X11/include', '/usr/X11/include', '/sw/lib', '/sw/include'); #for libpng

    my %lib_failures;
    while (my $name = each %libs){
	my @cmds = map {$_ =~ s/\s+$|^\s+//g; $_} split(/,/, $libs{$name});

	my $loc;
	foreach my $cmd (@cmds){
	    last if(($loc) = grep {-f $_} (DynaLoader::dl_findfile($cmd), $cmd));
	}

        if(! $loc){
	    $lib_failures{'lib_requires'}{$name}{have} = '<none>';
	    $lib_failures{'lib_requires'}{$name}{message} = "$name is not installed";
	    $lib_failures{'lib_requires'}{$name}{need} = $libs{$name};
	}
    }

    return (keys %lib_failures) ? \%lib_failures : undef;
}


#hidden connection entry between ACTION_??? algorithm install
#checks to see if already installed and then run the install method
sub _exe_action{
    my $self = shift;
    my $label = shift;
    my $script = shift;

    if(! $script && ! $self->exe_requires->{$label}){
	$script = $label;
    }

    my $fail = ($script) ? $self->exe_failures({$label => $script}) : $self->exe_failures();
    my @list = keys %{$fail->{exe_requires}} if($fail->{exe_requires});
    if(grep {/^$label$/i} @list){
	$self->_install_exe($label);
    }
    else{
	my $go = $self->y_n("WARNING: $label was already found on this system.\n".
			    "Do you still want me to install $label for you?", 'N');
	$self->_install_exe($label) if($go);
    }
}

#does actual installation of all external algorithms
sub _install_exe {
    my $self = shift;
    my $exe  = shift;
    my $base = $self->install_destination('exe');
    my $path = "$base/$exe";

    #get OS and architecture
    my %os_ok = (Linux_x86_64  => 1,
		 Linux_i386    => 1,
		 Darwin_i386   => 1,
		 Darwin_x86_64 => 1,
		 src           => 1); 
    my ($OS, $ARC) = (POSIX::uname())[0,4];

    #set all pentium achitectures to use i386 (safest choice)
    if($ARC =~ /^i.86$/){
	$ARC = 'i386';
    }

    ($OS, $ARC) = ('src', '') if(! $os_ok{"$OS\_$ARC"}); #use source code for unknown architectures

    #add fink paths just in case
    if($OS eq 'Darwin'){
	$ENV{C_INCLUDE_PATH} .= ":" if($ENV{C_INCLUDE_PATH});
	$ENV{LIBRARY_PATH} .= ":" if($ENV{LIBRARY_PATH});
	$ENV{C_INCLUDE_PATH} .= "/sw/include:/usr/X11/include";
	$ENV{LIBRARY_PATH} .= "/sw/lib:/usr/X11/lib";
    }

    #get url for exectuable to be installed
    my $data;
    open(LOC, '<', $self->base_dir()."/locations")
	or die "ERROR: Could not open locations to download dependencies\n";
    my $line = <LOC>;
    if($line =~ /^\#\#PACKAGE/){
	$data = join('', <LOC>);
	eval $data;
    }
    close(LOC);

    #prerequisite installation directory
    if(! -d $base){
	mkdir($base) ||
	    die "ERROR could not create directory $base for installing external program dependencies\n";
    }

    #install
    chdir($base);
    my @unlink;
    if($exe eq 'htslib'){
	#htslib
	File::Path::rmtree($path);
	my $file = "$base/$exe.tar.bz2"; #file to save to
        my $url = $data->{$exe}{"$OS\_$ARC"}; #url to blast for OS
	print "Downloading $exe...\n";
        $self->getstore($url, $file) or return $self->fail($exe, $path);
	print "Unpacking $exe tarball...\n";
        $self->extract_archive($file) or return $self->fail($exe, $path);
        push (@unlink, $file);
	my ($dir) = grep {-d $_} <htslib-*>;
	if(-d $dir){ #this is the source code and must be compiled
	    chdir($dir);
	    print "Configuring $exe...\n";
	    $self->do_system("./configure --prefix=$path") or return $self->fail($exe, $path);
	    $self->do_system("make install") or return $self->fail($exe, $path);
	}
        chdir($base);
	return $self->fail($exe, $path) if(! -f "$path/bin/tabix");

	#change back to proper directory or config_data doesn't work
	chdir($self->base_dir());
	$self->config_data(tabix   => "$path/bin/tabix");
	$self->config_data(HTS_INC => "$path/include");
	$self->config_data(HTS_LIB => "$path/lib");
	$self->config_data(include_extra => ["$path/include"]);
    }
    else{
	die "ERROR: No install method defined for $exe in INSTALL::Build::_install_exe.\n";
    }

    #report status
    print "Finished installing $exe.\n";

    #remove all tarballs/etc
    map{unlink($_)} @unlink;

    #change back to proper directory
    chdir($self->base_dir());
}

#fail/cleanup method for installing exes
sub fail {
    my $self = shift;
    my $exe = shift;
    my $path = shift;

    print "\n\nERROR: Failed installing $exe, now cleaning installation path...\n".
	"You may need to install $exe manually.\n\n";

    if(-f $path){
	unlink($path);
    }
    else{
	File::Path::rmtree($path);
    }

    chdir($self->base_dir()); #just in case
}

# install an external module using CPAN prior to testing and installation
# borrowed and modified from BioPerl, has flag for local install
sub cpan_install {
    my ($self, $desired, $local) = @_;

    unless(App::Which::which('make')){
	die "\n\n".
            "ERROR: Cannot find 'make' on your system. If this is a Mac you will need\n".
	    "to install Xcode developer tools before you can continue. This can be\n".
	    "done from the OS X installation disk or downloaded from the App Store.\n".
            "On some systems you must also select 'install command line tools' under\n".
	    "the 'Downloads' section of the 'Xcode >> Preferences' menu\n\n";
    }

    if(! $local){
	my $loc = $self->module_overide($desired);

	if($loc){
	    print "\n\nWARNING: There is another version of the module $desired\n".
		  "installed on this machine that will supercede a globally installed version.\n".
		  "I can only continue by installing a local package-only version. If you want a\n".
		  "global installation of the module, you will have to quit and delete the\n".
		  "offending module.\n".
                  "Location: $loc\n";

	    $local = $self->y_n("Do you want to continue with a local installation?", 'N');

	    die "\nWARNING: You will need to delete $loc before continuing.\n"if(! $local);
	}
    }

    #set up PERL5LIB environmental varable since CPAN doesn't see my 'use lib'
    if($local){
	my $PERL5LIB = $ENV{PERL5LIB} || '';
	$PERL5LIB = $self->base_dir."/../perl/lib:$PERL5LIB";
	$ENV{PERL5LIB} = $PERL5LIB;
    }
    
    #possible temporary install base location for any CPAN requirements
    my $base = $self->base_dir;
    if(-d "$base/inc/perl" && !scalar(grep {$_ eq "$base/inc/perl/lib"} @INC)){
	unshift(@INC, "$base/inc/perl/lib");
	my $PERL5LIB = $ENV{PERL5LIB} || '';
	$PERL5LIB = "$base/inc/perl/lib:$PERL5LIB";
	$ENV{PERL5LIB} = $PERL5LIB;
    }

    #if CPAN is too old, then install a newer one first
    if(!$self->check_installed_status('CPAN', '1.82')->{ok}){
	if(! -d "$base/inc/perl"){
	    mkdir("$base/inc/perl");
	    unshift(@INC, "$base/inc/perl/lib");
	    my $PERL5LIB = $ENV{PERL5LIB} || '';
	    $PERL5LIB = "$base/inc/perl/lib:$PERL5LIB";
	    $ENV{PERL5LIB} = $PERL5LIB;
	}

	#fill out environment variables
	my $perl_mb = $ENV{PERL_MB_OPT}; #backup
	$ENV{PERL_MB_OPT} = "--install_base $base/inc/perl/ --installdirs site".
            " --install_path libdoc=$base/inc/perl/man --install_path bindoc=$base/inc/perl/man".
            " --install_path lib=$base/inc/perl/lib --install_path arch=$base/inc/perl/lib".
            " --install_path bin=$base/inc/perl/lib/bin --install_path script=$base/inc/perl/lib/bin ";;
	my $perl_mm = $ENV{PERL_MM_OPT}; #backup
	$ENV{PERL_MM_OPT} = "DESTDIR=$base/inc/perl/ INSTALLDIRS=site INSTALLSITEMAN1DIR=man INSTALLSITEMAN3DIR=man".
            " INSTALLSITEARCH=lib INSTALLSITELIB=lib INSTALLSITEBIN=lib/bin INSTALLSITESCRIPT=lib/bin";
	my $prefer = $ENV{PERL_AUTOINSTALL_PREFER_CPAN}; #backup
	$ENV{PERL_AUTOINSTALL_PREFER_CPAN} = 1;
	my $mm_def = $ENV{PERL_MM_USE_DEFAULT}; #backup
	$ENV{PERL_MM_USE_DEFAULT} = 1;

	#CPAN config from local::lib's Makefile.PL
	my $cpan_config_command =
            'my $done; require ExtUtils::MakeMaker;
             my $orig = ExtUtils::MakeMaker->can("prompt");
             *ExtUtils::MakeMaker::prompt = sub ($;$) {
               if (!$done && $_[0] =~ /manual configuration/) {
                 $done++;
                 return "no";
               }
               return $orig->(@_);
             };
             $CPAN::Config->{prefer_installer} = "EUMM";
             CPAN::Config->load;
             unless ($done || -w $CPAN::Config->{keep_source_where}) {
               my $save = $CPAN::Config->{urllist};
               delete @{$CPAN::Config}{keys %$CPAN::Config};
               $CPAN::Config->{urllist} = $save;
               CPAN::Config->init;
             }';
	my $cpan_command = '';
	$cpan_command .= 'force("notest","install","ExtUtils::MakeMaker"); '
	    if(!$self->check_installed_status('ExtUtils::MakeMaker', '6.31')->{ok});
	$cpan_command .= 'force("notest","install","ExtUtils::Install"); '
	    if(!$self->check_installed_status('ExtUtils::Install', '1.43')->{ok});
	$cpan_command .= 'force("notest","install","CPAN"); ';

	#run CPAN via system call
	system($^X, '-MCPAN', '-e', $cpan_config_command);
	system($^X, '-MCPAN', '-e', $cpan_command);

	$ENV{PERL_MB_OPT} = $perl_mb; #restore
	$ENV{PERL_MM_OPT} = $perl_mm; #restore
	$ENV{PERL_MM_USE_DEFAULT} = $mm_def; #restore
	$ENV{PERL_AUTOINSTALL_PREFER_CPAN} = $prefer; #restore
    }

    # Here we use CPAN to actually install the desired module
    require CPAN;
    import MyModule;

    # Save this because CPAN will chdir all over the place.
    my $cwd = getcwd();

    #set up a non-global local module library for package
    my %bak;
    if($local){
	CPAN::HandleConfig->load;
	%bak = (makepl_arg => $CPAN::Config->{makepl_arg},
		mbuildpl_arg => $CPAN::Config->{mbuildpl_arg});
	$CPAN::Config->{makepl_arg} = "DESTDIR=$base/../perl/ INSTALLDIRS=site INSTALLSITEMAN1DIR=man INSTALLSITEMAN3DIR=man".
	    " INSTALLSITEARCH=lib INSTALLSITELIB=lib INSTALLSITEBIN=lib/bin INSTALLSITESCRIPT=lib/bin";
	$CPAN::Config->{mbuildpl_arg} = "--install_base $base/../perl/ --installdirs site".
	    " --install_path libdoc=$base/../perl/man --install_path bindoc=$base/../perl/man".
	    " --install_path lib=$base/../perl/lib --install_path  arch=$base/../perl/lib".
	    " --install_path bin=$base/../perl/lib/bin --install_path script=$base/../perl/lib/bin ";
	$CPAN::Config->{prefs_dir} = "$ENV{HOME}/.cpan/prefs" if(! -w $CPAN::Config->{prefs_dir});
	CPAN::Shell::setup_output();
	CPAN::Index->reload;
    }
    else{
	CPAN::HandleConfig->load;
	%bak = (makepl_arg => $CPAN::Config->{makepl_arg},
		mbuildpl_arg => $CPAN::Config->{mbuildpl_arg});
	$CPAN::Config->{makepl_arg} = "INSTALLDIRS=site";
	$CPAN::Config->{mbuildpl_arg} = "--installdirs site";
	CPAN::Shell::setup_output();
	CPAN::Index->reload;
    }

    #install YAML if needed to avoid other installation issues with prereqs
    CPAN::Shell->force('notest', 'install', 'YAML') if (! $self->check_installed_status('YAML', '0')->{ok});

    #CPAN::Shell->expand("Module", $desired)->cpan_version <= 2.16;
    #CPAN::Shell->install($desired);
    CPAN::Shell->force('notest', 'install', $desired);

    #restore old CPAN settings
    $CPAN::Config->{makepl_arg} = $bak{makepl_arg};
    $CPAN::Config->{mbuildpl_arg} = $bak{mbuildpl_arg};
    CPAN::Shell::setup_output();

    my $ok;
    my $expanded = CPAN::Shell->expand("Module", $desired);
    if ($expanded && $expanded->uptodate) {
	print "$desired installed successfully\n";
	$ok = 1;
    }
    else {
	print "$desired failed to install\n";
	$ok = 0;
    }
    
    chdir $cwd or die "Cannot chdir() back to $cwd: $!";
    return $ok;
}

#untars a package. Tries to use tar first then moves to the perl package untar Module.
sub extract_archive {
    my $self = shift;
    my $file = shift;

    return 0 if(! $file);
    
    if(App::Which::which('tar')){
	my $command;
	my $u = scalar getpwuid($>);
	my $g = scalar getgrgid($));
	if($file =~ /\.gz$|\.tgz$/){
	    $command = "tar -zxm -f $file";
	}
	elsif($file =~ /\.bz2?$|\.tbz2?$/){
	    $command = "tar -jxm -f $file";
	}
	else{
	    $command = "tar -xm -f $file";
	}
	$command .= " --owner $u --group $g" unless((POSIX::uname())[0] =~ /darwin/i);

	return $self->do_system($command); #fast
    }
    else{
	die "ERROR: Archive::Tar required to unpack missing executables.\n".
	    "Try running ./Build installdeps first.\n\n"
	    if(!$self->check_installed_status('Archive::Tar', '0')->{ok});

	eval 'require Archive::Tar';
	return (Archive::Tar->extract_archive($file)) ? 1 : 0; #slow
    }
}

#downloads files from the internet.  Tries to use wget, then curl,
#and finally LWP::Simple
sub getstore {
    my $self = shift;
    my $url = shift;
    my $file = shift;
    my $user = shift;
    my $pass = shift;

    if(App::Which::which('wget')){ #Linux
	my $command = "wget $url -c -O $file --no-check-certificate";
	$command .= " --user $user --password $pass" if(defined($user) && defined($pass));
	return $self->do_system($command); #gives status and can continue partial
    }
    elsif(App::Which::which('curl')){ #Mac
	my $command = "curl --connect-timeout 30 -f -L $url -o $file";
	$command .= " --user $user:$pass" if(defined($user) && defined($pass));
	my $continue = " -C -";

	#gives status and can continue partial
	my $stat = $self->do_system($command . $continue);
	#just redo if continue fails
	$stat = $self->do_system($command) if(! $stat);
	return $stat;
    }
    else{
	die "ERROR: LWP::Simple required to download missing executables\n".
	    "Try running ./Build installdeps first.\n\n"
	    if(!$self->check_installed_status('LWP::Simple', '0')->{ok});

	eval 'require LWP::Simple';
	$url =~ s/^([^\:]\;\/\/)/$1\:\/\/$user\:$pass\@/ if(defined($user) && defined($pass));
	return LWP::Simple::getstore($url, $file); #just gets the file with no features
    }
}



#prints a nice status message for package configuration and install
sub status {
    my $self = shift;

    my @perl = map {keys %{$_->{requires}}} $self->prereq_failures();
    my @exes = map {keys %{$_->{exe_requires}}} $self->exe_failures();
    my @libs = map {keys %{$_->{lib_requires}}} $self->lib_failures();

    my $dist_name = $self->dist_name;
    my $dist_version = $self->dist_version;

    my $mpi = ($self->feature('mpi_support')) ? 'ENABLED' : 'DISABLED';
    my $stat = 'CONFIGURATION OK';
    $stat = 'MISSING PREREQUISITES' if(@perl || @exes || @libs);
    $stat = 'INSTALLED' if($self->notes('install_done'));

    print "\n\n";
    print "==============================================================================\n";
    print "STATUS $dist_name $dist_version\n";
    print "==============================================================================\n";
    print "PERL Dependencies:\t";
    print ((@perl) ? 'MISSING' : 'VERIFIED');
    print"\n";
    print "\t\t  !  ". join("\n\t\t  !  ", @perl) ."\n\n" if(@perl);
    print "External Programs:\t";
    print ((@exes) ? 'MISSING' : 'VERIFIED');
    print "\n";
    print "\t\t  !  ". join("\n\t\t  !  ", @exes) ."\n\n" if(@exes);
    print "External C Libraries:\t";
    print ((@libs) ? 'MISSING' : 'VERIFIED');
    print "\n";
    print "\t\t  !  ". join("\n\t\t  !  ", @libs) ."\n\n" if(@libs);
    print $self->dist_name." PACKAGE:\t\t";
    print $stat;
    print "\n";

    print "\nImportant Commands:\n".
        "\t./Build installdeps\t\#installs missing PERL dependencies\n".
        "\t./Build installexes\t\#installs missing external programs\n".
        "\t./Build install\t\t\#installs this package, (".$self->dist_name.")\n".
        "\t./Build status\t\t\#Shows this status menu\n\n".
        "Other Commands:\n".
        "\t./Build htslib\t\#installs htslib\n";
}

#test if there is another version of the module overriding the CPAN install
sub module_overide {
    my $self = shift;
    my $desired = shift;

    my $mod = $desired; #holds expected .pm file name
    $mod =~ s/\:\:/\//g;
    $mod .= '.pm';
    
    my $test=  qq(\@INC = qw($Config{installsitelib}
			     $Config{installsitearch}
			     $Config{installvendorlib}
			     $Config{installvendorarch}
			     $Config{installprivlib}
			     $Config{installarchlib});
		  require $desired;
		  print \$INC{\"$mod\"};);
    
    my $ok = `$^X -e 'eval q{$test} or exit 1'`;
    my $loc = $self->module_loc($desired) if($ok);
    
    return ($loc && $loc ne $ok) ? $loc : undef;
}

#gets the location of a module
sub module_loc {
    my $self = shift;
    my $desired = shift;

    return if(! $desired);

    eval "require $desired"; #loads module into \%INC    

    $desired =~ s/\:\:/\//g;
    $desired .= ".pm";

    return $INC{$desired};
}

sub sync_dirs {
    my $self = shift;
    my @dirs = @_;
    my $sbase = $self->base_dir; #source base
    my $ibase = $self->install_base; #install base

    #find all files in src locations
    my %sfound;
    my @sfiles;
    my @sdirs = map {"$sbase/$_"} @dirs;
    foreach my $dir (@sdirs){
	foreach my $o (File::Glob::bsd_glob("$dir/*")){
	    next if($o =~ /(?:~|\/\.[^\/]+|\.PL|\/\.svn)\Z/);
	    if(-d $o){
		push(@sdirs, $o);
	    }
	    else{
		push(@sfiles, $o);
		$sfound{substr($o, length($sbase)+1)}++;
	    }
	}
    }

    #find files in install locations
    my @extra;
    my @idirs  = ("$ibase/bin", "$ibase/lib");
    foreach my $dir (@idirs){
	foreach my $o (File::Glob::bsd_glob("$dir/*")){
	    next if($o =~ /(?:~|\.PL|\/\.svn)\Z/);
	    if(-d $o){
		push(@idirs, $o);
		next;
	    }
	    elsif($sfound{my $frel = substr($o, length($ibase)+1)}){ #matching file
		my $sfile = "$sbase/$frel";
		my $ifile  = "$ibase/$frel";

		#write permission must have been set or no editing occured
		next unless(sprintf("%04o", (stat($ifile))[2] & 07777) =~ /[2367]/);

		#get file contents to compare for modifications (strip of shebang first)
		my $sheader_found;
		my $iheader_found;
		my $sdata = load_w_o_header($sfile, \$sheader_found);
		my $idata = load_w_o_header($ifile, \$iheader_found);
		next if($sdata eq $idata);

		#get modification times to see whiuch is newer
		my $smod = (stat($sfile))[9];
		my $imod = (stat($ifile))[9];
		next if($imod <= $smod);

		#copy file
		my $err;
		print "copying $ifile  -->  $sfile\n";
		File::Copy::move($sfile, "$sfile.bk~"); #backup incase of failure
		if(open(IN, "> $sfile")){
		    print IN "#!/usr/bin/perl\n\n" if($sheader_found || $iheader_found);
		    print IN $idata;
		    close(IN);
		}
		else{
		    $err = $!;
		}
		    
		#restore file on failure
		if(! -f $sfile && -f "$sfile.bk~"){
		    File::Copy::move("$sfile.bk~", $sfile);
		    die "ERROR: Could not copy $ifile to $sfile\n$err";
		}
		else{
		    unlink("$sfile.bk~");
		}
	    }
	    else{ #extra files not included in source directory
		next if($frel =~ /\/ConfigData\.pm\Z/);
		next if($frel =~ /\Alib\/auto\//);
		push(@extra, $frel);
	    }
	}
    }

    #identify which extra files to copy
    foreach my $frel (@extra){
	my $sfile = "$sbase/$frel";
	my $ifile  = "$ibase/$frel";

	my $go = $self->y_n("New file found that is not in the 'src' directory. Do you want to add ../$frel?", 'N');
	next unless ($go);

	print "copying $ifile  -->  $sfile\n";
	my $iheader_found;
	my $idata = load_w_o_header($ifile, \$iheader_found);
	File::Path::make_path(File::Basename::dirname($sfile)); #make containing directory
	open(IN, "> $sfile") or die "ERROR: Could not copy $ifile to $sfile\n$!";
	print IN "#!/usr/bin/perl\n\n" if($iheader_found);
	print IN $idata;
	close(IN);
    }
}

sub svn_w_args {
    my $self = shift;
    my @args = @_;

    my $svn = App::Which::which("svn");
    die "ERROR: Cannot find the executable 'svn' (subversion respository tool)\n" if(!$svn);
    die "ERROR: Failure to supply a subcommand to subversion\n" if(!@args);

    chdir($self->base_dir."/../"); #change to subversion root
    my $pid = open3(undef, \*CHILD_OUT, undef, $svn, @args);
    my $out = join('', <CHILD_OUT>);
    waitpid( $pid, 0 );
    chdir($self->base_dir);

    return $out;
}

sub check_update_version {
    my $self = shift;

    #get current subversion version
    my ($svn) = $self->svn_w_args('info') =~ /Revision\:\s*(\d+)/;
    die "ERROR: Could not query subversion repository\n" if(!$svn);

    #get old version information for last stable release
    my $file = $self->dist_version_from;
    open(IN, "< $file") or die "ERROR: Could not open version file\n";
    my $data = join('', <IN>);
    close(IN);
    my ($old_svn) = $data =~ /\$SVN[\s\t]*=[\s\t]*'?(\d+)/;
    $old_svn ||= 0;
    my ($old_version) = $data =~ /\$VERSION[\s\t]*=[\s\t]*'?([\d\.]+)/;
    $old_version ||= 0;

    #check if update is really needed
    if($old_svn == $svn){
        print $self->dist_name." is already up to date as stable release $old_version\n";
        return $old_version;
    }
    else{
        #set new version
	my $version;
        if(my @v = $old_version =~ /([\d]+)/g){
	    my ($m, $n, $s) = @v;
	    if(defined($s)){ #iterate subversion
		my $l = length($s);
		$s = ($l > 1 && $s =~ /^0/) ? sprintf('%0'.$l.'s', $s+1) : $s+1;
		$version = "$m.$n.$s";
	    }
	    elsif(defined($n)){ #iterate minor version
		my $l = length($n);
		$n = ($l > 1 && $n =~ /^0/) ? sprintf('%0'.$l.'s', $n+1) : $n+1;
		$version = "$m.$n";
	    }
	    elsif(defined($m)){ #iterate major version
		my $l = length($m);
		$m = ($l > 1 && $m =~ /^0/) ? sprintf('%0'.$l.'s', $m+1) : $m+1;
		$version = "$m";
	    }
	}
	$svn++; #iterate to next version

	#see if subversion repository is using trunk
	my ($trunk) = $self->svn_w_args('info') =~ /URL\: ([^\s\t\n]+\/trunk)[\s\t\n]/;
	$svn++ if($trunk); #because I will commit a second time

	#changing script version here
	$data =~ s/\$VERSION[\s\t]*=[\s\t]*'?[\d\.]+'?/\$VERSION = \'$version\'/;
	$data =~ s/\$SVN[\s\t]*=[\s\t]*'?\d+'?/\$SVN = \'$svn\'/;
	open(OUT, "> $file");
	print OUT $data;
	close(OUT);
	
	$self->svn_w_args('commit', '-m', "stable release version $version");
	$self->svn_w_args('update');

	#make tag of release if using trunk
	if($trunk){
	    (my $tag = $trunk) =~ s/trunk$/tags\/Version_$version\_r$svn/;
	    my $copy_message = "Adding tags/Version_$version\_r$svn";
	    $self->svn_w_args('copy', $trunk, $tag, '-m', $copy_message);
	}

        print $self->dist_name." has been updated to stable release $version\n";

        return $version;
    }
}

#trims modified shebang lines off of file
sub load_w_o_header {
    my $file = shift;
    my $flag = shift; #optional flag to signal that header was found

    my $header;
    my $data = '';
    open(IN, "< $file");
    while(my $line = <IN>){
        #strip of perl shebang portiion
        while($line =~ /^\#!.*perl/ ||
              $line =~ /exec\s+.*perl/ ||
              $line =~ /if\s*0\;/ ||
              $line =~ /^[\s\t\n]+$/
	    ){
	    $$flag = 1 if(ref($flag) eq 'SCALAR' || ref($flag) eq 'REF');
            $line = <IN>;
        }
        $data = join('', $line, <IN>);
    }
    close(IN);

    return $data;
}

sub safe_prompt {
    require Term::ReadKey;

    my $self = shift;
    my $m = shift;
    my $d = shift || '';

    print "$m [".('*'x(length($d)))." ]";

    my $key = 0;
    my $r = "";
    #Start reading the keys
    Term::ReadKey::ReadMode(4); #Disable the control keys (raw mode)

    while(ord($key = Term::ReadKey::ReadKey(0)) != 10) { #This will continue until the Enter key is pressed (decimal value of 10)
	if(ord($key) == 127 || ord($key) == 8) { #DEL/Backspace was pressed
	    if(length($r) > 0){
		#1. Remove the last char from the password
		chop($r);
		#2 move the cursor back by one, print a blank character, move the cursor back by one
		print "\b \b";
	    }
	} elsif(ord($key) < 32) {
	    # Do nothing with these control characters
	} else {
	    $r = $r.$key;
	    print "*";
	}
    }
    print "\n"; #because the user pressed enter
    Term::ReadKey::ReadMode(0); #Reset the terminal once we are done

    $r = $d if(length($r) == 0);
    return $r; #Return the response
}


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------


package App::Which;

use strict;
use warnings;
use File::Spec ();

# ABSTRACT: Perl implementation of the which utility as an API
our $VERSION = '1.27'; # VERSION

use constant IS_VMS => ($^O eq 'VMS');
use constant IS_MAC => ($^O eq 'MacOS');
use constant IS_WIN => ($^O eq 'MSWin32' or $^O eq 'dos' or $^O eq 'os2');
use constant IS_DOS => IS_WIN();
use constant IS_CYG => ($^O eq 'cygwin' || $^O eq 'msys');

our $IMPLICIT_CURRENT_DIR = IS_WIN || IS_VMS || IS_MAC;

# For Win32 systems, stores the extensions used for
# executable files
# For others, the empty string is used
# because 'perl' . '' eq 'perl' => easier
my @PATHEXT = ('');
if ( IS_WIN ) {
  # WinNT. PATHEXT might be set on Cygwin, but not used.
    if ( $ENV{PATHEXT} ) {
	push @PATHEXT, split /;/, $ENV{PATHEXT};
    } else {
    # Win9X or other: doesn't have PATHEXT, so needs hardcoded.
	push @PATHEXT, qw{.com .exe .bat};
    }
} elsif ( IS_VMS ) {
    push @PATHEXT, qw{.exe .com};
} elsif ( IS_CYG ) {
  # See this for more info
  # http://cygwin.com/cygwin-ug-net/using-specialnames.html#pathnames-exe
    push @PATHEXT, qw{.exe .com};
}


sub which {
    my ($exec) = @_;

    return undef unless defined $exec;
    return undef if $exec eq '';

    my $all = wantarray;  ## no critic (Freenode::Wantarray)
    my @results = ();

  # check for aliases first
    if ( IS_VMS ) {
	my $symbol = `SHOW SYMBOL $exec`;
	chomp($symbol);
	unless ( $? ) {
	    return $symbol unless $all;
	    push @results, $symbol;
	}
    }
    if ( IS_MAC ) {
	my @aliases = split /\,/, $ENV{Aliases};
	foreach my $alias ( @aliases ) {
      # This has not been tested!!
      # PPT which says MPW-Perl cannot resolve `Alias $alias`,
      # let's just hope it's fixed
	    if ( lc($alias) eq lc($exec) ) {
		chomp(my $file = `Alias $alias`);
		last unless $file;  # if it failed, just go on the normal way
		return $file unless $all;
		push @results, $file;
        # we can stop this loop as if it finds more aliases matching,
        # it'll just be the same result anyway
		last;
	    }
	}
    }

  return $exec  ## no critic (ValuesAndExpressions::ProhibitMixedBooleanOperators)
      if !IS_VMS and !IS_MAC and !IS_WIN and $exec =~ /\// and -f $exec and -x $exec;

    my @path;
    if($^O eq 'MSWin32') {
    # File::Spec (at least recent versions)
    # add the implicit . for you on MSWin32,
    # but we may or may not want to include
    # that.
	@path = split /;/, $ENV{PATH};
    s/"//g for @path;
    @path = grep length, @path;
  } else {
    @path = File::Spec->path;
  }
  if ( $IMPLICIT_CURRENT_DIR ) {
    unshift @path, File::Spec->curdir;
  }

  foreach my $base ( map { File::Spec->catfile($_, $exec) } @path ) {
    for my $ext ( @PATHEXT ) {
      my $file = $base.$ext;

      # We don't want dirs (as they are -x)
      next if -d $file;

      if (
        # Executable, normal case
        -x _
        or (
          # MacOS doesn't mark as executable so we check -e
          IS_MAC  ## no critic (ValuesAndExpressions::ProhibitMixedBooleanOperators)
          ||
          (
            ( IS_WIN or IS_CYG )
            and
	   grep {   ## no critic (BuiltinFunctions::ProhibitBooleanGrep)
              $file =~ /$_\z/i
	   } @PATHEXT[1..$#PATHEXT]
          )
          # DOSish systems don't pass -x on
          # non-exe/bat/com files. so we check -e.
          # However, we don't want to pass -e on files
          # that aren't in PATHEXT, like README.
          and -e _
        )
	  ) {
	  return $file unless $all;
	  push @results, $file;
      }
    }
  }

    if ( $all ) {
	return @results;
    } else {
	return undef;
    }
}


sub where {
  # force wantarray
    my @res = which($_[0]);
    return @res;
}

1;
