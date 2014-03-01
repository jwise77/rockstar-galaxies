#!/usr/bin/perl -w

my $snapshot = $ARGV[0];
die "Usage: $0 snapnum\n" unless (defined($snapshot) and $snapshot >= 0);

our $rockstar_path = "/path/to/rockstar";
our $configfile = "test.cfg";
our $server_output = "server_$snapshot.dat";
our $client_output = "clients_$snapshot.dat";
our $output = "-o";
our $queue = "-q kipac-ibq";
our $numproc = "-n";
our $total_numproc = 256;
our $machine_procs = 8;
our $batch_command = "bsub";
our $batch_jobs = "bjobs";

unlink("auto-rockstar.cfg");
system("$batch_command $queue $output $server_output ".
       "$rockstar_path -c $configfile -s $snapshot");
sleep 1 while (!(-e "auto-rockstar.cfg")); #wait for server to start

for (1..($total_numproc / $machine_procs)) { #Launch analysis processes
  system("$batch_command $queue $numproc $machine_procs $output $client_output".
	 "$rockstar_path -c auto-rockstar.cfg");
}

sleep 1 while (grep { /rockstar/ } `bjobs`); #Wait for jobs to finish or time out

$success = 0;
open SERVER, "<", $server_output;
my @output = <SERVER>;
close SERVER;
$success = 1 if (grep {/Success/} @output);
$success = 0 if (grep {/Error/} @output);

open CLIENTS, "<", $client_output;
@output = <CLIENTS>;
close CLIENTS;
$success = 0 if (grep {/Error/} @output);

if ($success) {
    #Do operations for successful finish
    system("./rockstar_success.sh $snapshot");
} else {
    system("./rockstar_failure.sh $snapshot");
}

