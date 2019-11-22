#!/usr/bin/env perl
#
# Author: petr.danecek@sanger
#

use strict;
use warnings;
use Carp;
use List::Util 'shuffle';
use File::Temp qw/ tempdir /;

my $opts = parse_params();
for (my $i=0; $i<$$opts{ntests}; $i++)
{
    run_test($opts);
}

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print 
        "Usage: test-splicing.pl [OPTIONS]\n",
        "Options:\n",
        "   -e, --perm-test-exec <file>     \n",
        "   -n, --ntests <int>              Number of tests to run\n",
        "   -s, --random-seed <int>         \n",
        "   -t, --temp-dir <dir>            \n",
        "   -h, -?, --help                  This help message.\n",
        "\n";
    exit -1;
}
sub parse_params
{
    my $opts = { ntests=>100, exec=>'perm-test' };
    while (defined(my $arg=shift(@ARGV)))
    {
        if ( $arg eq '-e' || $arg eq '--perm-test-exec' ) { $$opts{exec}=shift(@ARGV); next }
        if ( $arg eq '-n' || $arg eq '--ntests' ) { $$opts{ntests}=shift(@ARGV); next }
        if ( $arg eq '-s' || $arg eq '--random-seed' ) { $$opts{seed}=shift(@ARGV); next }
        if ( $arg eq '-t' || $arg eq '--temp-dir' ) { $$opts{tmp_dir}=shift(@ARGV); next }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    if ( !exists($$opts{tmp_dir}) )
    {
        $$opts{tmp_dir} = tempdir( CLEANUP => 1 );
        print STDERR "Using tempdir $$opts{tmp_dir}\n";
    }
    else
    {
        `mkdir -p $$opts{tmp_dir}`;
        if ( $? ) { error("mkdir -p $$opts{tmp_dir}: $!"); }
    }
    if ( !exists($$opts{seed}) )
    {
        $$opts{seed} = time;
        print STDERR "Using random seed $$opts{seed}\n";
    }
    srand($$opts{seed});
    return $opts;
}
sub run_test
{
    my ($opts) = @_;
    my $nregs = int(rand(1000)) + 1;
    my @exp  = ();
    my @regs = ();
    my @prev = ();
    my $end  = 0;
    my @types = ('SKIP','BG','TGT');
    for (my $ireg=0; $ireg<$nregs; $ireg++)
    {
        my $type = int(rand(3));    # 0:skip, 1:bg, 2:tgt
        if ( @prev && $prev[-1] eq $type ) { $type = ++$type % 3; }
        push @prev, $type;
        if ( @prev>3 ) { shift @prev; }
        my $len = int(rand(10000)) + 1;
        if ( $type eq 0 ) { $end += $len; next; }
        my $beg = $end + 1;
        $end += $len;
        push @exp,{ type=>$types[$type], beg=>$beg, end=>$end, len=>$end-$beg+1 };
        push @regs,{ type=>$type, beg=>$beg, end=>$end };
        my $nembed = int(rand(100)) + 1;
        for (my $iembed=0; $iembed<$nembed; $iembed++)
        {
            my $ebeg = $beg + int(rand($len));
            my $eend = $ebeg + int(rand($end-$ebeg+1));
            push @regs,{ type=>$type, beg=>$ebeg, end=>$eend };
        }
        if ( join('',@prev) eq '121' )
        {
            # let the target region splice out the background region
            my $beg0 = $exp[-3]{beg};
            my $len0 = $exp[-3]{len} + $exp[-2]{len} + $exp[-1]{len};
            $nembed = int(rand(100)) + 1;
            for (my $iembed=0; $iembed<$nembed; $iembed++)
            {
                my $ebeg = $beg0 + int(rand($len0));
                my $eend = $ebeg + int(rand($end-$ebeg+1));
                push @regs,{ type=>$type, beg=>$ebeg, end=>$eend };
            }
        }
    }
    @regs = shuffle(@regs);

    open(my $fh_exp,'>',"$$opts{tmp_dir}/exp.txt") or error("$$opts{tmp_dir}/exp.txt: $!");
    for my $exp (@exp)
    {
        print $fh_exp "$$exp{type}\t$$exp{beg}\t$$exp{end}\n";
    }
    close($fh_exp) or error("close failed: $$opts{tmp_dir}/exp.txt");

    open(my $fh_fai,'>',"$$opts{tmp_dir}/ref.fai") or error("$$opts{tmp_dir}/ref.fai: $!");
    print $fh_fai "1\t".($end+1)."\n";
    print $fh_fai "22\t".($end+1)."\n";
    print $fh_fai "333\t".($end+1)."\n";
    close($fh_fai) or error("close failed: $$opts{tmp_dir}/ref.fai");

    open(my $fh_tgt,'>',"$$opts{tmp_dir}/tgt.txt") or error("$$opts{tmp_dir}/tgt.txt: $!");
    open(my $fh_bg,'>',"$$opts{tmp_dir}/bg.txt") or error("$$opts{tmp_dir}/bg.txt: $!");
    open(my $fh_call,'>',"$$opts{tmp_dir}/call.txt") or error("$$opts{tmp_dir}/call.txt: $!");
    for my $reg (@regs)
    {
        if ( defined $fh_call && $$reg{type} eq 2 )
        {
            print $fh_call "1\t$$reg{beg}\t$$reg{end}\n";
            close($fh_call) or error("close failed: $$opts{tmp_dir}/call.txt");
            $fh_call = undef;
        }
        my $fh = $$reg{type} eq 1 ? $fh_bg : $fh_tgt;
        print $fh "1\t$$reg{beg}\t$$reg{end}\n";
        print $fh "22\t$$reg{beg}\t$$reg{end}\n";
        print $fh "333\t$$reg{beg}\t$$reg{end}\n";
    }
    close($fh_tgt) or error("close failed: $$opts{tmp_dir}/tgt.txt");
    close($fh_bg) or error("close failed: $$opts{tmp_dir}/bg.txt");

    my $cmd = qq[$$opts{exec} -t $$opts{tmp_dir}/tgt.txt -b $$opts{tmp_dir}/bg.txt -c $$opts{tmp_dir}/call.txt -f $$opts{tmp_dir}/ref.fai -o $$opts{tmp_dir}/out.txt -d];
    print STDERR "$cmd\n";
    `$cmd`;
    if ( $? ) { error("The program exited with an error: $!\n\t$cmd\n"); }

    my @reg1 = ();
    my @reg2 = ();
    my @reg3 = ();
    open(my $fh,'<',"$$opts{tmp_dir}/out.txt") or error("$$opts{tmp_dir}/out.txt: $!");
    while (my $line=<$fh>)
    {
        if ( !($line=~/^BG/) and !($line=~/^TGT/) ) { next; }
        my @cols = split(/\t/,$line);
        chomp($cols[-1]);
        if ( $cols[1] eq '1' ) { push @reg1,"$cols[0]\t$cols[2]\t$cols[3]"; }
        elsif ( $cols[1] eq '22' ) { push @reg2,"$cols[0]\t$cols[2]\t$cols[3]"; }
        elsif ( $cols[1] eq '333' ) { push @reg3,"$cols[0]\t$cols[2]\t$cols[3]"; }
    }
    close($fh) or error("close failed: $$opts{tmp_dir}/out.txt");
    if ( @exp != @reg1 ) { error(sprintf("Different number of regions: %d vs %d\n", scalar @exp,scalar @reg1)); }
    if ( @exp != @reg2 ) { error(sprintf("Different number of regions: %d vs %d\n", scalar @exp,scalar @reg2)); }
    if ( @exp != @reg3 ) { error(sprintf("Different number of regions: %d vs %d\n", scalar @exp,scalar @reg3)); }
    for (my $i=0; $i<@exp; $i++)
    {
        my $exp = "$exp[$i]{type}\t$exp[$i]{beg}\t$exp[$i]{end}";
        if ( $exp ne $reg1[$i] ) { error(sprintf("Different regions: $exp[$i] vs $reg1[$i]\n")); }
        if ( $exp ne $reg2[$i] ) { error(sprintf("Different regions: $exp[$i] vs $reg2[$i]\n")); }
        if ( $exp ne $reg3[$i] ) { error(sprintf("Different regions: $exp[$i] vs $reg3[$i]\n")); }
    }
}


