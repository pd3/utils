#!/usr/bin/env perl
#
#   Copyright (C) 2018 Genome Research Ltd.
#
#   Author: Petr Danecek <pd3@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

use strict;
use warnings;
use Carp;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;
use Cwd qw/ abs_path /;

my $opts = parse_params();
test_annot_regs($opts,src=>'src.1.txt',dst=>'dst.1.txt',out=>'out.1.1.txt',args=>'-t smpl:overlap --allow-dups');
test_annot_regs($opts,src=>'src.1.txt',dst=>'dst.1.txt',out=>'out.1.2.txt',args=>'-t smpl:overlap');
test_annot_regs($opts,src=>'src.1.txt',dst=>'dst.1.txt',out=>'out.1.2.txt',args=>'-t smpl:overlap -c chr,beg,end');
test_annot_regs($opts,src=>'src.1.txt',dst=>'dst.1.txt',out=>'out.1.3.txt',args=>'-t smpl,value:overlap,value');
test_annot_regs($opts,src=>'src.1.txt',dst=>'dst.1.txt',out=>'out.1.4.txt',args=>'-t smpl:overlap -o 0.5');
test_annot_regs($opts,src=>'src.1.txt',dst=>'dst.1.txt',out=>'out.1.5.txt',args=>'-t smpl:overlap -ro 0.5');
test_annot_regs($opts,src=>'src.2.txt',dst=>'dst.2.txt',out=>'out.2.1.txt',args=>'-c 1,2,3:1,2,3 -t 4:5 --allow-dups');
test_annot_regs($opts,src=>'src.2.txt',dst=>'dst.2.txt',out=>'out.2.2.txt',args=>'-c 1,2,3:1,2,3 -t 4:5');
test_annot_regs($opts,src=>'src.2.txt',dst=>'dst.2.txt',out=>'out.2.3.txt',args=>'-c 1,2,3:1,2,3 -t 4,value:5,value');
test_annot_regs($opts,src=>'src.2.txt',dst=>'dst.2.txt',out=>'out.2.4.txt',args=>'-c 1,2,3:1,2,3 -t value,4:value,5');
test_annot_regs($opts,src=>'src.2.txt',dst=>'dst.2.txt',out=>'out.2.5.txt',args=>'-c 1,2,3:1,2,3 -t value,4:value,5 -a nbp,frac');
test_annot_regs($opts,src=>'src.3.txt',dst=>'dst.3.txt',out=>'out.3.1.txt',args=>'-t smpl:overlap -a nbp,frac');
test_annot_regs($opts,src=>'src.4.txt',dst=>'dst.4.txt',out=>'out.4.1.txt',args=>'-c 2,3,4:2,3,4 -m 1:1 -t 1:1 -a nbp,frac');
test_annot_regs($opts,src=>'src.5.txt',dst=>'dst.5.txt',out=>'out.5.1.txt',args=>'-c 2,3,4:2,3,4 -a nbp,frac');
test_annot_regs($opts,src=>'src.6.txt',dst=>'dst.6.txt',out=>'out.6.1.txt',args=>'-c 1,2,2:1,2,2 -a nbp');

print "\nNumber of tests:\n";
printf "    total   .. %d\n", $$opts{nok}+$$opts{nfailed};
printf "    passed  .. %d\n", $$opts{nok};
printf "    failed  .. %d\n", $$opts{nfailed};
print "\n";

exit ($$opts{nfailed} != 0);

#--------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print
        "About: htslib consistency test script\n",
        "Usage: test.pl [OPTIONS]\n",
        "Options:\n",
        "   -p, --plugins                   Test also plugins, requires libhts.so.\n",
        "   -r, --redo-outputs              Recreate expected output files.\n",
        "   -t, --temp-dir <path>           When given, temporary files will not be removed.\n",
        "   -h, -?, --help                  This help message.\n",
        "\n";
    exit -1;
}
sub parse_params
{
    my $opts = { bgzip=>"bgzip", keep_files=>0, nok=>0, nfailed=>0, tabix=>"tabix", plugins=>0 };
    my $help;
    Getopt::Long::Configure('bundling');
    my $ret = GetOptions (
            'e|exec=s' => sub { my ($tool, $path) = split /=/, $_[1]; $$opts{$tool} = $path if $path },
            't|temp-dir:s' => \$$opts{keep_files},
            'p|plugins' => \$$opts{test_plugins},
            'r|redo-outputs' => \$$opts{redo_outputs},
            'h|?|help' => \$help
            );
    if ( !$ret or $help ) { error(); }
    $$opts{tmp} = $$opts{keep_files} ? $$opts{keep_files} : tempdir(CLEANUP=>1);
    if ( $$opts{keep_files} ) { cmd("mkdir -p $$opts{keep_files}"); }
    $$opts{path} = $FindBin::RealBin;
    $$opts{bin}  = $FindBin::RealBin;
    $$opts{bin}  =~ s{/test/?$}{};
    return $opts;
}
sub _cmd
{
    my ($cmd) = @_;
    my $kid_io;
    my @out;
    my $pid = open($kid_io, "-|");
    if ( !defined $pid ) { error("Cannot fork: $!"); }
    if ($pid)
    {
        # parent
        @out = <$kid_io>;
        close($kid_io);
    }
    else
    {
        # child
        exec('/bin/bash', '-o','pipefail','-c', $cmd) or error("Cannot execute the command [/bin/sh -o pipefail -c $cmd]: $!");
    }
    return ($? >> 8, join('',@out));
}
sub cmd
{
    my ($cmd) = @_;
    my ($ret,$out) = _cmd($cmd);
    if ( $ret ) { error("The command failed: $cmd\n", $out); }
    return $out;
}
sub test_cmd
{
    my ($opts,%args) = @_;
    if ( !exists($args{out}) )
    {
        if ( !exists($args{in}) ) { error("FIXME: expected out or in key\n"); }
        $args{out} = "$args{in}.out";
    }
    my ($package, $filename, $line, $test)=caller(1);
    $test =~ s/^.+:://;

    print "$test:\n";
    print "\t$args{cmd}\n";

    my ($ret,$out) = _cmd("$args{cmd}");
    if ( $ret ) { failed($opts,$test,"Non-zero status $ret"); return; }
    if ( $$opts{redo_outputs} && -e "$$opts{path}/$args{out}" )
    {
        rename("$$opts{path}/$args{out}","$$opts{path}/$args{out}.old");
        open(my $fh,'>',"$$opts{path}/$args{out}") or error("$$opts{path}/$args{out}: $!");
        print $fh $out;
        close($fh);
        my ($ret,$out) = _cmd("diff -q $$opts{path}/$args{out} $$opts{path}/$args{out}.old");
        if ( !$ret && $out eq '' ) { unlink("$$opts{path}/$args{out}.old"); }
        else
        {
            print "\tthe expected output changed, saving:\n";
            print "\t  old .. $$opts{path}/$args{out}.old\n";
            print "\t  new .. $$opts{path}/$args{out}\n";
        }
    }
    my $exp = '';
    if ( exists($args{exp}) ) { $exp = $args{exp}; }
    elsif ( open(my $fh,'<',"$$opts{path}/$args{out}") )
    {
        my @exp = <$fh>;
        $exp = join('',@exp);
        close($fh);
    }
    else
    {
        open(my $fh,'>',"$$opts{path}/$args{out}.new") or error("$$opts{path}/$args{out}.new: $!");
        print $fh $out;
        close($fh);
        if ( !$$opts{redo_outputs} ) { failed($opts,$test,"$$opts{path}/$args{out}.new"); return; }
    }

    if ( $exp ne $out )
    {
        open(my $fh,'>',"$$opts{path}/$args{out}.new") or error("$$opts{path}/$args{out}.new");
        print $fh $out;
        close($fh);
        if ( !-e "$$opts{path}/$args{out}" && !exists($args{exp}) )
        {
            rename("$$opts{path}/$args{out}.new","$$opts{path}/$args{out}") or error("rename $$opts{path}/$args{out}.new $$opts{path}/$args{out}: $!");
            print "\tthe file with expected output does not exist, creating new one:\n";
            print "\t\t$$opts{path}/$args{out}\n";
        }
        else
        {
            if ( exists($args{exp}) )
            {
                open(my $fh,'>',"$$opts{path}/$args{out}") or error("$$opts{path}/$args{out}");
                print $fh $exp;
                close($fh);
            }
            failed($opts,$test,"The outputs differ:\n\t\t$$opts{path}/$args{out}\n\t\t$$opts{path}/$args{out}.new");
        }
        return;
    }
    passed($opts,$test);
}
sub failed
{
    my ($opts,$test,$reason) = @_;
    $$opts{nfailed}++;
    if ( defined $reason ) { print "\n\t$reason"; }
    print "\n.. failed ...\n\n";
}
sub passed
{
    my ($opts,$test) = @_;
    $$opts{nok}++;
    print ".. ok\n\n";
}
sub is_file_newer
{
    my ($afile,$bfile) = @_;
    my (@astat) = stat($afile) or return 0;
    my (@bstat) = stat($bfile) or return 0;
    if ( $astat[9]>$bstat[9] ) { return 1 }
    return 0;
}



# The tests --------------------------

sub test_annot_regs
{
    my ($opts,%args) = @_;
    my $args  = exists($args{args}) ? $args{args} : '';
    test_cmd($opts,%args,cmd=>"$$opts{bin}/annot-regs $args -s $$opts{path}/$args{src} -d $$opts{path}/$args{dst}");
}

