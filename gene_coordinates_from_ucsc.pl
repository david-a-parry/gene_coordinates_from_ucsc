#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use DBI;
use LWP::Simple;
use XML::Simple; 
use Pod::Usage;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/lib";
use SortGenomicCoordinates;
my @genes = ();
my %opts = (gene => \@genes);
GetOptions(\%opts,
            "gene=s{,}",
            "list=s",
            "cds",
            "merge",
            "bed", #output bed rather than regions
            "ensembl", #search ensembl transcripts rather than refGene
            "ucsc", #search ucsc known gene transcripts rather than refGene
            "db=s", #genome database
            "no_chr",
            "help",
            "manual",
        ) or pod2usage(-exitval => 2, -message => "Syntax error");
pod2usage (-verbose => 2) if $opts{manual};
pod2usage (-verbose => 1) if $opts{help};
pod2usage(-exitval => 2, -message => "--gene or --list argument required - for help use --help") if not @{$opts{gene}} and not $opts{list};
pod2usage(-exitval => 2, -message => "--ensembl or --ucsc arguments are mutually exclusisve, please choose only one") if $opts{ensembl} and $opts{ucsc};
if ($opts{list}){
    open (my $LIST, $opts{list}) or die "Can't open $opts{list} for reading: $!\n";
    while (<$LIST>){
        chomp;
        my @split = split(/\s/);
        push @genes, $split[0];
    }
}

my $db = 'hg19';#default genome
if ($opts{db}){
    $db = $opts{db};
}
my $transcript_db = 'refGene';#default to RefSeq genes
my $name_search = 'name2';#search by symbol using name2
if ($opts{ensembl}){
    $transcript_db = 'ensGene';
    $name_search = 'name';#use name - have to get name using symbol from ensemblToGeneName db
}elsif($opts{ucsc}){
    $transcript_db = 'knownGene';
    $name_search = 'name';#use name - have to get name using symbol from ensemblToGeneName db
}
#default to transcript start and end sites
my $start_fetch = "txStart";
my $end_fetch = "txEnd";
if ($opts{cds}){
    $start_fetch = "cdsStart";
    $end_fetch = "cdsEnd";
}
#FOR ENSEMBL QUERIES
#we need to get transcript names in a separate step
my %names = ();
my $dbh = DBI->connect_cached("dbi:mysql:$db:genome-mysql.cse.ucsc.edu", "genome", '')
  or die "Error connecting to 'dbi:mysql:$db:genome-mysql.cse.ucsc.edu' genome.\n";
foreach my $gene (@genes){
    if ($opts{ensembl}){
        my $command = "SELECT name FROM $db.ensemblToGeneName  WHERE value=?";
        my $sth = $dbh->prepare($command);
        $sth->execute($gene);
        while (my @row = $sth->fetchrow_array ) { 
            push @{$names{$gene}}, $row[0];
        }
    }elsif($opts{ucsc}){
        my $command = "SELECT kgID FROM $db.kgAlias WHERE alias=?";
        my $sth = $dbh->prepare($command);
        $sth->execute($gene);
        while (my @row = $sth->fetchrow_array ) { 
            push @{$names{$gene}}, $row[0];
        }
    }else{
    #for RefSeq we can just use the symbol
        push @{$names{$gene}}, $gene;
    }
}
my @all_regions = ();#collect regions here and merge if necessarry

foreach my $gene (@genes){
    if (not exists $names{$gene} ){
        print STDERR "Warning - no transcript found for gene $gene.\n";
        next;
    }
    if (not @{$names{$gene}}){
        print STDERR "Warning - no transcript found for gene $gene.\n";
        next;
    }
    foreach my $name (@{$names{$gene}}){
        my $found = 0;
        my $command = "SELECT chrom, $start_fetch, $end_fetch, strand, name FROM $db.$transcript_db WHERE $name_search=?";  
        my $sth = $dbh->prepare($command);
        $sth->execute($name);
        while (my  @row = $sth->fetchrow_array ) {
            my $region = join("\t", (@row[0..4])) ."[$gene]";
            if ($opts{no_chr}){
                $region =~ s/^chr//;
            }
            push @all_regions, $region;
            $found++;
        }
        if (not $found){
            print STDERR "Warning - no transcript found for $name";
            if ($gene ne $name){
                print "/$gene";
            }
            print ".\n";
        }
    }
}

if (not @all_regions){
    die "No regions identified from user input.\n";
}

#SORT AND MERGE OUR EXONS+FLANKS ARRAY
if ($opts{merge}){
    my $sort_obj = SortGenomicCoordinates -> new(array => \@all_regions, type => "bed");
    $sort_obj -> order();
    my $merged = $sort_obj -> merge();
    $sort_obj -> DESTROY();
    #print Dumper $merged;
    foreach my $reg(@$merged){
        my @names = ();
        my %strand = ();
        foreach my $inf (@{$reg->{info}}){
            my @i = split("\t", $inf);
            push @names,  $i[4];
            $strand{$i[3]}++;
        }
        my $str;
        if (keys %strand > 1){
            print STDERR "WARNING - MORE THAN ONE VALUE FOR STRAND IDENTIFIED\n";
            my $n = 0;
            foreach my $k (keys %strand){
                $strand{$k} > $n ? $str = $k : $str = $str;
                $n = $strand{$k};
            }
        }else{
            $str = (keys %strand)[0];
        }
        
        my $out = join ("\t", ( "chr$reg->{chrom}", "$reg->{start}", "$reg->{end}", "$str")) ."\t" 
            .join("/", @names) ."\n";
        if (not $opts{bed}){
            $out =~ s/\t/:/;
            $out =~ s/\t/-/;
        }
        if ($opts{no_chr}){
            $out =~ s/^chr//;
        }
        print $out;
    }
}else{
    if ($opts{bed}){
        print join("\n", @all_regions) ."\n";
    }else{
        foreach my $reg (@all_regions){
            $reg =~ s/\t/:/;
            $reg =~ s/\t/-/;
            print "$reg\n";
        }
    }
}



=head1 NAME

gene_coordinates_from_ucsc.pl - get transcript coordinates for RefSeq, Ensembl or UCSC genes

=head1 SYNOPSIS

        gene_coordinates_from_ucsc.pl -g [gene symbol] [options]
        gene_coordinates_from_ucsc.pl -l [gene_list.txt] [options]
        gene_coordinates_from_ucsc.pl -h (display help message)
        gene_coordinates_from_ucsc.pl --manual (display manual page)


=cut

=head1 ARGUMENTS

=over 8

=item B<-g    --gene>

One or more gene symbols to search for. By default RefSeq transcripts will be identified.

=item B<-l    --list>

Text file containing a list of gene symbols to search for. Expects each symbol to be on a new line and only reads text on each line up to the first whitespace encountered. By default RefSeq transcripts will be identified.

=item B<-e    --ensembl>

Search ensembl transcripts rather than RefSeq.

=item B<-u    --ucsc>

Search UCSC known gene transcripts rather than RefSeq.

=item B<-c    --cds>

Retrieve CDS start and end coordinates rather than transcript start and end coordinates.

=item B<--merge>

Merge coordinates for overlapping transcripts.

=item B<-b    --bed>

Output in bed style format.

=item B<-d    --db>

Genome database to search. Default is hg19.

=item B<-n    --no_chr>

Don't print leading 'chr' for chromosomes.

=item B<-h    --help>

Display help message.

=item B<--manual>

Show manual page.


=back 

=cut


=head1 DESCRIPTION

For a given gene symbol(s), this program prints the genomic coordinates of the transcripts as retrieved from the UCSC genome browser. For example:

    perl gene_coordinates_from_ucsc.pl -g ADAM9 -d hg38



=cut

=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2014  David A. Parry

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

