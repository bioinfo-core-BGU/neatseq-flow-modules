#Script by Vered Chalifa-Caspi, Ben-Gurion University
#This script joins the per-sample output files of kaiju2table script.
#It creates one file per taxonomic level, where each file shows the data for all samples in the experiment
#Each row in the file is for one taxon_id, and the rows are sorted by total number of reads (in all samples), in descending order
#The lines for "cannot be assigned" and "unclassified" appear on top,
#and the line for viruses appear the end.

use strict;
use warnings;

# print "\n\n-->@ARGV<--\n\n";
# exit;
my @samples = split(/\,/,$ARGV[0]);
my @files = split(/\,/,$ARGV[1]);
my $out_file = $ARGV[2];
# print "\n\n---> @files\n\n";
# exit;

# my $nr_samples = scalar @samples;
# my @levels = ('phylum','class','order','family','genus','species');

my $sample;
# my $level;

# foreach $level (@levels) {
	# #data structure for tax levels data
	my $struct;
	
	#read and store sample data for this tax level
    my $in_file;
    my @samples_to_spend = @samples;
	foreach $in_file (@files) {
        $sample = shift @samples_to_spend;

		# my $in_file = "$sample/$sample.kaiju_summary_by_$level.tsv";
        print "$in_file";
		open (IN, $in_file) or die $!;
		my $heading = <IN>;
        
		while (<IN>) {
			chomp;
			# print;
			#read
			my @data = split (/\t/);
			my ($file, $percent, $reads, $taxon_id, $taxon_path) = split (/\t/);
			
            # print "\n\n$file, $percent, $reads, $taxon_id, $taxon_path\n\n";
            # exit;
			#treat tax id of "cannot be assigned" and "unclassified"
			if ($taxon_id == 0) {
				if ($taxon_path =~ /^cannot be assigned/) {
					$taxon_id = -1;
				} elsif ($taxon_path eq 'unclassified') {
					$taxon_id = -2;
				}
			}
			
			#store in data structure
			$struct->{$taxon_id}->{'taxon_path'} = $taxon_path;
			$struct->{$taxon_id}->{$sample}->{'percent'} = $percent;
			$struct->{$taxon_id}->{$sample}->{'reads'} = $reads;			
		}		
		close (IN);
	}
	# use Data::Dumper;
	# print Dumper $struct; exit;
	# my $out_file = "kaiju.summary_by_$level.tsv";
	open (OUT, ">$out_file") or die $!;

	#generate and print headings line for summary table
	my @headings = ('taxon_id', 'taxon_path', 'total_reads');
	foreach $sample (@samples) {
        print "$sample\n";
		push (@headings, "reads $sample");
	}
	foreach $sample (@samples) {
		push (@headings, "percent $sample");
	}
	print OUT join ("\t", @headings), "\n";
	
	#summarize data
	my $struct1;
	my $unclassified;
	my $unassigned;
	my $viruses;
	my @data;
	my $taxon_id;
	
	#foreach $taxon_id (sort {$a <=> $b} keys %{$struct}) {
	foreach $taxon_id (keys %{$struct}) {
		@data = ($taxon_id, $struct->{$taxon_id}->{'taxon_path'});
		my $total_reads = 0;
		my @reads_data;
		foreach $sample (@samples) {
			my $reads = 0;
			if ($struct->{$taxon_id}->{$sample}->{'reads'}) {
				$reads = $struct->{$taxon_id}->{$sample}->{'reads'};
			}
			push (@reads_data, $reads);
			$total_reads += $reads;
		}
		push (@data, $total_reads, @reads_data);
		foreach $sample (@samples) {
			my $percent = 0;
			if ($struct->{$taxon_id}->{$sample}->{'percent'}) {
				$percent = $struct->{$taxon_id}->{$sample}->{'percent'};
			}
			push (@data, $percent);
		}
		#exclude unassigned, unclassified and viruses from the sort, and then create a data structure for sorting the rest of the tax IDs
		if ($taxon_id == -1) {
			$data[0] = "";  #don't print tax ID as "-1"
			$unassigned = [@data];
		} elsif ($taxon_id == -2) {
			$data[0] = "";  #don't print tax ID as "-2"
			$unclassified = [@data];
		} elsif ($taxon_id == 10239) {
			$viruses = [@data];
		} else {
			$struct1->{$taxon_id} = [@data];
		}		
	}
	
	#sort and print data
	print OUT join ("\t", @{$unclassified}), "\n";
	print OUT join ("\t", @{$unassigned}), "\n";
	my $tax_id;
	foreach $tax_id (sort {$struct1->{$b}->[2] <=> $struct1->{$a}->[2]} keys %{$struct1}) {
		print OUT join ("\t", @{$struct1->{$tax_id}}), "\n";
	}
	print OUT join ("\t", @{$viruses}), "\n";	
	close (OUT);
	


	
print "\nDone\n";

#sub by_total_reads {
#	$struct1->{b}->[2] <=> $struct1->{a}->[2];
#}