#!/usr/bin/perl

use strict;
use warnings;

use SWISS::Entry;
use SWISS::OCs;
use SWISS::GNs;

# usage: -i uniprot_sprot.dat -t swissTable.txt -f swissProtein.fasta -d discard_directory
my $usage = "./extractSwiss -i <.dat file> -t <output table> ";

# input file, .dat for swissprot data, including annotation and a lot of other information
my $input = "uniprot_sprot.dat";

# output files:
my $table = "swissTable.txt";

#foreach my $arg (@ARGV){
#    print "$arg\t";
#}

# Parse the command line arguments.
for(my $i = 0; $i < @ARGV ; $i++) {
    if($ARGV[$i] =~ /^\-i/) { $input = $ARGV[$i+1]; $i++; }
    elsif($ARGV[$i] =~ /^\-t/) { $table = $ARGV[$i+1]; $i++; }
    elsif($ARGV[$i] =~ /^\-h/) { die("Usage:\n $usage \n"); }
}


#die("end here\n");

# echo input arguments:
print "input dat file: $input\n";
print "output table file: $table\n";

# Open files for reading and writing
open INPUT, "<$input" or die("couldn't open $input, $! \n");
open TABLE, ">$table" or die("couldn't open $table, $! \n");


# set separator
# proteins records are separated by \n//\n instead of newline
local $/ = "\n//\n";

#print TABLE "ac;\tOC0;\tOC1;\tRecNameFull;\tRecNameEC;\tGeneName;\tFastaPosition;\tProteinLength;\tGO\n";
## use formatted printf for readability
printf TABLE "%-8s  \t  %-15s  \t %-40s \t %-40s \t    %-40s \t  %-40s \t  %-40s \t  %-70s \t %-100s \t %-120s \n",
             "AC", "Organism",  "Species", "RecNameEC","GeneName", "OLN", "ORF", "GO", "Keywords", "RecNameFull";


# Read an entire record at a time
while (<INPUT>){
    # Read the entry, using SwissKnife module, $_: default or implicit variable
    my $entry = SWISS::Entry->fromText($_);
    

    #next if $entry->isFragment;     # skip entries that are fragments
    
    my $ac = "NA"; # accession number
    my $OCAll = "NA";
    my $OSAll = "NA";
    my @RecNameFull = ();
    my $RecNameFullAll = "NA"; # recommended full names
    my @RecNameEC = ();
    my $RecNameECAll = "NA"; # recommendec enzyme commission number
    my @GeneName = ();
    my $GeneNameAll = "NA"; # gene names
    my @OLN = ();
    my $OLNAll = "NA";
    my @ORF = ();
    my $ORFAll = "NA";
    my @KW = ();
    my $KWAll = "NA";
    my $GoTerms = "NA"; # go terms
    
    # get the primary ac of each entry and the sequence
    $ac = $entry->AC;
    
    ############# get GO terms from DRs #############
    my @terms = ();
    foreach my $DR ($entry->DRs->list ){
        my @DRarray = @$DR;
        foreach my $xref (@DRarray){
            my @itemArray = @$xref;
            if($itemArray[0] eq "GO" ){
                $itemArray[1] =~ /GO:(\d+)/;
                push(@terms, $1);
            }
        }
        
    }
    
    if(scalar(@terms) > 0 ){
        $GoTerms = join(':', @terms);
    }

    ############# get RecNames from DEs #############
    foreach my $de ($entry->DEs->elements) {
        if( $de->category eq "RecName" && $de->type eq  "Full")
        {
            push(@RecNameFull, $de->text);
        }
        if( $de->category eq "RecName" && $de->type eq  "EC")
        {
            push(@RecNameEC, $de->text);
        }
    }
    # Includes: protein contains multiple functional domains
    my $includes = $entry->DEs->Includes;
    # Contains: protein cleaved into multiple functional components
    my $contains = $entry->DEs->Contains;

    if($includes->size > 0){
        foreach my $include ($includes->elements){ # $inclue : $DEs
            foreach my $element ($include->elements){
                if( $element->category eq "RecName" && $element->type eq  "Full")
                {
                    push(@RecNameFull, $element->text);
                }
                if( $element->category eq "RecName" && $element->type eq  "EC")
                {
                    push(@RecNameEC, $element->text);
                }
            }
            
        }
    }
    
    if($contains->size > 0){
        foreach my $include ($contains->elements){ # $inclue : $DEs
            foreach my $element ($include->elements){
                if( $element->category eq "RecName" && $element->type eq  "Full")
                {
                    push(@RecNameFull, $element->text);
                }
                if( $element->category eq "RecName" && $element->type eq  "EC")
                {
                    push(@RecNameEC, $element->text);
                }
            }
            
        }
    }

    if(scalar(@RecNameFull) > 0 ){
        $RecNameFullAll = join(':', @RecNameFull);
    }
    if(scalar(@RecNameEC) > 0 ){
        $RecNameECAll = join(':', @RecNameEC);
    }

    ############# get gene Names, OLNs and ORFs from GNs #############
    foreach my $geneGroup ($entry->GNs->elements) {
        
        foreach my $gn ($geneGroup->Names->elements) {
            push(@GeneName, $gn->text);
        }
        foreach my $gn ($geneGroup->OLN->elements) {
            push(@OLN, $gn->text);
        }
        foreach my $gn ($geneGroup->ORFNames->elements) {
            push(@ORF, $gn->text);
        }
    }
    if(scalar(@GeneName) > 0 ){
        $GeneNameAll = join(':', @GeneName);
    }
    if(scalar(@OLN) > 0 ){
        $OLNAll = join(':', @OLN);
    }
    if(scalar(@ORF) > 0 ){
        $ORFAll = join(':', @ORF);
    }
    
    ############# get keywords from KWs #############
    foreach my $keywords ($entry->KWs->elements) {
#        print $keywords->text,"\n";
        push(@KW, $keywords->text);

    }
    if(scalar(@KW) > 0 ){
        $KWAll = join(':', @KW);
    }

    ############# get organism info from KWs #############
    foreach my $organism ($entry->OCs->elements) {
        $OCAll = $organism;
        last;
        
    }
    
    foreach my $species ($entry->OSs->elements){
#        print $species->text, "\n";
        $OSAll = $species->text;
        last;
    }
#    foreach my $species ($entry->OS->elements) {
#        $OCAll = $organism;
#        last;
#        
#    }
    printf TABLE "%-8s  \t  %-15s  \t %-40s \t %-40s \t    %-40s \t  %-40s \t  %-40s \t  %-70s \t %-100s \t %-120s \n",
                $ac,$OCAll, $OSAll, $RecNameECAll,$GeneNameAll, $OLNAll, $ORFAll, $GoTerms, $KWAll, $RecNameFullAll;

}

close TABLE;
close INPUT;




