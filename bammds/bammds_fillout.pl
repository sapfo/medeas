#!/usr/bin/perl

# Copyright (C) 2013,2014 Ole Tange, Mike DeGiorgio, Anna-Sapfo
# Malaspinas, Jose Victor Moreno-Mayar, Yong Wang and Free Software
# Foundation, Inc.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License,
# or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

use local::lib;
use Text::CSV;

my $file = shift;
# Assume we can guess the CSV format from the first 100KB.
my $sniff_length = 100_000;

my @rows;

my @csv_file_types = 
  ( { binary => 1 },
    { binary => 1, sep_char => "\t" },
    { binary => 1, sep_char => "," },
  );

#    $csv = Text::CSV->new ({
#            quote_char          => '"',
#            escape_char         => '"',
#                        => ',',
#            eol                 => $\,
#            always_quote        => 0,
#            quote_space         => 1,
#            quote_null          => 1,
#            binary              => 0,
#            keep_meta_info      => 0,
#            allow_loose_quotes  => 0,
#            allow_loose_escapes => 0,
#            allow_whitespace    => 0,
#            blank_is_undef      => 0,
#            empty_is_undef      => 0,
#            verbatim            => 0,
#            auto_diag           => 0,
#            });

open(my $fh, "<:encoding(utf8)", $file) or die "Cannot open $file: $!";
my @lines = <$fh>;
close $fh;

my $succesful_csv_type;
my $csv;
for my $csv_file_type (@csv_file_types) {
  $csv = Text::CSV->new ( $csv_file_type )  # should set binary attribute.
    or die "Cannot use CSV: ($csv_file_type) ".Text::CSV->error_diag ();
  $succesful_csv_type = $csv_file_type;
  for my $line (@lines) {
    if(not $csv->parse($line)) {
      $succesful_csv_type = 0;
      last;
    }
  }
}
if(not $succesful_csv_type) {
  $csv->error_diag();
}

$csv = Text::CSV->new ( $succesful_csv_type )  # should set binary attribute.
  or die "Cannot use CSV: ".Text::CSV->error_diag ();
my @table;
my $max_columns;
my $i = 0;
for my $line (@lines) {
  $csv->parse($line);
  my @columns = $csv->fields();
  for(my $j=0; $j<= $#columns; $j++) {
    $table[$i][$j] = $columns[$j];
  }
  $max_columns = max($max_columns,$#columns);
  $i++;
}

# Remove lines starting with '#'
@table = grep { $$_[0] !~ /^#/ } @table;
# Remove empty lines
@table = grep { grep { /\S/ } @$_ and $_  } @table;

# Empty value = see above (Get value from previous line)
my @above = @{$table[0]};
for(my $i = 0; $i <= $#table; $i++) {
  for(my $j = 0; $j <= $max_columns; $j++) {
    if($table[$i][$j] =~ /^\s*$/) {
      $table[$i][$j] = $above[$j];
    }
  }
  @above = @{$table[$i]};
}

my @predefined_colnames = qw(Population pop_label Individual indv_label sample order Color pch cex);
my @colnames = @{$table[0]};
my %col;
@col{@colnames} = 0..$#colnames;

# Get the columns number of the predefined names
my($Population, $pop_label, $Individual, $indv_label, $sample, $order, $Color, $pch, $cex) =
  @col{@predefined_colnames};

# Find population default values
for my $row (@table) {
  if($row->[$Individual] eq "*") {
    my %values;
    @values{@predefined_colnames} =  @$row[@col{@predefined_colnames}];
    $default{$row->[$Population]} = \%values;
  }
}

$default->{'*'}{'cex'} = 0.8;
$default->{'*'}{'sample'} = 0;

# Replace population default values
for my $row (@table) {
  if($row->[$Individual] eq "*") {
    # This is a population line
    if($row->[$pop_label] eq "*") {
      # pop_label == *  => pop_label = Population
      $default->{$row->[$Population]}{'pop_label'} = $row->[$Population];
      $row->[$pop_label] = $default->{$row->[$Population]}{'pop_label'};
    } else {
      $default->{$row->[$Population]}{'pop_label'} = $row->[$pop_label];
    }

    if($row->[$sample] eq "*") {
      # sample == * => sample = sample[*]
      $default->{$row->[$Population]}{'sample'} = $default->{'*'}{'sample'};
      $row->[$sample] = $default->{$row->[$Population]}{'sample'};
    } else {
      $default->{$row->[$Population]}{'sample'} = $row->[$sample];
    }

    if($row->[$Color] eq "*") {
      # Color == * 
      # Skip: Will be done later
    }

    if($row->[$pch] eq "*") {
      # pch == * => pch = FirstChar of pop_label
      $default->{$row->[$Population]}{'pch'} = substr($row->[$pop_label],0,1);
      $row->[$pch] = $default->{$row->[$Population]}{'pch'};
      if($row->[$sample] eq 1) {
	# Sample => Default to [n]
	$row->[$pch] = "[".$row->[$pch]."]";
      }
    } else {
      $default->{$row->[$Population]}{'pch'} = $row->[$pch];
    }
    
    if($row->[$cex] eq "*") {
      # cex == * => cex = cex[*];
      $default->{$row->[$Population]}{'cex'} = $default->{'*'}{'cex'};
      $row->[$cex] = $default->{$row->[$Population]}{'cex'};
    } else {
      $default->{$row->[$Population]}{'cex'} = $row->[$cex];
    }
  }
}

# Replace values
my ($uniq_pop, $uniq_individual);

for my $row (@table) {
  # Population = * => Done in finding default values

  if($row->[$pop_label] eq "*") {
    # pop_label == *  => pop_label[Population]
    $row->[$pop_label] = $default->{$row->[$Population]}{'pop_label'};
  }
  $uniq_pop{$row->[$pop_label]}++;

  # Individual = * => Done in finding default values

  if($row->[$indv_label] eq "*") {
    # indv_label == *  => indv_label = Individual
    $row->[$indv_label] = $row->[$Individual];
  }
  $uniq_individual{$row->[$indv_label]}++;
}

my (%known_position,$pop_id,$indv_id);
my $i = 1;
for my $row (@table) {
  if($row->[$sample] eq "*") {
    # sample == * => sample = sample[Population]
    $row->[$sample] = $default->{$row->[$Population]}{'sample'} || $default->{'*'}{'sample'};
  }

  if($row->[$order] eq "*") {
    # order == *  => order = lineno
    $row->[$order] = $i*10;
  }

  if($row->[$Color] eq "*") {
    # Color == * => color = hash(indv_label or pop_label)
    my ($name,$position);
    if($row->[$Individual] eq "*") {
      # This is a population line
      $name = $row->[$pop_label];
      if(not $known_position{$name}) {
	  # If we have not seen this name before: Get new color
	  $known_position{$name} = $pop_id++/(keys %uniq_pop);
      }
      $position = $known_position{$name};
    } else {
      # This is an individual line
      $name = $row->[$indv_label];
      $position = $indv_id++/(keys %uniq_individual);
    }
    if($row->[$sample] eq 1) {
      # This is a sample. We make that black
      $row->[$Color] = color_from_position(0,0,0);
    } else {
      # Not a sample -> pastel color
      $row->[$Color] = color_from_position($position,0.6,1);
    }
  }

  if($row->[$pch] eq "*") {
    if(not $row->[$Individual] eq "*") {
      # This is an individual line
      # pch == * => pch = hash(indv_label);
      $row->[$pch] = pch_from_name($row->[$indv_label]);
      if($row->[$sample] eq 1) {
	# Sample => Default to [n]
	$row->[$pch] = "[".$row->[$pch]."]";
      }
    }
  }

  if($row->[$cex] eq "*") {
    # cex == * => cex = cex[Population];
    $row->[$cex] = $default->{$row->[$Population]}{'cex'} || $default->{'*'}{'cex'};
  }
  $i++;
}


$csv = Text::CSV->new ({ binary => 1 })  # should set binary attribute.
  or die "Cannot use CSV: ".Text::CSV->error_diag ();
for my $row (@table) {
  my @columns = @$row[@col{@predefined_colnames}];
  $status = $csv->combine(@columns);    # combine columns into a string
  $line   = $csv->string();             # get the combined string
  print $line,"\n";
}

sub color_from_position {
  # Compute the HSV color for the value
  # Input:
  #   h = 0.0 .. 1.0
  #   s = 0.0 .. 1.0
  #   v = 0.0 .. 1.0
  # Output:
  #   $hex = "#xxxxx" hex color value
  my $h = shift;
  my $s = shift;
  my $v = shift;
  my ($r,$g,$b) = hsv2rgb($h*360,$s,$v);
  return sprintf('#%02x%02x%02x', $r*255, $g*255, $b*255);
}  

sub hsv2rgb {
  my ( $h, $s, $v ) = @_;

  use POSIX;
  
  if ( $s == 0 ) {
    return $v, $v, $v;
  }
  
  $h /= 60;
  my $i = floor( $h );
  my $f = $h - $i;
  my $p = $v * ( 1 - $s );
  my $q = $v * ( 1 - $s * $f );
  my $t = $v * ( 1 - $s * ( 1 - $f ) );
  
  if ( $i == 0 ) {
    return $v, $t, $p;
  }
  elsif ( $i == 1 ) {
    return $q, $v, $t;
  }
  elsif ( $i == 2 ) {
    return $p, $v, $t;
  }
  elsif ( $i == 3 ) {
    return $p, $q, $v;
  }
  elsif ( $i == 4 ) {
    return $t, $p, $v;
  }
  else {
    return $v, $p, $q;
  }
}

sub color_from_name {
  # Input:
  #   $name = name to give a color
  # Output:
  #   $hex = "#xxxxx" hex color value
  my $name = shift;
  use Digest::MD5 qw(md5);
  my $str = substr( md5($name), 0, 3 );
  my ($b_R,$b_G,$b_B) = unpack('CCC', $str);
  return sprintf('#%02x%02x%02x', $b_R, $b_G, $b_B);
}


sub pch_from_name {
  my $name = shift;
  my @pch_values=("A".."Z","a".."z",0..9);
  use Digest::MD5 qw(md5 md5_hex);
  my ($val64bit) = unpack('Q', md5($name));
  my $mod = $val64bit % ($#pch_values+1);
  return $pch_values[$mod];
}

sub my_dump {
    # Returns:
    #   ascii expression of object if Data::Dump(er) is installed
    #   error code otherwise
    my @dump_this = (@_);
    eval "use Data::Dump qw(dump);";
    if ($@) {
        # Data::Dump not installed
        eval "use Data::Dumper;";
        if ($@) {
            my $err =  "Neither Data::Dump nor Data::Dumper is installed\n".
                "Not dumping output\n";
            print $Global::original_stderr $err;
            return $err;
        } else {
            return Dumper(@dump_this);
        }
    } else {
	# Create a dummy Data::Dump:dump as Hans Schou sometimes has
	# it undefined
	eval "sub Data::Dump:dump {}";
        eval "use Data::Dump qw(dump);";
        return (Data::Dump::dump(@dump_this));
    }
}

sub max {
    # Returns:
    #   Maximum value of array
    my $max;
    for (@_) {
        # Skip undefs
        defined $_ or next;
        defined $max or do { $max = $_; next; }; # Set $_ to the first non-undef
        $max = ($max > $_) ? $max : $_;
    }
    return $max;
}

