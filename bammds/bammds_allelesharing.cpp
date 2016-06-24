/*
 Copyright (C) 2013,2014 Ole Tange, Mike DeGiorgio, Anna-Sapfo
 Malaspinas, Jose Victor Moreno-Mayar, Yong Wang and Free Software
 Foundation, Inc.
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License,
 or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cstdlib>

using namespace std;

typedef INNERINT innerint;

void _output(int num_indv, vector<int> asd, vector<int> count) {
  cout << num_indv << endl;
  cout << "ASD ";
  for(int j = 0; j < num_indv; j++) {
    for(int k = 0; k < num_indv; k++) {
      int idx=j*num_indv+k;
      cout << asd[idx] << " ";
    }
  }
  cout << endl;
  
  cout << "Count ";
  for(int j = 0; j < num_indv; j++) {
    for(int k = 0; k < num_indv; k++) {
      int idx=j*num_indv+k;
      cout << count[idx] << " ";
    }
  }
  cout << endl;
}

void packed_output(int num_indv, vector<int> asd, vector<int> count) {
  cout << num_indv << endl;
  cout << "ASD   ";
  for(int j = 0; j < num_indv; j++) {
    for(int k = j; k < num_indv; k++) {
      int packedidx = (num_indv*(num_indv-1) - (num_indv-j)*((num_indv-j)-1))/2 + k;
      cout << asd[packedidx] << " ";
    }
  }
  cout << endl;

  cout << "Count ";
  for(int j = 0; j < num_indv; j++) {
    for(int k = j; k < num_indv; k++) {
      int packedidx = (num_indv*(num_indv-1) - (num_indv-j)*((num_indv-j)-1))/2 + k;
      cout << count[packedidx] << " ";
    }
  }
  cout << endl;
}

void output_marker(vector< vector<bool> > markertable) {
  for(unsigned int i = 0; i < markertable.size(); i++) {
    cout << "Marker ";
    for(unsigned int j = 0; j < markertable[i].size(); j++) {
      cout << markertable[i][j] << " ";
    }
    cout << endl;
  }
}


vector<int> readline() {
    vector<int> alleles;

    std::string line;
    std::getline(cin, line);
    std::istringstream is(line);
    
    string dummy;
    int t=0;
    while (is) {
	// Ignore the first 5 columns
	if(t++ < 5) {
	    is >> dummy;
	} else {
	    int data;
	    is >> data;
	    alleles.push_back(data);
	}
    }
    return(alleles);
}


vector<bool> process_line(int num_indv, vector<int> &alleles, vector<innerint> &asd, vector<innerint> &missing) {
  vector<bool> marker_defined;
  for(int j = 0; j < num_indv; j++) {
    for(int k = j; k < num_indv; k++) {
      int idx = (num_indv*(num_indv-1) - (num_indv-j)*((num_indv-j)-1))/2 + k;
      // If alleles are not missing
      if(alleles[j] + alleles[k] >= 0) {
	// If the alleles are different: absolute value of difference
	if(alleles[j] - alleles[k] != 0) {
	  asd[idx]++;
	}
      } else {
	missing[idx]++;
      }
    }

    if(alleles[j] == -9) {
      // This marker is undefined
      marker_defined.push_back(false);
    } else {
      marker_defined.push_back(true);
    }
  }
  return(marker_defined);
}


int main() {
    int num_indv;
    vector<int> alleles;
    // Read the first line to get the num_indv
    alleles = readline();
    if(!cin) {
      // no input - just stop
      exit(0);
    }
    num_indv = alleles.size() - 1;

    // Initialize asd + count + missing matrix
    // innerint is a smaller type than int, is faster to use,
    // but can only contain MAX as value.
    vector<innerint> asd(num_indv*num_indv);
    vector<innerint> missing(num_indv*num_indv);
    vector<int> asd_int(num_indv*num_indv);
    vector<int> missing_int(num_indv*num_indv);
    vector<int> count_int(num_indv*num_indv);
    vector< vector<bool> > markertable;

    // Process the first line
    markertable.push_back(process_line(num_indv,alleles,asd,missing));

    int lineno=1;
    // Read and process the rest of the file
    while (cin) {
      alleles = readline();
      if(!cin) {
	// We read the last line
	break;
      }
      markertable.push_back(process_line(num_indv,alleles,asd,missing));

      if(!(lineno%MAX)) {
	// A char can contain up to 127, so empty the shorts into ints
	// A uchar can contain up to 255, so empty the shorts into ints
	// A short can contain up to 32767, so empty the shorts into ints
	// A ushort can contain up to 65535, so empty the shorts into ints
	// A int can contain up to 2.1b, so this will never happen
	for(int i = 0; i < num_indv*num_indv; i++) {
	  asd_int[i] += asd[i];
	  asd[i] = 0;
	  missing_int[i] += missing[i];
	  missing[i] = 0;
	}
      }
      lineno++;
    }

    // Empty the innerints into ints
    for(int i = 0; i < num_indv*num_indv; i++) {
      asd_int[i] += asd[i];
      asd[i] = 0;
      missing_int[i] += missing[i];
      missing[i] = 0;
    }

    // Compute count_int
    for(int i = 0; i < num_indv*num_indv; i++) {
      count_int[i] = lineno-missing_int[i];
    }

    // Copy packed upper triangle to full matrix
    vector<int> asd_print(num_indv*num_indv);
    vector<int> count_print(num_indv*num_indv);
    for(int j = 0; j < num_indv; j++) {
      for(int k = j; k < num_indv; k++) {
	int packedidx = (num_indv*(num_indv-1) - (num_indv-j)*((num_indv-j)-1))/2 + k;
	int idx = j*num_indv+k;       
	int symidx = k*num_indv+j;       
	asd_print[idx] = asd_int[packedidx];
	count_print[idx] = count_int[packedidx];
	asd_print[symidx] = asd_int[packedidx];
	count_print[symidx] = count_int[packedidx];
      }
    }

    // Print asd and count matrices
    // output(num_indv, asd_print,count_print);
    packed_output(num_indv, asd_int, count_int);

    output_marker(markertable);
    cout << "END" << endl;
    exit(0);
}
