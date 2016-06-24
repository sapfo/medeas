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

typedef vector<int> Row;
typedef vector<Row> Table;


void output(int num_indv, Table asd, Table count, Table marker) {
  vector<int> keep(num_indv);
  // Print (sharing distance, count) for each pair of individuals
  for(int i = 0; i < num_indv; i++) {
    for(int j = 0; j < num_indv; j++) {
      if(count[i][j] > 0) {
	cout << asd[i][j] << " " << count[i][j] << " ";
      } else {
	cout << -9 << " " << count[i][j] << " ";
      }
    }
    cout << endl;
  }
  // Print (count marker, 0) for each marker and individual
  // The 0 is needed to maintain 2*indv columns as above
  for(unsigned int m = 0; m < marker[0].size(); m++) {
    for(int i = 0; i < num_indv; i++) {
      cout << marker[i][m] << " ";
      cout << 0 << " ";
    }
    cout << endl;
  }
}


Row readline() {
  // Read until EOL
  string stdin_line;
  getline(cin, stdin_line);
  istringstream is(stdin_line);
  Row row;
  string ASD_or_count;
  // Read the string 'ASD ' or 'Count ' and ignore it.
  is >> ASD_or_count;
  while(is) {
    int data;
    is >> data;
    row.push_back(data);
  }

  return(row);
}


Table unpack_table (int num_indv, Row r) {
  // convert from a list representing a triangle to a full matrix
  Table t;
  for(int i = 0; i < num_indv; i++) {
    Row empty(num_indv);
    t.push_back(empty);
  }
  for(int i = 0; i < num_indv; i++) {
    Row empty(num_indv);
    t.push_back(empty);
    for(int j = i; j < num_indv; j++) {
      int packedidx = (num_indv*(num_indv-1) - (num_indv-i)*((num_indv-i)-1))/2 + j;
      t[i][j] = r[packedidx];
      t[j][i] = r[packedidx];
    }
  }
  return(t);
}


int main() {
  // Input format:
  //    num_indv
  //    ASD n n n n n ...
  //    Count n n n n n ...
  //    Marker n*num_indv
  //    Marker n*num_indv
  //    ...
  //    END
  

  string stdin_line;
  getline(cin, stdin_line);
  if(!cin) {
    // no input - just stop
    exit(0);
  }
  istringstream is(stdin_line);
  int num_indv;
  // Read: num_indv
  is >> num_indv;

  Row sum_asd,sum_count;
  Row asd,count;
  string asd_count;
  // Read:   ASD n n n n n ...
  sum_asd = readline();
  // Read:   Count n n n n n ...
  sum_count = readline();

  Table marker;
  for(int i = 0; i < num_indv; i++) {
    Row empty;
    marker.push_back(empty);
  }

  string line;
  string end = "END";
  while (getline(cin, line)) {
    if(line == end) { break; }
    // Read a marker line
    // {0,1} for each indv
    istringstream is(line);
    string marker_string;
    // Ignore the header: Marker
    is >> marker_string;
    for(int i = 0; i < num_indv; i++) {
      int data;
      is >> data;
      marker[i].push_back(data);
    }
  }

  while(cin) {
    // Next record
    getline(cin, stdin_line);
    if(!cin) {
      // no input - just stop
      break;
    }
    istringstream is(stdin_line);
    is >> num_indv;
    asd = readline();
    for(unsigned int i = 0; i < asd.size(); i++) {
      sum_asd[i] += asd[i];
    }
    count = readline();
    for(unsigned int j = 0; j < count.size(); j++) {
      sum_count[j] += count[j];
    }
    string line;
    string end = "END";
    while (getline(cin, line)) {
      if(line == end) { break; }
      // Read a marker line
      // {0,1} for each indv
      istringstream is(line);
      string marker_string;
      // Ignore the header: Marker
      is >> marker_string;
      for(int i = 0; i < num_indv; i++) {
	int data;
	is >> data;
	marker[i].push_back(data);
      }
    }
  }

  Table asd_unpacked = unpack_table(num_indv,sum_asd);
  Table count_unpacked = unpack_table(num_indv,sum_count);
  output(num_indv, asd_unpacked, count_unpacked, marker);
}
