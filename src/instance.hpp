#pragma once

#include <vector>
#include <iomanip>
#include <cmath>

#include <boost/multi_array.hpp>

#include "rng.hpp"

struct Instance : public boost::multi_array<double,2> {
  unsigned n;
  unsigned desks, deskDistances, colors, uncolored;
  std::vector<unsigned> deg;
  unsigned nnz;
  boost::multi_array<double,2> penaltiesColors;
  boost::multi_array<double,2> distancesDesks;
  boost::multi_array<double,2> q;
  boost::multi_array<double,2> A;
  boost::multi_array<double,2> At;
  boost::multi_array<double,2> AtA;
  vector<double> b;
  vector<double> twoDiagAtb;
  double P;
  double btb;

  void setSize(unsigned _n) {
    n = _n;
    checkPlausible();
    resize(extents[n][n]);
    q.resize(extents[n][n]);
    deg.resize(n);
    nnz = 0;
  }

  void checkPlausible() {
    if (n>10000) {
      cerr << "Instance has " << n << " variables. Refused." << endl;
      exit(1);
    }
  }

  void readInstance(istream& ins) {
    ins >> desks;
    ins >> deskDistances;
    ins >> colors;
    ins >> uncolored;

    nnz = 0;
    setSize(desks * colors);
    penaltiesColors.resize(extents[colors][colors]);
    distancesDesks.resize(extents[desks][desks]);

    vector<string> deskNames;
    for(int i = 0; i < desks; i++) {
      string deskName;
      ins >> deskName;
      deskNames.push_back(deskName);
    }

    for(int i = 0; i < desks; i++) {
      for(int j = 0; j < desks; j++) {
        distancesDesks[i][j] = 0;
      }
    }

    for(int i = 0; i < deskDistances; i++) {
      string name1, name2;
      double value;
      ins >> name1;
      ins >> name2;
      ins >> value;

      unsigned index1, index2;
      for(int j = 0; j < deskNames.size(); j++) {
        if(deskNames[j] == name1) index1 = j;
        if(deskNames[j] == name2) index2 = j;
      }

      distancesDesks[index1][index2] = value;
      distancesDesks[index2][index1] = value;
    }

    for(int i = 0; i < colors; i++) {
      for(int j = 0; j < colors; j++) {
        penaltiesColors[i][j] = 0;
      }
    }

    unsigned penalties = (colors) * (colors - 1);
    for(int i = 0; i < penalties; i++) {
      unsigned color1, color2;
      double penalty;
      ins >> color1;
      ins >> color2;
      ins >> penalty;
      
      penaltiesColors[color1][color2] = penalty;
      penaltiesColors[color2][color1] = penalty;
    }

    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        q[i][j] = 0;
      }
    }

    for(int c1 = 0; c1 < colors; c1++) {
      for(int d1 = 0; d1 < desks; d1++) {
        for(int c2 = 0; c2 < colors; c2++) {
          for(int d2 = 0; d2 < desks; d2++) {
            unsigned index1 = (colors * d1) + c1;
            unsigned index2 = (colors * d2) + c2;
            q[index1][index2] = distancesDesks[d1][d2] * penaltiesColors[c1][c2];
            deg[index1]++;
            deg[index2]++;
            if(q[index1][index2] != 0) nnz++;
          }
        }
      }
    }

    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        q[i][j] = q[i][j] / 2;
      }
    }    

    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        (*this)[i][j] = q[i][j];
      }
    }
    computeTransformation();
  }

  void computeTransformation() {
    A.resize(extents[desks + 1][n]);
    for(int d = 0; d < desks; d++) {
      int start = d * colors;
      int finish = start + colors;
      for(int i = 0; i < n; i++) {
        if(i >= start && i < finish) {
          A[d][i] = 1;
        } else {
          A[d][i] = 0;
        }
      }
    }
    int last = desks;
    int color = 0;
    for(int i = 0; i < n; i++) {
      color++;
      if(color == 1) {
        A[last][i] = 1;
      } else {
        A[last][i] = 0;
      }
      if(color == colors) {
        color = 0;
      }
    }

    At.resize(extents[n][desks + 1]);
    for(int i = 0; i < desks + 1; i++) {
      for(int j = 0; j < n; j++) {
        At[j][i] = A[i][j];
      }
    }

    AtA.resize(extents[n][n]);
    for(int i = 0; i < n; i ++) {
      for(int j = 0; j < n; j++) {
        int result = 0;
        for(int k = 0; k < desks + 1; k++) {
          result += At[i][k] * A[k][j];
        }
        AtA[i][j] = result;
      }
    }

    b.clear();
    for(int i = 0; i < desks; i++) {
      b.push_back(1);
    }
    b.push_back(uncolored);

    vector<unsigned> twoDiagAtb;
    twoDiagAtb.clear();
    for(int i = 0; i < n; i++) {
      unsigned value = 0;
      for(int j = 0; j < desks + 1; j++) {
         value += At[i][j] * b[j];
      }
      twoDiagAtb.push_back(2 * value);
    }

    P = 0;
    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        P += q[i][j];
      }
    }

    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        if(i == j) {
          AtA[i][j] = AtA[i][j] - twoDiagAtb[i];
        }
      }
    }

    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        AtA[i][j] = P * AtA[i][j];
      }
    }

    for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
        (*this)[i][j] = q[i][j] + AtA[i][j];
      }
    }

    btb = desks + (uncolored * uncolored);
  }
};
