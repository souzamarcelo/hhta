#pragma once

#include "bqp.hpp"
#include "elite.hpp"

namespace recombine {

  struct Fn {
    unsigned i;
    bool rr;
    Fn(bool _rr) : i(0), rr(_rr) {}

    unsigned operator()(const Solution& S, const vector<unsigned>& NC) {
      pair<int,unsigned> best = make_pair(numeric_limits<int>::max(),S.x.size()); // delta, index
      unsigned nbest = 0;
      unsigned ncand = NC.size();
      i = i%NC.size();
      while (ncand-->0 && S.value <= get<0>(best)) {
        if (S.flipvalue(NC[i])<get<0>(best)) {
          best = make_pair(S.flipvalue(NC[i]),i);
          nbest = 1;
        } else if (S.flipvalue(NC[i])==get<0>(best)) {
          assert(get<1>(best) != S.x.size());
          nbest++;
          if ((nbest)*random01(rng)<=1.0)
          best = make_pair(S.flipvalue(NC[i]),i);
        }
        i = (i+1)%NC.size();
      }

      if(!rr) i = 0;

      return get<1>(best);
    }
  };

  template<typename Improver>
  void recombine(const Instance& I, Solution& S, const Solution& T, Improver improve, double gamma) {
    const unsigned n = S.x.size();
    assert(n==T.x.size());
    assert(0<=gamma && gamma<0.5);

    Solution B(S.I);
    B.value = numeric_limits<int>::max();

    vector<unsigned> NC;
    for(unsigned i=0; i<n; i++)
      if (S.x[i]!=T.x[i])
    NC.push_back(i);

    unsigned dmin = gamma*NC.size();
    unsigned dmax = NC.size()-dmin;
    
    unsigned d = 0;
    while (NC.size() > 0) {
      unsigned i = improve(S,NC);
      S.flip(NC[i]);
      d++;

      swap(NC[i],NC.back());
      NC.pop_back();

      if (S.value < B.value && dmin <= d && d <= dmax) {
        B = S;
        B.time = chrono::system_clock::now();
      }
    }

    S = B;
  }

  template <typename Improver, typename Recombiner>
  void recombiner(const Instance& I, Solution& S, unsigned b, Recombiner recombine, Improver improve, chrono::system_clock::time_point start, int target) {
    chrono::system_clock::time_point last_report = chrono::system_clock::now();
    unsigned steps = 0;
    vector<Solution> current;
    Solution B(S.I);
    B = S;

    Elite e(I,b);

    while (true) {
      e.clear();
      for(unsigned i=0; i<b; i++) {
        Solution S(I);
        improve(S);
        if (verbose)
          cout << "R" << i << " " << setw(7) << S.value << endl;
        e.add(S);
        if (termination(e[0].value,e.back().value,steps,target,start,last_report))
          goto done;
      }

      current = e;
      for(Solution& s : e)
      s.novel = false;

      unsigned num_novel = b;

      while (num_novel>0) {
        num_novel = 0;

        for(unsigned i=0; i<current.size(); i++)
        for(unsigned j=i+1; j<current.size(); j++) {
          if (!current[i].novel && !current[j].novel)
          continue;
          steps++;

          Solution S = current[i];
          recombine(I,S,current[j]);
          improve(S);
          S.novel=true;
          if (e.add(S))
          num_novel++;

          if (termination(e[0].value,e.back().value,steps,target,start,last_report))
          goto done;

          S = current[j];
          recombine(I,S,current[i]);
          improve(S);
          S.novel=true;
          if (e.add(S))
          num_novel++;

          if (termination(e[0].value,e.back().value,steps,target,start,last_report))
          goto done;
        }

        current = e;
        for(Solution& s : e)
        s.novel = false;
      }
      if (e[0].value<B.value)
      B = e[0];
    }
    done:
    if (e[0].value<B.value)
    B = e[0];
    S = B;
  }

}