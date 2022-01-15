#pragma once

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>

struct BTS {
  BTS() {}

  template <typename Excluder>
  unsigned bi(const Solution& S, Excluder exclude, int bkv) {
    pair<int,unsigned> best = make_pair(numeric_limits<int>::max(),S.x.size());
    unsigned nbest = 0;
    for(unsigned i=0; i<S.x.size(); i++) {
      if (S.flipvalue(i)+S.value<bkv)
      return i;
      if (exclude(i))
      continue;
      if (S.flipvalue(i)<get<0>(best)) {
        best = make_pair(S.flipvalue(i),i);
        nbest = 1;
      } else if (get<1>(best) != S.x.size() && S.flipvalue(i)==get<0>(best)) {
        nbest++;
        if (nbest*random01(rng)<=1)
        best = make_pair(S.flipvalue(i),i);
      }
    }
    return get<1>(best);
  }

  template <typename Excluder>
  unsigned step(const Solution& S, Excluder exclude, int bkv) {
    return bi(S,exclude,bkv);
  }
};

template <typename Improvement, typename Duration>
unsigned tabusearch(Solution& S, Improvement improve, Duration dgen, chrono::system_clock::time_point start, int target, const unsigned maxstagnate = numeric_limits<unsigned>::max(), const unsigned maxsteps = numeric_limits<unsigned>::max()) {
  chrono::system_clock::time_point last_report = chrono::system_clock::now();
  unsigned i, steps = 0, notimproved = 0;
  Solution B(S.I);
  B = S;

  vector<unsigned> tabu(S.x.size()+1,0);

  #if defined(EXP_STATISTICS)
  unsigned currentdistance = 0;
  Solution S0(S.I);
  S0 = S;
  accumulators::accumulator_set<double, accumulators::stats<accumulators::tag::mean,accumulators::tag::max> > acc;
  vector<unsigned> numdiff(S.x.size(),0);
  #endif

  do {
    if (chrono::system_clock::now()-last_report>chrono::duration<int>(3)) {
      last_report = chrono::system_clock::now();
      if (verbose) {
        cerr << "T " << setw(7) << steps << " " << setw(3) << chrono::duration_cast<chrono::seconds>(last_report-start).count() << " " << setw(7) << S.value << " " << setw(7) << B.value;
        #if defined(EXP_STATISTICS)
        cerr << " " << accumulators::mean(acc) << " " << accumulators::max(acc);
        #endif
        cerr << endl;
      }
    }
    if (S.value <= target)
    break;
    if (chrono::system_clock::now()-start>chrono::duration<int>(timelimit))
    break;
    if (notimproved > maxstagnate)
    break;
    if (steps > maxsteps)
    break;
    steps++;
    i = improve.step(S,[&steps,&tabu](unsigned j) { return steps<tabu[j]; },B.value);
    if (!isValid(i,S))
    continue;
    S.flip(i);
    #if defined(EXP_STATISTICS)
    if (S.x[i]!=S0.x[i]) {
      currentdistance++;
      numdiff[i]++;
    } else
    currentdistance--;
    acc(currentdistance);
    #endif

    unsigned d = dgen();
    if (d==S.x.size())
    tabu[i]=steps+S.I.deg[i]+1;
    else if (d<S.x.size())
    tabu[i]=steps+d+1;

    if (S.value < B.value) {

      B = S;
      B.time = chrono::system_clock::now();
      if (record)
      cerr << " " << setw(7) << steps << " " << setw(6) << chrono::duration_cast<chrono::milliseconds>(B.time-start).count() << " " << setw(7) << S.value << " " << setw(7) << B.value << endl;
      notimproved = 0;
    } else
    notimproved++;
  } while (true);
  steps--;
  S = B;
  return steps;
}
