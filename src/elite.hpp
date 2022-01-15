#pragma once

#include <vector>

#include "solution.hpp"
#include "rng.hpp"

struct Elite : public std::vector<Solution> {
  typedef std::vector<Solution> Base;

  unsigned maxsize;
  std::vector<short> freq;

  Elite(const Instance& I, unsigned _n) : maxsize(_n), freq(I.n,0) {
    assert(maxsize>0);
    reserve(_n);
  }

  bool add(const Solution& S) {
    assert(S.x.size()==freq.size());
    const unsigned n = S.x.size();

    if (size()<maxsize || S.value<back().value) {
      auto point = lower_bound(begin(),end(),S,[](const Solution& S, const Solution& T) { return S.value < T.value || (S.value == T.value && S.ones < T.ones); });
      while (point != end() && point->value == S.value && point->ones == S.ones)
	if (*point == S)
	  return false;
	else
	  point++;
      insert(point,S);
      for(unsigned i=0; i<n; i++)
	if (S.x[i])
	  freq[i]++;
      if (size()>maxsize) {
	for(unsigned i=0; i<n; i++)
	  if (back().x[i])
	    freq[i]--;
	pop_back();
      }
      assert(size()<=maxsize);
      return true;
    } else
      return false;
  }

  Solution& getRandom() {
    return (*this)[uniform_int_distribution<>(0,size()-1)(rng)];
  }
};
