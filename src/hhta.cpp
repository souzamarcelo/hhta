#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <string>
#include <cassert>
#include <algorithm>
#include <unordered_set>
#include <tuple>
#include <chrono>
#include <iterator>
#include <iomanip>
#include <boost/multi_array.hpp>
#include <boost/program_options.hpp>
#include <sys/ioctl.h>

namespace po = boost::program_options;
using namespace boost;
using namespace std;

#include "instance.hpp"
#include "solution.hpp"
#include "rng.hpp"
#include "rng.cpp"
#include "bqp.hpp"
#include "bqp.cpp"
#include "tabusearch.hpp"
#include "recombine.hpp"

Instance I;
chrono::system_clock::time_point start_time;
int target;

unsigned setupRandom(unsigned seed) {
  if (seed == 0) {
    seed = time(0);
    ifstream f("/dev/urandom");
    if (f.good()) {
      f.read((char *)(&seed), sizeof(unsigned int));
    }
  }
  rng.seed(seed);
  srand48(seed);
  return seed;
}

void hhta1(Solution& S, chrono::system_clock::time_point start_time, unsigned td, unsigned sm, unsigned iv, unsigned e, double gamma) {
  recombine::Fn fn(true);
  unsigned maxstagnate = sm * I.n;
  unsigned maxsteps = iv;
  BTS bt;

  recombine::recombiner(I, S, e, 
                        [&](const Instance& I, Solution& S, const Solution& T) { return recombine::recombine(I, S, T,fn, gamma); },
                        [&](Solution& S) { return tabusearch(S, bt, [&]() { return I.n/td; }, start_time, target, maxstagnate, maxsteps);},
                        start_time, target);
}

void hhta2(Solution& S, chrono::system_clock::time_point start_time, unsigned sm, unsigned e, double gamma) {
  recombine::Fn fn(false);
  unsigned maxstagnate = sm * I.n;
  unsigned maxsteps = I.n>5000?15000:(I.n>3000?12000:10000);
  BTS bt;

  recombine::recombiner(I, S, e,
                        [&](const Instance& I, Solution& S, const Solution& T) { return recombine::recombine(I, S, T,fn, gamma); },
                        [&](Solution& S) { return tabusearch(S, bt, [&]() { return I.n; }, start_time, target, maxstagnate, maxsteps);},
                        start_time, target);

}

int main(int argc, char *argv[]) {

  struct winsize wsize;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &wsize);

  po::options_description desc("General options",wsize.ws_col);
  desc.add_options()
    ("help", "show help")
    ("ins", po::value<string>(), "instance")
    ("seed", po::value<unsigned>()->default_value(0), "random seed")
    ("target", po::value<int>()->default_value(numeric_limits<int>::min()), "target value")
    ("timelimit", po::value<unsigned>()->default_value(10), "time limit (s)")
    ("variant", po::value<unsigned>()->default_value(1), "1: HHTA1; 2: HHTA2")
    ("td", po::value<unsigned>()->default_value(1), "denominator for tabu tenure (HHTA1 only)")
    ("sm", po::value<unsigned>()->default_value(81), "multiplier for max stagnation of tabu search")
    ("iv", po::value<unsigned>()->default_value(4204), "max steps of tabu search (HHTA1 only)")
    ("e", po::value<unsigned>()->default_value(20), "elite set size")
    ("gamma", po::value<double>()->default_value(0.32), "distance scale")
    ;

  po::positional_options_description pod;
  pod.add("ins", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
  po::notify(vm);

  if(vm.count("help") || !vm.count("ins")) {
    cout << desc << endl;
    return 0;
  }

  setupRandom(vm["seed"].as<unsigned>());
  string instance_name = vm["ins"].as<string>();

  target = vm["target"].as<int>();
  timelimit = vm["timelimit"].as<unsigned>();

  ifstream ins(instance_name);
  if (!ins.good()) {
    cout << "Can't open instance file " << instance_name << endl;
    return 1;
  }
  
  I.readInstance(ins);

  Solution S(I);

  start_time = chrono::system_clock::now();

  unsigned variant = vm["variant"].as<unsigned>();
  unsigned td = vm["td"].as<unsigned>();
  unsigned sm = vm["sm"].as<unsigned>();
  unsigned iv = vm["iv"].as<unsigned>();
  unsigned e = vm["e"].as<unsigned>();
  double gamma = vm["gamma"].as<double>();

  if(variant == 1)
    hhta1(S, start_time, td, sm, iv, e, gamma);
  else
    hhta2(S, start_time, sm, e, gamma);

  S.value += (I.btb * I.P);
  cout << S.value << endl;
}