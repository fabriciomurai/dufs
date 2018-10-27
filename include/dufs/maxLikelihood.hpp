#include <map>
#include <set>
//#include "RandomLib/Random.hpp"
#include <gsl/gsl_rng.h>
using namespace std;

#define _FLOAT long double
typedef pair<string, _FLOAT> LabelPi;

namespace dufs
{

#define RANDOM_VERTEX_SAMPLING 0
#define RANDOM_WALK_SAMPLING 1


class Counts {
    public:
        int rv_samples; // number of random vertex samples
        int rw_samples; // number of random walk samples
        map<string,_FLOAT> data_e;
        map<string,_FLOAT> data_v;
        map<string,_FLOAT> data_mu;
        set<string> labels;

        void push(string label, _FLOAT pi, _FLOAT counts, int type) {
            labels.insert(label);
            if (type == RANDOM_WALK_SAMPLING) {
                rw_samples += counts;
                data_e[label] += counts;
                data_mu[label] += counts/pi;
            } else if (type == RANDOM_VERTEX_SAMPLING) {
                rv_samples += counts;
                data_v[label] += counts;
            } else
                throw 1;
        }

        int size() { return labels.size(); }
	Counts(): rv_samples(0), rw_samples(0) {}
};

map<string,_FLOAT> maxLikelihoodFS(Counts c, gsl_rng *rnd = NULL);
map<string,_FLOAT> maxLikelihoodFS(map<string,_FLOAT> n, map<string,_FLOAT> m,
        map<string,_FLOAT> mu, gsl_rng *rnd = NULL);
} // namespace dufs
