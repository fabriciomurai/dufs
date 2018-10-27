/*
 * Copyright (c) 2018 Fabricio Murai (<email>)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */

#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dufs/maxLikelihood.hpp"
#include <gsl/gsl_multimin.h>
//#include <gsl/gsl_randist.h>
using namespace std;

#define RANDOM_STARTS 1
#define MAX(a,b) ( (a)>(b) ? (a) : (b) )

//const gsl_rng_type * rng_type;
//gsl_rng * r;

namespace dufs
{

void initialize_estimate_FS(gsl_vector *alpha, Counts *n) {
    
    double norm;
    set<string>::reverse_iterator rit = n->labels.rbegin();
    if(n->data_e[*rit] == 0){
        fprintf(stderr, "#WARNING: no random walk sample was observed for last element in n. Skipping initialization...\n");
        return;
    } else {
        norm = log(n->data_mu[*rit]);
    }

    // initialize estimate to p_k/p_{LEN}, i.e. alpha = log(p_k)-log(p_{LEN})
    set<string>::iterator it;
    int k;
    for(k=0, it=n->labels.begin(); it != n->labels.end(); k++, it++) {
        if( n->data_e[*it])
            gsl_vector_set(alpha, k, log(n->data_mu[*it]) - norm);
    }
}

void initialize_estimate_rwsamples(gsl_vector *alpha, Counts *n) {
    set<string>::iterator it;
    _FLOAT eps = numeric_limits<_FLOAT>::max();

    int k;
    for(k=0, it=n->labels.begin(); it != n->labels.end(); k++, it++) {
        _FLOAT mu = n->data_mu[*it];
        if( mu > 0.0 ) {
            gsl_vector_set(alpha, k, log(mu));
            eps = min(eps, mu);
        }
    }

    // Replace -Inf by log(eps) (smallest finite number)
    for(k=0, it=n->labels.begin(); it != n->labels.end(); k++, it++) {
        _FLOAT mu = n->data_mu[*it];
        if( mu == 0.0 )
            gsl_vector_set(alpha, k, log(eps));
    }

    // alpha^prime_i <- alpha_i - alpha_W
    gsl_vector_add_constant(alpha,-gsl_vector_get(alpha,k-1));
}

void initialize_estimate_random(gsl_vector *alpha, Counts *n, gsl_rng *rnd) {
    
    set<string>::iterator it;
    int k;
    for(k=0, it=n->labels.begin(); it != n->labels.end(); k++, it++)
        gsl_vector_set(alpha, k, log(-log(gsl_rng_uniform(rnd))));

    // alpha^prime_i <- alpha_i - alpha_W
    gsl_vector_add_constant(alpha,-gsl_vector_get(alpha,k-1));
}

double my_f(const gsl_vector *v, void *params) {
    Counts *n = static_cast<Counts *>(params);
    double sum = 0.0;
    double sum_exp = 0.0;
    double sum_kexp = 0.0;
    int k;

    set<string>::iterator it;
    for(k=0, it=n->labels.begin(); it != n->labels.end(); k++, it++) {
        double alpha_k = gsl_vector_get(v, k);
        sum +=(n->data_v[*it]+n->data_e[*it])*alpha_k;
        sum_exp += exp(alpha_k);
        if( n->data_e[*it])
            sum_kexp += exp(alpha_k)*n->data_e[*it]/n->data_mu[*it];
    }

    return n->rw_samples*log(sum_kexp)+n->rv_samples*log(sum_exp)-sum;
}

void my_df(const gsl_vector *v, void *params, gsl_vector *df) {
    Counts *n = static_cast<Counts *>(params);
    double sum_exp = 0.0;
    double sum_kexp = 0.0;
    int k;

    set<string>::iterator it;
    for(k=0, it=n->labels.begin(); it != n->labels.end(); k++, it++) {
        double alpha_k = gsl_vector_get(v, k);
        sum_exp += exp(alpha_k);
        if( n->data_e[*it])
            sum_kexp += exp(alpha_k)*n->data_e[*it]/n->data_mu[*it];
    }
    for(k=0, it=n->labels.begin(); it != n->labels.end(); k++, it++) {
        double exp_alpha_k = exp(gsl_vector_get(v, k));
        double grad;
        if( n->data_e[*it]) {
            grad = -n->rv_samples*exp_alpha_k/sum_exp + n->data_v[*it] -
                    n->rw_samples*exp_alpha_k*n->data_e[*it]/n->data_mu[*it]/sum_kexp + n->data_e[*it];
        } else {
            grad = -n->rv_samples*exp_alpha_k/sum_exp + n->data_v[*it];
        }
        gsl_vector_set(df, k, -grad );
    }
    gsl_vector_set(df,n->size()-1,0.0);
}

void my_fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df) {
    Counts *n = static_cast<Counts *>(params);
    double sum = 0.0;
    double sum_exp = 0.0;
    double sum_kexp = 0.0;
    int k;

    set<string>::iterator it;
    for(k=0, it=n->labels.begin(); it != n->labels.end(); k++, it++) {
        double alpha_k = gsl_vector_get(v, k);
        sum +=(n->data_v[*it]+n->data_e[*it])*alpha_k;
        sum_exp += exp(alpha_k);
        if( n->data_e[*it])
            sum_kexp += exp(alpha_k)*n->data_e[*it]/n->data_mu[*it];
    }
    for(k=0, it=n->labels.begin(); it != n->labels.end(); k++, it++) {
        double exp_alpha_k = exp(gsl_vector_get(v, k));

        double grad;
        if( n->data_e[*it]) {
            grad = -n->rv_samples*exp_alpha_k/sum_exp + n->data_v[*it] -
                    n->rw_samples*exp_alpha_k*n->data_e[*it]/n->data_mu[*it]/sum_kexp + n->data_e[*it];
        } else {
            grad = -n->rv_samples*exp_alpha_k/sum_exp + n->data_v[*it];
        }
        gsl_vector_set(df, k, -grad );
    }
    gsl_vector_set(df,n->size()-1,0.0);
    *f = n->rw_samples*log(sum_kexp)+n->rv_samples*log(sum_exp)-sum;
}

/*double my_f(const gsl_vector *v, void *params) {
    Counts *n = static_cast<Counts *>(params);
    int len=n->data.size();
    double avg_deg = 0.0;
    for(int k=0; k < len; k++) {
        avg_deg += n->data[k].degree*gsl_vector_get(v, k);
        assert(gsl_vector_get(v, k) > 0);
    }

    double g = -n->rw_samples*log(avg_deg)-n->rv_samples;
    for(int k=0; k < len; k++) {
        g += (n->data[k].rv_counts+n->data[k].rw_counts)*log(gsl_vector_get(v, k)) +
            n->data[k].rw_counts*log(n->data[k].degree);
    }
    g *= -1;
    cerr << "my_f:" << g << endl;

    return g;
}

void my_df(const gsl_vector *v, void *params, gsl_vector *df) {
    Counts *n = static_cast<Counts *>(params);
    int len=n->data.size();
    double avg_deg = 0.0;
    for(int k=0; k < len; k++) {
        avg_deg += n->data[k].degree*gsl_vector_get(v, k);
    }
    for(int k=0; k < len; k++) {
        double grad = -n->rw_samples*(n->data[k].degree/avg_deg)+
            (n->data[k].rv_counts+n->data[k].rw_counts)/gsl_vector_get(v, k)-
            n->rv_samples;
        gsl_vector_set(df, k, -grad );
    }
}

void my_fdf(const gsl_vector *v, void *params, double *f, gsl_vector *df) {
    Counts *n = static_cast<Counts *>(params);
    int len=n->data.size();
    double avg_deg = 0.0;
    for(int k=0; k < len; k++) {
        avg_deg += n->data[k].degree*gsl_vector_get(v, k);
    }

    *f = -n->rw_samples*log(avg_deg)-n->rv_samples;
    for(int k=0; k < len; k++) {
        *f += (n->data[k].rv_counts+n->data[k].rw_counts)*log(gsl_vector_get(v, k)) +
            n->data[k].rw_counts*log(n->data[k].degree);
        double grad = -n->rw_samples*(n->data[k].degree/avg_deg)+
            (n->data[k].rv_counts+n->data[k].rw_counts)/gsl_vector_get(v, k)-
            n->rv_samples;
        gsl_vector_set(df, k, -grad );
    }
    *f *= -1;

    cerr << "my_fdf:" << *f << endl;
}*/


/*map<int,_FLOAT> fixedPointFS(Counts c) {
    int dim = c.size();
    vector<_FLOAT> theta(dim);
    _FLOAT abs_tol = 1e-6;
    int max_it = 10000;

    _FLOAT sum=0.0;
    for(int k=0; k < dim; k++) {
        theta[k] = c.data[k].rw_counts*1.0/c.data[k].degree;
        sum += theta[k];
    }
    _FLOAT avg_deg = 0.0;
    for(int k=0; k < dim; k++) {
        theta[k] /= sum;
        avg_deg += c.data[k].degree*theta[k];
        //cout << "--> theta[" << k << "]=" << theta[k] << endl;
    }

    _FLOAT diff = 0;
    int it = 0;
    do {
        it++;
        diff = 0.0;
        sum = 0.0;
        for(int k=0; k < dim; k++) {
            _FLOAT old_theta = theta[k];
            theta[k] = (c.data[k].rv_counts+c.data[k].rw_counts)/
                (c.rv_samples+c.rw_samples*c.data[k].degree/avg_deg);
            sum += theta[k];
            diff += (theta[k]-old_theta)*(theta[k]-old_theta);
        }

        avg_deg = 0.0;
        for(int k=0; k < dim; k++) {
            theta[k] /= sum;
            avg_deg += c.data[k].degree*theta[k];
        }
    } while( sqrt(diff) > abs_tol && it < max_it );

    map<int,_FLOAT> est_pmf;
    for(int i=0; i < dim; i++) {
        est_pmf[c.data[i].degree] = theta[i];
    }
    return est_pmf;
}*/

map<string,_FLOAT> maxLikelihoodFS(Counts c, gsl_rng *rnd) {
    if (!rnd) {
        gsl_rng_env_setup();
        rnd = gsl_rng_alloc(gsl_rng_default);
    }
    int dim = c.size();
    gsl_vector *alpha = gsl_vector_calloc(dim);
    gsl_vector *df = gsl_vector_alloc(dim);

    size_t iter = 0;
    int status;
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s = NULL;

    gsl_multimin_function_fdf my_func;
    my_func.n = dim;
    my_func.f = my_f;
    my_func.df = my_df;
    my_func.fdf = my_fdf;
    my_func.params = &c;

    //T = gsl_multimin_fdfminimizer_conjugate_fr;
    T = gsl_multimin_fdfminimizer_vector_bfgs2;

    gsl_rng_env_setup();
    //rng_type = gsl_rng_mt19937;
    //r = gsl_rng_alloc (rng_type);

    for( int j=0; j < RANDOM_STARTS; j++) {
        //initialize_estimate_random(alpha, &c, rnd);
        initialize_estimate_rwsamples(alpha, &c);
        if(s) gsl_multimin_fdfminimizer_free(s);
        s = gsl_multimin_fdfminimizer_alloc(T, dim);

        gsl_multimin_fdfminimizer_set(s, &my_func, alpha, 0.01, 1e-3);
        iter = 0;
        do {
          iter++;
          status = gsl_multimin_fdfminimizer_iterate(s);
          if(status)
              break;

          //fprintf (stdout, "iteration %lu (f=%10.5e).\n", iter, s->f);
          status = gsl_multimin_test_gradient(s->gradient, 1e-3);
          //if(status == GSL_SUCCESS)
          //    break;

        } while(status == GSL_CONTINUE && iter < 10000);

        // teste envolvendo random starts
        //fprintf (stdout, "Stopped at iteration %lu (f=%10.5f).\n", iter, s->f);

        //fprintf(stdout, "#%03d",j);
        //for(int i=0; i < dim; i++)
        //    fprintf(stdout, " %e", gsl_vector_get(s->x,i));
        //fprintf(stdout, "\n");

    }

    //gsl_rng_free (r);

    _FLOAT sum_exp_alpha = 0.0;
    for(int i=0; i < dim; i++)
        sum_exp_alpha += exp(gsl_vector_get(s->x,i));

    map<string,_FLOAT> est_pmf;
    set<string>::iterator it;
    int i;
    for(i=0, it=c.labels.begin(); i < dim; i++, it++) {
        _FLOAT val = exp(gsl_vector_get(s->x, i))/sum_exp_alpha;
        est_pmf[*it] = val;
        //cout << "theta[" << c.data[i].degree << "]=" << val << endl;
    }

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(alpha);
    gsl_vector_free(df);

    return est_pmf;
}


map<string,_FLOAT> maxLikelihoodFS(map<string,_FLOAT> rv_obs, map<string,_FLOAT> rw_obs,
        map<string,_FLOAT> mu, gsl_rng *rnd) {
    if (!rnd) {
        gsl_rng_env_setup();
        rnd = gsl_rng_alloc(gsl_rng_default);
    }
    map<int,int> all_degrees;
    Counts n;
    n.data_v = rv_obs;
    n.data_e = rw_obs;
    n.data_mu = mu;

    for(map<string,_FLOAT>::iterator it=rv_obs.begin(); it != rv_obs.end(); it++) {
        n.labels.insert(it->first);
        n.rv_samples += it->second;
    }

    for(map<string,_FLOAT>::iterator it=rw_obs.begin(); it != rw_obs.end(); it++) {
        n.labels.insert(it->first);
        n.rw_samples += it->second;
    }


    //return fixedPointFS(n);
    return maxLikelihoodFS(n, rnd);
}

/*int test(void) {
    Counts c;
    int d, rv_counts, rw_counts;
    while( scanf("%d %d %d", &d, &rv_counts, &rw_counts) == 3 ) {
        //TODO: create the push function in c
        c.push(Triple (d,rv_counts,rw_counts));
    }
    map<int,_FLOAT> res = maxLikelihoodFS(c);
    for(map<int,_FLOAT>::iterator it=res.begin(); it !=res.end(); it++)
        cout << "theta[" << it->first << "]=" << it->second << endl;

    return 0;
}*/

} // namespace dufs
