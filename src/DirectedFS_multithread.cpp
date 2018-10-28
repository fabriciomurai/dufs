/*
 * Copyright (c) 2018 Fabricio Murai (<email>)
 * Distributed under the MIT License.
 * See accompanying file LICENSE.md or copy at http://opensource.org/licenses/MIT
 */
#include <iostream>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <inttypes.h>
#include <limits.h>
#include <cmath>
#include <stack>
#include <mutex>
#include <stdio.h>
#include "dufs/maxLikelihood.hpp"
#include "dufs/ThreadPool.h"
#include <execinfo.h>
#include <signal.h>
#include <getopt.h>
#include <cassert>
#include <numeric>


#define _FLOAT long double
#define _PTYPE int
#define _FLAGTYPE uint8_t
#define STDERR_FILENO 2

#define VERTEXSAMPLING 1
#define FRONTIERSAMPLING 2
#define EDGESAMPLING 3
#define RWSAMPLING 5
#define SNOWBALLSAMPLING 4


#define ADDED_EDGE_FLAG 2
#define UNDIRECTED_GRAPH_EDGE_FLAG 1
#define NO_UNDIRECTED_GRAPH_EDGE_FLAG 0
#define INVALID_EDGE_FLAG -1

#define PRINT_PMF
//#define PRINT_CCDF

//largest positive value an int32_t can hold.
//#define INT32_MAX 0x7fffffffL

//largest negative value an int32_t can hold.
//#define INT32_MIN (-INT32_MAX - 1L) 

#define GET_OUTDEGREE 0
#define GET_INDEGREE 1
#define GET_BOTH 2

// define errors
#define AVERAGE_UNDEFINED_WHEN_INEDGES_HIDDEN 10

void myexit(char *msg, int code=2) {
    cout << "ERROR:" << msg << endl;
    exit(code);
}

// typedefs
typedef map<string,_FLOAT> PDF;

// declare global variables
char *filename;
int32_t independence_cost;
int32_t marginal;
int nthreads;
int see_incoming_edges;

// initialize global flags
int avg_given = 0;
int not_recursive = 0;
int reduce_var = 0;
int use_hybrid = 0;
int use_durw = 0;
int use_mrw = 0;
int propto_degree = 0;
int inv_propto_degree = 0;
int COST_REPEATED_VISIT = 1;

//
//int32_t maxOutdegree = 0;
PDF mu_star;



void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}




using namespace std;
int njumps_limit;
//RandomLib::Random *rnd;        // r created with random seed
//mutex edge_iterator_mutex;

//#include <stdio>

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)<(b)?(b):(a))
#define ABSMAX(a,b) (abs(a)<abs(b)?(b):(a))

#define MAXGTOHIST (10)
#define NOPRINTORIG
//#define DEBUGFRONT2
//#define DEBUGDURW
//#define DEBUG
//#define TESTDIST


struct edge {
  int32_t u,v;
  _PTYPE p;
};

struct orderededge {
  edge *e;
  int32_t order;
};


struct vertex {
    int32_t  deg, sampledNums, undirected_deg;
    int32_t  indeg, outdeg;
    int32_t  neigsize;
    int32_t id;
    _FLAGTYPE flag;
    edge** edges;
    edge* virtualedge;
    _FLOAT dv;
};

void *zrealloc(void *ptr, size_t newsize, size_t oldsize) {
  int8_t *p = (int8_t *) realloc(ptr,newsize);
  if (p == NULL) {
      cerr << "Could not allocate "<<((_FLOAT)(newsize))/(1024*1024)<<"Kbytes. Aborting."<<endl;
      abort();
  }
  memset(p+oldsize,0,newsize - oldsize);
  return (void*)p;
}

void *zalloc(size_t newsize) {
  int8_t *p = (int8_t *) malloc(newsize);
  if (p == NULL) {
      cerr << "Could not allocate "<<((_FLOAT)(newsize))/(1024*1024)<<"Kbytes. Aborting."<<endl;
      abort();
  }
  memset(p,0,newsize);
  return (void*)p;
}

class Gdb {
protected:
  int32_t reservedspaceE;
  int32_t reservedspaceV;
  int32_t edgeindex;
  int32_t N;
  int32_t _MAXID;
  int32_t _MAXDEGREE;
  int32_t _MAXDEGREEID;
  int32_t noedges;
  
  int32_t i_idE, i_idV;
  gsl_rng *rnd;        // r created with random seed
  
  vertex *V;
  edge *E;

  void updatevertex(int32_t u) {
//    cerr << "add="<<pe<<endl;
    if (u >= reservedspaceV) {
      cerr << "Fatal Error: More vertices than previously allocated = "<<reservedspaceV<<endl;
      abort();
    }
    _MAXID = MAX(_MAXID,u);
    if (V[u].edges == NULL) {
      N++;
      V[u].neigsize = 4;
      V[u].edges = (edge **) zalloc(V[u].neigsize*sizeof(edge*));
    }
    if (V[u].deg >= V[u].neigsize) {
      V[u].edges = (edge **) zrealloc(V[u].edges, (V[u].neigsize*2)*sizeof(edge**), V[u].neigsize*sizeof(edge**));
      V[u].neigsize *= 2;
    }
    V[u].id = u;
  }


  
public:
  Gdb() {
  }


  void setRNG(gsl_rng *rng) {
      rnd = rng;
  }

  int myrand(int max) {
      return gsl_rng_uniform_int(rnd, max);
  }

  double myUnifrand() {
      // r.Fixed() is in the interval (0, 1)
      return gsl_rng_uniform(rnd);
  }

  void init(int32_t p_reservedspaceE, int32_t p_reservedspaceV) {
    edgeindex = 0;
    N = 0;
    _MAXDEGREE = 0;
    _MAXID = 0;
    noedges = 0;
    reservedspaceE = p_reservedspaceE;
    reservedspaceV = p_reservedspaceV;
    E = (edge *) zalloc(reservedspaceE * sizeof(edge));
    V = (vertex *) zalloc(reservedspaceV * sizeof(vertex));
  }
  
  inline _FLOAT samplingFreq(vertex *) {
      return 1.0/((_FLOAT)size());
  }
  
  inline vertex* getvertex(int32_t u) {
    return (exists(u))? &V[u] : NULL;
  }
  
  inline void edge_iterator_start() {
    i_idE = 0;
  }
  
  inline bool exists_e(int32_t id) {
      return ((id < edgeindex)? ((E[id].v > 0)? true : false ): false);
  }
  
  bool exists(edge *e) {
      return ((e != NULL)? ((e->v > 0)? true : false ): false);
  }
  
  edge* next_edge() {
      while ((!exists_e(i_idE)) && (i_idE < edgeindex)) {
          i_idE++;
      }
      return (exists_e(i_idE)? &E[i_idE++] : NULL );
  }
  
  inline void vertex_iterator_start() {
    i_idV = 0;
  }
  
  vertex* next_vertex() {
    while ((!exists(i_idV)) && (i_idV <= MAXID()))
      i_idV++;
    if (i_idV > MAXID())
      return NULL;
    return &V[i_idV++];
  }
  
  inline edge* getedge(int32_t id) {
      return (exists(&E[id])? &E[id] : NULL);
  }
  
/*  vertex* getneighbor(int32_t u, int32_t k) {
    if (exists(u)) {
      if (k >= V[u].deg)
        return NULL;
      return &V[(V[u].edges[k])->v];
    }
    return NULL;
  }
  */
  vertex* getneighbor(vertex *u, int32_t k) {
    if (u != NULL) {
        if (!exists(u->edges[k])) return NULL;
        return &V[(u->edges[k])->v];
    }
    return NULL;
  }
  
  edge* getedge(vertex *u, int32_t k) {
    if (u != NULL) {
        if (k >= u->deg) return NULL;
        if (!exists(u->edges[k])) return NULL;
        return u->edges[k];
    }
    return NULL;
  }
  
  vertex* getrandomneighbor(vertex *u) {
    int32_t i;
    if (u != NULL) {
      if (u->deg == 0) return u;
      i = myrand(u->deg);
      if (!exists(u->edges[i])) return NULL;
      return &V[(u->edges[i])->v];
    }
    return NULL;
  }
  
  edge* randomedge() {
      int32_t r;
      
      do {
          r = myrand(edgeindex);
      } while (!exists(&E[r]));
      return &E[r];
  }
  
  edge* randomedge(vertex *u) {
      edge *e;
      if (u != NULL) {
          if (u->deg > 0) {
              do {
                  e = u->edges[myrand(u->deg)];
              } while (!exists(e));
              return e;
          }
      }
      return NULL;
  }
    
  vertex *mostpopularvertex() {
      vertex *v;
      do {
          v = randomvertex();
      } while (v->deg < 10000);
      return (v);
  }
  
  vertex* randomvertex() {
    int u = (myrand(_MAXID+1));
    while (!exists(u)) 
      u = (myrand(_MAXID+1));
    return &V[u];
  }
  
  inline int32_t indegree(vertex *u) {
      return ((u != NULL)? (u->indeg) : (-1));
  }
  
  inline int32_t outdegree(vertex *u) {
      return ((u != NULL)? (u->outdeg) : (-1));
  }
  
  inline int32_t degree(vertex *u) {
      return ((u != NULL)? (u->deg) : (-1));
  }
  
  inline int32_t degree(int32_t u) {
      return degree(getvertex(u));
  }
  
  inline int32_t indegree(int32_t u) {
      return indegree(getvertex(u));
  }
  
  inline int32_t outdegree(int32_t u) {
      return outdegree(getvertex(u));
  }
  
  inline int32_t MAXID() {
    return _MAXID;
  }

/*  inline int32_t MAXgValue() {
    return _MAXgValue;
  }
*/

 inline int32_t MAXDEGREE() {
    return _MAXDEGREE;
 }

  inline int32_t size() {
    return N;
  }

  inline int32_t no_edges() {
    return noedges;
  }

  inline _FLAGTYPE flag(int32_t n) {
    return ( exists(n)? (V[n].flag) : (-1)  );
  }

  inline _FLAGTYPE flag(vertex *v) {
      return ( exists(v)? (v->flag) : (-1)  );
  }
  
    inline void setflag(int32_t n, _FLAGTYPE flag) {
        if (exists(n)) {
            V[n].flag = flag;
        }
    }

    inline void setflag(vertex *v, _FLAGTYPE flag) {
        if (exists(v)) {
            v->flag = flag;
        }
    }

   void zeroflags() {
    zeroflags(~0x00);
   }

 void zeroflags(int8_t mask) {
    int32_t i;
    for (i=0; i <= MAXID(); i++) {
        setflag(i,(flag(i) & (~mask)));
    }
 }

 void reset_undirected_graph() {
    int32_t i,j;
    vertex *u; 
    edge *e;
    
    for (i=0; i <= MAXID(); i++) {
        if (exists(i)) {
            u = getvertex(i);
            u->undirected_deg = -1;
            u->dv = 0;
            for (j = 0; j < u->deg; j++) {
                e = getedge(u,j);
                if (e == NULL) {
                    cerr << "EDGE e == NULL!!!, vertex "<<i<<" vertex number "<<j<<" and degree "<<u->deg<<endl;
                } else {
                    if (getedge(u,j)->p == ADDED_EDGE_FLAG)  
                        getedge(u,j)->p = INVALID_EDGE_FLAG;
                    if (getedge(u,j)->p == UNDIRECTED_GRAPH_EDGE_FLAG)  
                        getedge(u,j)->p = NO_UNDIRECTED_GRAPH_EDGE_FLAG;
                }
            }
        }
    }
 }

 void zeroedgep() {
  int32_t i;
  for (i=0; i < edgeindex; i++) {
    E[i].p = 0;
  }
}

  void setAlledgep(_FLOAT p) {
    int32_t i;
    for (i=0; i < edgeindex; i++) {
      E[i].p = p;
    }
  }

  inline bool exists(int32_t n) {
    return ((n <= MAXID())? (V[n].edges != NULL) : false);
  }
  
  inline bool exists(vertex *v) {
      return (v != NULL);
  }
  
  edge * reverseEdge(edge *e) {
      vertex *v;
      int32_t k;
      
      if (exists(e)) {
          v = getvertex(e->v);
          if (exists(v)) {
              for (k =0; k < v->deg; k++) {
                  if (v->edges[k]->v == e->u) 
                      return v->edges[k];
              }
          }
      }
      return NULL;
  }
  
  void trimGraphDifferent(_FLAGTYPE flag) {
      int32_t i,j;
      int32_t N_MAXDEGREE = 0;
      int32_t N_MAXID = 0;

      for (i=0; i <= MAXID(); i++) {
          if (exists(i)) {
              if (V[i].flag != flag) {
                  // edge is removed
                  N--;
                  if (V[i].edges != NULL) {
                      for (j = 0; j < V[i].deg; j++) {
                          if (exists(V[i].edges[j])) {
                              noedges--;
                              V[i].edges[j]->u = -1;
                              V[i].edges[j]->v = -1;
                          }
                      }
                      free(V[i].edges);
                      V[i].edges = NULL;
                  }
                  V[i].deg = 0;
              } else {
                  // edge stays
                  if (N_MAXDEGREE < V[i].deg) {
                      N_MAXDEGREE = V[i].deg;
                      _MAXDEGREEID = i;
                  }
                  N_MAXID = MAX(N_MAXID,V[i].id);
              }
          }
      }
      _MAXDEGREE = N_MAXDEGREE;
      _MAXID = N_MAXID;
  }
  

  void addedge(int32_t u, int32_t v, _PTYPE p, bool trueedge) {
    vertex *up;
    int32_t k;

    if (edgeindex >= reservedspaceE) {
      cerr << "Fatal Error: More edges than previously allocated = " << reservedspaceE << endl;
      abort();
    }
    // see if edge is repeated
    up = getvertex(u);
    if (up != NULL) {
      for (k=0; k < up->deg; k++) {
        if (up->edges[k] != NULL) {
          // cerr << "q " << up->edges[k] << "  deg = "<<up->deg<<endl;
          if (up->edges[k]->v == v) {
          // repeated edge... update and get out
                if (trueedge) {
                    if( up->edges[k]->p == NO_UNDIRECTED_GRAPH_EDGE_FLAG ) {
                        cerr << "#Trying to add repeated edge " << u << "," << v << endl;
                        exit(1);
                    }
                    V[u].outdeg++;
                    //maxOutdegree = maxOutdegree>V[u].outdeg? maxOutdegree:V[u].outdeg;
                }
                if (trueedge) V[v].indeg++;
                if (trueedge) (up->edges[k])->p = p;
                return;
          }
        }
      }
    }
    
    E[edgeindex].u = u;
    E[edgeindex].v = v;
    E[edgeindex].p = p;
    if (trueedge) updatevertex(u);
    if (trueedge) updatevertex(v);
    if (V[u].deg >= V[u].neigsize) {
        cerr << "Error: Dynamic edge allocation full. Cannot allocate more dynamic edges for vertex "<<u<<"; it already has "<<V[u].deg<<" edges (its maximum)"<<endl;
        exit(0);
    }
    V[u].edges[V[u].deg] = &E[edgeindex];
    V[u].deg++;
    _MAXDEGREE = MAX(_MAXDEGREE,V[u].deg);
    if (trueedge){
        V[u].outdeg++;
        //maxOutdegree = maxOutdegree>V[u].outdeg? maxOutdegree:V[u].outdeg;
    } 
    if (trueedge) V[v].indeg++;
    if (trueedge) V[u].flag = 0;
    if (trueedge) V[v].flag = 0;
    edgeindex++;
    noedges++;

    //cerr << "current " << u << " indeg:" << V[u].indeg << " outdeg:" << V[u].outdeg << endl;
  }
  
  void print_edge_probabilities() {
      int32_t i;

      for (i=0; i < edgeindex; i++) {
          cerr<< E[i].u <<"\t"<< E[i].v << "\t" << E[i].p << endl;
      }
  }
  
};

class MySampler {
    private:
        vector<pair<_FLOAT,int32_t> > cdf2node;
        Gdb *pG;
    public:
        MySampler(Gdb *pG) {
            vertex *v;
            int i;
            _FLOAT norm, sum;

            // get proportions
            norm = 0.0;
            i=0;
            pG->vertex_iterator_start();
            while ((v = pG->next_vertex()) != NULL) {
                _FLOAT recip = 1.0/v->deg;
                cdf2node.push_back(make_pair(recip,v->id));
                norm += recip;
                i++;
            }
            //cout << "Norm:" << norm << endl;
            // accumulate and normalize proportions
            sum = 0.0;
            for(auto & item: cdf2node) {
                sum += item.first;
                item.first = sum/norm;
                //cout << "cdf:" << item.first << ", node:" << item.second << endl;
            }

            this->pG = pG;
        }

        int32_t draw() {
            _FLOAT u = pG->myUnifrand();
            auto comp = []( _FLOAT b, pair<_FLOAT,int32_t> a){
                        return b < a.first;
                    };
            auto item = upper_bound(cdf2node.begin(), cdf2node.end(), u, comp);
            //cout << "u:" << u << ", item->second:" << item->second << endl;
            return item->second;
        }
};


class Graph : public Gdb {

  private:

      void read_file(const char *file)
    {
        FILE *f;
        int32_t i,j;
        int32_t maxid = 0;
        int32_t p_reservedspaceE;
        char st[1000];
        char *line = NULL;
        
        p_reservedspaceE = 0;
        sprintf(st,"gunzip -c %s",file);
        f = popen(st,"r");
        if (f == NULL) {
            cerr << "# Error opening file "<<file<<endl;
            exit(1);
        }
        int read;
        size_t len;
        while ((read = getline(&line, &len, f)) != -1) {
           if (line[0] != '#' && sscanf(line,"%d %d\n",&i,&j) == 2){
              if (i != j) {
                  maxid = MAX(MAX(maxid,i+1),j+1);
                  p_reservedspaceE+=2;
              }
          }
        }

        pclose(f);
        if (maxid == 0) {
            cerr << "# Error with graph file "<<file<<endl;
            exit(1);
        }
        //cerr << "#Edges = "<<p_reservedspaceE<<endl;
        //cerr << "#Vertices = "<<maxid+1<<endl;
        init(p_reservedspaceE,maxid+1);
        f = popen(st,"r");
        while ((read = getline(&line, &len, f)) != -1) {
          if (line[0] != '#' && sscanf(line,"%d %d\n",&i,&j) == 2){
              if (i != j) {
                //cerr << "edge ("<<i+1<<","<<j+1<<")"<<endl;
                addedge(i+1,j+1,NO_UNDIRECTED_GRAPH_EDGE_FLAG,true);
                // creates edge j->i (because the graph is directed! and we need an undirected graph)
                addedge(j+1,i+1,INVALID_EDGE_FLAG,false);
              }
          }
        }
        //cerr <<"#Maximum degree = "<<MAXDEGREE()<<endl;
        pclose(f);
    }
 
  public: 
    
    Graph(const char *filename) 
    :Gdb()
    { 
      read_file(filename);
    }

    // constructor for deep copy
    Graph(Graph *src) : Gdb() {
        //reservedspaceE = src.reservedspaceE;
        //reservedspaceV = src.reservedspaceV;
        edgeindex = src->edgeindex;
        N = src->N;
        _MAXID = src->_MAXID;
        _MAXDEGREE = src->_MAXDEGREE;
        _MAXDEGREEID = src->_MAXDEGREEID;
        noedges = src->noedges;
        i_idE = src->i_idE;
        i_idV = src->i_idV;

        init(src->reservedspaceE,src->reservedspaceV);
        src->edge_iterator_start();

        edge *e;
        while ((e = src->next_edge()) != NULL) {
            //fprintf( stdout, "Adding SRC %d %d %d\n", e->u, e->v, e->p );
            if(e->p!= INVALID_EDGE_FLAG) {
                addedge(e->u,e->v,e->p,true);
                addedge(e->v,e->u,INVALID_EDGE_FLAG,false);
            }
        }
        //this->edge_iterator_start();
        //while ((e = this->next_edge()) != NULL)
        //    fprintf( stdout, "DST %d %d %d\n", e->u, e->v, e->p );
              
    }
    
    ~Graph() {
    }
    
};

int32_t compare_ints (const void *a, const void *b) {
  if (*((int32_t*)a) > *((int32_t*)b))
    return 1;
  else if (*((int32_t*)a) < *((int32_t*)b))
    return -1;
  else
    return 0;
}    

class NodeSet {
  private:
    Graph *G;
    vertex **vertices;
    int32_t ndidx;
    int32_t reserved;
    int32_t totaldeg;
    int32_t it;
    double jump_weight;
    
    
  public:
    NodeSet(Graph *pG, double jump_weight=0.0) {
      G = pG;
      reserved = 10000;
      vertices = (vertex **) zalloc(sizeof(vertex *)*reserved);
      ndidx = 0;
      totaldeg = 0;
      this->jump_weight = jump_weight;
    }
    
    ~NodeSet() {
      free(vertices);
    }
    
    void clear(){
      totaldeg = 0;
      ndidx = 0;
    }
    
    inline void iterator_start() {
        it = 0;
    }
    
    inline vertex *next_vertex() {
        return element(it++);
    }
    
    void add(vertex *u) {
      if (ndidx >= reserved) {
        vertices = (vertex **) zrealloc(vertices,sizeof(vertex *)*(reserved+10000),sizeof(vertex *)*reserved);
        reserved += 10000;
      }
      vertices[ndidx++] = u;
      if( see_incoming_edges )
          totaldeg += G->degree(u) + jump_weight;
      else
          totaldeg += u->undirected_deg + jump_weight;
    }
    
    void eraseall(vertex *u) {
        int32_t i,k;
        
        for (i=0; i < ndidx; i++) {
            if (vertices[i] == u) {
                ndidx--;
                if( see_incoming_edges )
                    totaldeg -= G->degree(u) + jump_weight;
                else
                    totaldeg -= u->undirected_deg + jump_weight;
                for (k=i; k < ndidx; k++) {
                    vertices[k] = vertices[k+1];
                }
                vertices[ndidx] = 0;
            }
        }
    }
    
    void eraseone(vertex *u) {
      int32_t i,k;
      
      for (i=0; i < ndidx; i++) {
        if (vertices[i] == u) {
          ndidx--;
          if( see_incoming_edges )
              totaldeg -= G->degree(u) + jump_weight;
          else
              totaldeg -= u->undirected_deg + jump_weight;
          for (k=i; k < ndidx; k++) {
            vertices[k] = vertices[k+1];
          }
          vertices[ndidx] = 0;
          break; // should remove just one
        }
      }
    }
    
    vertex* element(int32_t k) {
      if (k >= size()) {
        return NULL;
      }
      return vertices[k];
    }
    
    vertex* randomvertex_fromedge() {
      int32_t e,i;
      
      if(size() == 1)
          return vertices[0];

      if(totaldeg > 0) {
          if( see_incoming_edges ) {
              #ifdef DEBUGFRONT3
              cerr << "Current totaldeg: " << totaldeg << endl; 
              #endif
              e = (G->myrand(totaldeg)) + 1;
              #ifdef DEBUGFRONT3
              cerr << "Chosen edge "<< e << " from " << totaldeg << " edges" << endl;
              #endif
              i = 0;
              while (e > 0) {
              #ifdef DEBUGFRONT3
              fprintf(stderr, "i = %d, vertices.size = %d, e = %d, G->degree(v) = %d\n", i, ndidx, e, G->degree(vertices[i])+jump_weight);
              #endif
                e -= G->degree(vertices[i++]) + jump_weight;
              }
             // cerr << "Chosen vertex "<< vertices[i-1]->id << " with degree " <<  vertices[i-1]->deg << endl;
          } else {
              e = (G->myrand(totaldeg)) + 1;
              i = 0;
              while (e > 0) {
              #ifdef DEBUGFRONT3
              fprintf(stderr, "totaldeg = %d, i = %d, vertices.size = %d, e = %d, id = %d, deg = %d\n",
                      totaldeg, i, ndidx, e, vertices[i]->id,
                      vertices[i]->undirected_deg+jump_weight);
              #endif
                e -= (vertices[i++])->undirected_deg + jump_weight;
              }
          }
      } else {
          i = G->myrand(size())+1;
      }
      return vertices[i-1];
    }
    
/*    void freq(NodeStatistics *freq) {
      int32_t i;
      
      freq->clear();
      for (i=0; i < idx; i++) {
        freq->inc(G->g(vertices[i]));
      }
    }
    */
/*    void print(const string str) {
      int32_t i;
      
      for (i=0; i < size(); i++) {
        cerr << str << "\t" << G->g(vertices[i]) << endl;
      }
    }
    */
    inline int32_t size() {
        return ndidx;
    }

    inline int32_t no_edges() {
        return totaldeg;
    }
};

class EdgeSet {
  private:
    Graph *G;
    int32_t reservedspace;
    int32_t edgeindex;
    int32_t it;
    orderededge *edges;
    
    inline int32_t space_left() {
        return (reservedspace - edgeindex);
    }
    
    public:
    
    EdgeSet(Graph *pG) {
      G = pG;
      // Vector with all vertex ids that were uniformly sampled
      reservedspace = 10000;
      edges = (orderededge *) zalloc(sizeof(orderededge)*reservedspace);
      edgeindex = 0;
    }
    
    ~EdgeSet() {
      free(edges);
    }
    
    void clear(){
      edgeindex = 0;
    }
    
    /*void print(char *str) {
      int32_t i;
      
      for (i = 0; i < edgeindex; i++) {
        cerr << str << "\t" << G->g(G->getvertex(edges[i].e->u)) << "\t" << G->g(G->getvertex(edges[i].e->v)) << "\t" << edges[i].order << endl;
      }
    }*/
    
    void iterator_start() {
        it = 0;
    }
    
    edge *next_edge() {
        if (it != edges[it].order) {
            cerr << "Wrong edge order, it = "<< it << " and edge order = "<<edges[it].order<<endl;
        }
        return (edges[it++].e);
    }
    
    
    void add(edge *e, int32_t order) {
      if (space_left() <= 0) {
        edges = (orderededge *) zrealloc(edges,sizeof(orderededge)*(reservedspace+1000),sizeof(orderededge)*reservedspace);
        reservedspace += 1000;
      }
      edges[edgeindex].e = e;
      edges[edgeindex].order = order;
      edgeindex++;
    }
    
    orderededge* element(int32_t k) {
      if (k < size()) {
        return &edges[k];
      }
      return NULL;
    }

    inline int32_t size() {
      return edgeindex;
    }
    
};

class FSResults {
private:
    map<LabelPi,_FLOAT> val_e;
    map<string,_FLOAT> val_v;
    Graph *Gr;
    int32_t sizeofEprime, sizeofVprime;
    inline Graph *G() {
        return Gr;
    }

    map<int,int> count_vsamples_without_esamples(PDF m) {
        // example input:
        // val_v = {1:1, 2:1,      4:2, 5:3, 6:3}
        // val_e = {   , 2:1, 3:2,              , 7:1}
        //
        // optional intermediate step (vertex labels without edge samples)
        //         {1:1,           4:2, 5:3, 6:3}
        //
        // output:
        // counts = {1:1, 2:1, 3:2} // '1'-> 1, '4' -> 2,  '5', '6' -> 3
        map<int,int> counts;
        for(auto const v_it: val_v)
            if((int)m[v_it.first] == 0)
                counts[(int)v_it.second] += 1;

        return counts;
    }

public:
    PDF stats;

    // constructor
    FSResults(Graph *pG) : Gr(pG), sizeofEprime(0), sizeofVprime(0) {}

    void add(vertex *v, _FLOAT pi) {
        if (pi > 0.0) {
            stringstream s;
            if( marginal == GET_OUTDEGREE )
                s << v->outdeg;
            else if( marginal == GET_INDEGREE )
                s << v->indeg;
            else if( marginal == GET_BOTH )
                s << v->indeg << '\t' << v->outdeg;

            val_v[s.str()] += 1.0;
            sizeofVprime++;
        }
    }

    void add(edge *e, _FLOAT pi) {
        if (pi > 0.0) {
            stringstream s;
            if( marginal == GET_OUTDEGREE )
                s << G()->outdegree(e->v);
            else if( marginal == GET_INDEGREE )
                s << G()->indegree(e->v);
            else if( marginal == GET_BOTH )
                s << G()->indegree(e->v) << '\t' << G()->outdegree(e->v);

            auto tmp = make_pair(s.str(),pi);
            val_e[tmp] += 1.0;
            sizeofEprime++;
        }

    }

    pair<PDF,PDF> consolidate(map<LabelPi,_FLOAT> fvalue, gsl_rng *rnd) {
        PDF est_pdf;
        //map<LabelPi,_FLOAT> fvalue = trueSum(); // original distribution
        //_FLOAT TrueSumg = normalize_TrueSum();  // normalization constant original distribution

        _FLOAT TrueSumg = 0.0;
        for(auto it=fvalue.begin(); it != fvalue.end(); it++)
            TrueSumg += it->second;

        PDF mu;  // estimated distribution based on random-walk samples
        PDF m;   // number of samples associated with given label
        for(map<LabelPi,_FLOAT>::iterator it = val_e.begin(); it != val_e.end(); it++) {
            m[it->first.first] += it->second;
            mu[it->first.first] += it->second/it->first.second; // mu[i] = \sum_j m_{i,j}/j
            //fprintf(stdout, "MU %s\t%Lf\t%Lf\n", it->first.first.c_str(), it->second, it->first.second);
        }

        // compute how many labels have val_v[label] >= 1, 2 and val_e[label] == 0
        if(!use_durw) {
            map<int,int> counts = count_vsamples_without_esamples(m);
            assert(counts[0] == 0);
            stats[string("e=0,v>0")] =
                std::accumulate(std::begin(counts), std::end(counts), 0,
                        [](const int previous, const std::pair<int,int>& p) { return previous+p.second; });
            stats[string("e=0,v>1")] = (int)(stats[string("e=0,v>0")]) - counts[1];
        }

        // original estimator (IMC 2010)
        if(!use_hybrid || use_durw) {
            return make_pair(mu,stats);
        }

        // hybrid estimator (TKDD 2016)
       cout << "Hybrid estimator";

        //for(auto it=mu_star.begin(); it != mu_star.end(); it++)
        //if(mu.find(it->first) != mu.end())
        //    fprintf(stdout, "%s\t%Lf\t%Lf\n", it->first.c_str(), mu[it->first], it->second);

        if(reduce_var) {
            cout << " +reduce_var";
              for( map<LabelPi,_FLOAT>::iterator it = fvalue.begin(); it != fvalue.end(); it++ ) {
                  string label = it->first.first;
                  if(mu.find(label) == mu.end()) {
                      sizeofVprime -= val_v[label];
                      val_v.erase(label);
                  }
              }
        }

        // compute estimator assuming average is known
        if(avg_given || not_recursive) {

          _FLOAT mu0 = 0.0;
          if(avg_given) {
              cout << " +avg_given";
              for( map<LabelPi,_FLOAT>::iterator it = fvalue.begin(); it != fvalue.end(); it++ ) {
                  string label = it->first.first;
                  if(mu.find(label) != mu.end()) {
                      //fprintf(stdout, "%s: %Lf * %Lf / %Lf\n", label.c_str(), it->second/TrueSumg, m[label], mu[label]);
                      mu0 += it->second/TrueSumg * m[label] / mu[label];
                  }
              }
          }
          if(not_recursive) {
              cout << " +not_recursive";
              mu0 = sizeofEprime / std::accumulate(
                      std::begin(mu), std::end(mu), (_FLOAT) 0.0,
                      [](const _FLOAT previous, const std::pair<string,_FLOAT>& p) {
                        return previous+p.second;
                      }
                    );

          }

          // note: OK if same label is seen more than once, same result will be obtained all times
          for( map<LabelPi,_FLOAT>::iterator it = fvalue.begin(); it != fvalue.end(); it++ ) {
              string label = it->first.first;
              if(mu.find(label) == mu.end()) {
                  est_pdf[label] = val_v[label]/sizeofVprime;
              } else {
                  est_pdf[label] = (val_v[label]+m[label])/
                      (sizeofVprime + sizeofEprime/mu0 *m[label]/mu[label] );
              }
              //fprintf(stdout, "%s: (%Lf+%Lf)/(%d + %d/%Lf * %Le/%Le)\n",
              //        label.c_str(), val_v[label], m[label], sizeofVprime, sizeofEprime, mu0, m[label], mu[label]);
          }

          // Normalization.
          // Step 1: compute constant
          long double norm = 0.0;
          for( map<string,_FLOAT>::iterator it = est_pdf.begin(); it != est_pdf.end(); it++ ){
            //cout << "A " << it->first << "\t" << it->second << endl;
            norm += it->second;
          }
          #ifdef DEBUG
          cout << "Norm:" << norm << endl;
          #endif

          // Step 2: divide each point by constant
          for( auto && pdf: est_pdf )
              pdf.second /= norm;
              //pdf.second /= 2;
        }
        else {
          // call gradient descent method to compute estimates
          est_pdf = dufs::maxLikelihoodFS(val_v, m, mu, rnd);

          //#ifdef DEBUG
          //_FLOAT mup = 0.0;
          //for( map<string,_FLOAT>::iterator it = est_pdf.begin(); it != est_pdf.end(); it++ ) {
          //    string label = it->first;
          //    if(mu.find(label) != mu.end()) {
          //        //fprintf(stdout, "%s: %Lf * %Lf / %Lf\n", label.c_str(), it->second/TrueSumg, m[label], mu[label]);
          //        mup += it->second * m[label] / mu[label];
          //    }
          //}
          //fprintf(stdout, "mu0: %Lf, mup: %Lf\n", mu0, mup);
          //#endif
        }
        cout << endl;
        return make_pair(est_pdf, stats);
    }

};

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
class F {
    private:
        Graph *Gr;
        
    public:
        F(Graph *pG) {
            Gr = pG;
        }
        
        virtual ~F() {
        }
        
        inline Graph *G() {
            return Gr;
        }
        
        virtual void reset() {
        }
        
        virtual void prepare_nextrun() {
        }
        
        virtual void add(edge *, int32_t , int8_t, _FLOAT ) {
        }

        virtual void add(vertex *, int32_t , int8_t, _FLOAT ) {
        }
        
        virtual void add_all(edge *, int8_t , _FLOAT , _FLOAT) {
        }
            
        virtual void  consolidate_TrueSum() {
        }
        

        virtual _FLOAT normalize_val(map<LabelPi,_FLOAT> &) {
            cerr << "F::normalize_val()" << endl;
            return 0.0;
        }
        
        virtual _FLOAT normalize_TrueSum() {
            cerr << "F::normalize_TrueSum()" << endl;
            return 0.0;
        }
        
        //virtual  _FLOAT* trueSum() {
        //    return NULL;
        //}
        virtual  map<LabelPi,_FLOAT> trueSum() {
            map<LabelPi,_FLOAT> tmp;
            return tmp;
        }

        virtual  map<LabelPi,_FLOAT> reallyTrueSum() { // workaround
            map<LabelPi,_FLOAT> tmp;
            return tmp;
        }
        
        virtual void consolidate() {
        }

        virtual void consolidate(PDF, PDF) {
        }
        
        virtual void print(const char *) {
        }
        virtual void printrealvalue(const char *preamble) {
        }
};

void call_from_thread() {
    std::cout << "Hello, World" << std::endl;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
class FrontierSampling {
    private:
        Graph *G;
        
    public:
        FrontierSampling(Graph *pG) {
            G = pG;
        }



        _FLOAT runsampled(int64_t samplebudget, int32_t front_size, int32_t initial_nodes, int64_t semiruns, bool VarDist, double jump_weight, F *Fsum, int64_t *){
            map<LabelPi,_FLOAT> trueSum = Fsum->reallyTrueSum();

            
            #ifdef TESTDIST
            int64_t t;
            int32_t degdist[100];
            memset(degdist,0,100*sizeof(int32_t));
            #endif

            //G->zeroedgep();
            if (Fsum != NULL) Fsum->reset();

            if (samplebudget % front_size != 0) {
                //cerr << "Warning in multiple random walks: Sample budget must be a multiple of the random walkers";
            }
            if (front_size % initial_nodes != 0) {
                //cerr << "Warning in multiple walkers starting on the same node: Number of walkers must be a multiple of the number of initial nodes" << endl;
            }

            //G->zeroflags();
            //G->reset_undirected_graph();
            //if (Fsum != NULL) Fsum->prepare_nextrun();

            ThreadPool tpool(nthreads);
            map<thread::id, Graph *> gpool;

            std::vector< std::future<int> > return_vals;
            for (thread::id &id: tpool.get_ids()) {
                return_vals.emplace_back(tpool.enqueue([&gpool,id]() {
                    gpool[id] = new Graph(filename);
                    cout << "Done " << id << " by " << this_thread::get_id() << endl;
                    return 0;
                }));
            }
            for(auto && return_val: return_vals)
                return_val.get();

            time_t start_time, current_time;
            time(&start_time);
            std::vector< std::future< pair<PDF,PDF> > > results;
            for(int run = 0; run < semiruns; ++run) {
                results.emplace_back(
                   tpool.enqueue([&,run] {

                       // get graph from pool
                       Graph *pG = gpool[this_thread::get_id()];
                       pG->zeroflags();
                       pG->reset_undirected_graph();

                       vertex *r, *v;
                       edge *e, *rev_e, *sampled_e;
                       int64_t sampled, total_steps;
                       int njumps = 0;
                       int32_t i,k;
                       NodeSet *samplednodes = new NodeSet(pG, jump_weight);
                       FSResults fsr(pG);
                       fsr.stats[string("budget")] = samplebudget;
                       fsr.stats[string("revisits")] = 0.0;
                       fsr.stats[string("walkers")] = initial_nodes;
                       fsr.stats[string("bpw")] = (samplebudget - independence_cost*initial_nodes)/(1.0*initial_nodes);

                       gsl_rng_env_setup();
                       gsl_rng *rnd = gsl_rng_alloc(gsl_rng_default);
                       gsl_rng_set(rnd, run);
                       pG->setRNG(rnd);

                       MySampler *mysampler = NULL;
                       if (inv_propto_degree)
                           mysampler = new MySampler(pG);

                       for (k =0; k < initial_nodes; k++) {
                           if (propto_degree) {
                               e = pG->randomedge();
                               r = pG->getvertex(e->v);
                               cout << "r->id: " << r->id << endl;
                           } else if (inv_propto_degree) {
                               r = pG->getvertex(mysampler->draw());
                           } else {
                               r = pG->randomvertex();
                           }
                           #ifdef DEBUGFRONT2
                           cout << "initial sample: " << r->id << endl;
                           #endif
                           r->undirected_deg = 0;
                           r->flag = 1;
                           for (i =0; i < pG->degree(r); i++) {
                               e = pG->getedge(r,i);
                               if ( pG->flag(e->v) == 0 ) {
                                   bool check_rev = false;
                                   if( e->p == NO_UNDIRECTED_GRAPH_EDGE_FLAG ) {
                                       e->p = UNDIRECTED_GRAPH_EDGE_FLAG;
                                       check_rev = true;
                                   } else if( see_incoming_edges && e->p == INVALID_EDGE_FLAG ) {
                                       e->p = ADDED_EDGE_FLAG;
                                       check_rev = true;
                                   }
                                   if( check_rev ) {
                                       rev_e = pG->reverseEdge(e);
                                       if (rev_e == NULL) {
                                           cerr << " reverse edge does not exist (" << e->v << "," << e->u << ")!!"<<endl;
                                           exit(0);
                                       }
                                       if( rev_e->p == NO_UNDIRECTED_GRAPH_EDGE_FLAG ) {
                                           rev_e->p = UNDIRECTED_GRAPH_EDGE_FLAG;
                                       } else if( rev_e->p == INVALID_EDGE_FLAG ) {
                                           rev_e->p = ADDED_EDGE_FLAG;
                                       }
                                   }
                               }
                               if ((e->p == UNDIRECTED_GRAPH_EDGE_FLAG) || (e->p == ADDED_EDGE_FLAG)) {
                                   #ifdef DEBUGFRONT2
                                   cout << "\tis connected to " << e->v << endl;
                                   #endif
                                   r->undirected_deg++;
                               }
                           }
                           for( int j = 0; j < front_size/initial_nodes; j++) {
                               samplednodes->add(r);
                           }
                           fsr.add(r,1.0);

                           #ifdef TESTDIST
                           //        degdist[r->deg]++;
                           #endif
                         }
                         if(mysampler)
                         delete mysampler;

                           // std::cout << "hello " << run << std::endl;
                           // std::this_thread::sleep_for(std::chrono::seconds(20));
                           // std::cout << "world " << run << std::endl;

                           sampled_e = NULL;
                           sampled = 0;
                           total_steps = 0;
                           edge virtualedge;
                           vertex *jumpvertex;
                           while (sampled < (samplebudget - independence_cost*initial_nodes) && njumps < njumps_limit) {
                               // randomly sample walker in proportion to degree
                               r = samplednodes->randomvertex_fromedge();
                               samplednodes->eraseone(r);

                               // compute undirected degree
                               r->undirected_deg = 0;
                               for (i =0; i < pG->degree(r); i++) {
                                   e = pG->getedge(r,i);
                                   if ((e->p == UNDIRECTED_GRAPH_EDGE_FLAG) || (e->p == ADDED_EDGE_FLAG))
                                       r->undirected_deg++;
                               }

                               // random jump or regular walk
                               if (pG->myUnifrand() <= (jump_weight/(jump_weight + r->undirected_deg))) {
                                   jumpvertex = pG->randomvertex();
                                   //cout << " random jump! : "<<jumpvertex->id;
                                   virtualedge.u = r->id;
                                   virtualedge.v = jumpvertex->id;
                                   virtualedge.p = INVALID_EDGE_FLAG;
                                   sampled_e = &virtualedge;
                                   if(jumpvertex->flag !=1)
                                   {
                                       sampled+=independence_cost;
                                   } else {
                                       sampled += COST_REPEATED_VISIT*independence_cost;
                                   }
                                   njumps++;
                                   fsr.stats[string("pjump")] += 1.0/(samplebudget- independence_cost*initial_nodes);
                               } else {
                                   do {
                                       e = pG->randomedge(r);
                                   } while ((e->p == NO_UNDIRECTED_GRAPH_EDGE_FLAG) || (e->p == INVALID_EDGE_FLAG));
                                   sampled_e = e;
                                   //cerr << " selected edge (" << e->r << "," << e->v << ")"<<endl;
                                   if(pG->flag(e->v) != 1)
                                   {
                                       sampled++;
                                   } else {
                                       sampled += COST_REPEATED_VISIT;
                                       fsr.stats[string("revisits")] += 1.0;
                                   }
                               }

                               v = pG->getvertex(sampled_e->v);
                               #ifdef DEBUGFRONT2
                               cout << "Selected edge: " << sampled_e->u << " -> " << sampled_e->v << endl; 
                               #endif

                               // update v' edges flags, compute undirected degree and store sample
                               v->undirected_deg = 0;
                               for (i =0; i < pG->degree(v); i++) {
                                   e = pG->getedge(v,i);
                                   // cerr << " bu";
                                   if ( pG->flag(e->v) == 0 ) {
                                       bool check_rev = false;
                                       if( e->p == NO_UNDIRECTED_GRAPH_EDGE_FLAG ) {
                                           e->p = UNDIRECTED_GRAPH_EDGE_FLAG;
                                           check_rev = true;
                                       } else if( see_incoming_edges && e->p == INVALID_EDGE_FLAG ) {
                                           e->p = ADDED_EDGE_FLAG;
                                           check_rev = true;
                                       }
                                       if( check_rev ) {
                                           rev_e = pG->reverseEdge(e);
                                           if (rev_e == NULL) {
                                               cerr << " reverse edge does not exist (" << e->v << "," << e->u << ")!!"<<endl;
                                               exit(0);
                                           }
                                           if( rev_e->p == NO_UNDIRECTED_GRAPH_EDGE_FLAG ) {
                                               rev_e->p = UNDIRECTED_GRAPH_EDGE_FLAG;
                                           } else if( rev_e->p == INVALID_EDGE_FLAG ) {
                                               rev_e->p = ADDED_EDGE_FLAG;
                                           }
                                       }
                                   }
                                   if ((e->p == UNDIRECTED_GRAPH_EDGE_FLAG) || (e->p == ADDED_EDGE_FLAG)) {
                                       //#ifdef DEBUGFRONT2
                                       //    cout << "\tis connected to " << e->v << endl;
                                       //#endif
                                       v->undirected_deg++;
                                   }
                               }
                               #ifdef DEBUGFRONT2
                               cout << "Adding " << v->id << ", weight" << v->undirected_deg+jump_weight << endl;
                               #endif
                               v->flag = 1;
                               fsr.add(sampled_e,(_FLOAT)v->undirected_deg+jump_weight);

                               //if (v->undirected_deg != debug_deg) cerr << "# deg("<<v->id<<")="<<v->undirected_deg<<" , old = "<<debug_deg << endl;

                               // Now select next vertex
                               //do {
                               //    e = pG->randomedge(v);
                               //} while ((e->p == NO_UNDIRECTED_GRAPH_EDGE_FLAG) || (e->p == INVALID_EDGE_FLAG));
                               
                               samplednodes->add(v);

                           #ifdef DEBUGFRONT
                           cerr << "Edge sampled: "<< u->id << "," << v->id<< " sampled so far "<<sampled<<endl;
                           #endif
				total_steps += 1;
				if( total_steps >= 2*samplebudget )
                                    break;
                           }
                           fsr.stats[string("sampled")] = sampled;
                           fsr.stats[string("revisits")] /= samplebudget;


                       delete samplednodes;
                       //delete pG;
                       return fsr.consolidate(trueSum, rnd);
                       //return fsr.est_pdf;
                   })
                );
            }
            time(&current_time);
            for(auto && result: results) {
                auto ret = result.get();
                Fsum->consolidate(ret.first, ret.second);
            };
            cout << "Time spent: " << current_time-start_time << endl;


            #ifdef TESTDIST
            cerr << "Dist = " << endl;
            t =0;
            for (k = 0; k < 14; k++) {
                t+=degdist[k];
            }
            for (k = 0; k < 14; k++) {
                cerr << "   d[" << k << "] = " << (_FLOAT)degdist[k]/t << endl;
            }
            #endif

            // has to return something
            return -1;
        }
};

class RandomWalkSampling {
 private:
  Graph *G;
  
 public:
  RandomWalkSampling(Graph *pG) {
    G = pG;
  }
  
  _FLOAT runsampled(int64_t samplebudget, int32_t initial_nodes, int64_t semiruns, bool VarDist, F *Fsum, int64_t *sedges) {

            map<LabelPi,_FLOAT> trueSum = Fsum->reallyTrueSum();

            
            #ifdef TESTDIST
            int64_t t;
            int32_t degdist[100];
            memset(degdist,0,100*sizeof(int32_t));
            #endif

            //G->zeroedgep();
            if (Fsum != NULL) Fsum->reset();

            if (samplebudget % initial_nodes != 0) {
                //cerr << "Warning in multiple random walks: Sample budget must be a multiple of the random walkers";
            }

            //G->zeroflags();
            //G->reset_undirected_graph();
            //if (Fsum != NULL) Fsum->prepare_nextrun();

            ThreadPool tpool(nthreads);
            map<thread::id, Graph *> gpool;

            std::vector< std::future<int> > return_vals;
            for (thread::id &id: tpool.get_ids()) {
                return_vals.emplace_back(tpool.enqueue([&gpool,id]() {
                    gpool[id] = new Graph(filename);
                    cout << "Done " << id << " by " << this_thread::get_id() << endl;
                    return 0;
                }));
            }
            for(auto && return_val: return_vals)
                return_val.get();

            time_t start_time, current_time;
            time(&start_time);
            std::vector< std::future< pair<PDF,PDF> > > results;
            for(int run = 0; run < semiruns; ++run) {
                results.emplace_back(
                   tpool.enqueue([&,run] {

                       // get graph from pool
                       Graph *pG = gpool[this_thread::get_id()];
                       pG->zeroflags();
                       pG->reset_undirected_graph();

                       vertex *r;
                       edge *e, *sampled_e;
                       int64_t sampled, total_steps;
                       int njumps = 0;
                       int32_t k;
                       FSResults fsr(pG);
                       fsr.stats[string("budget")] = samplebudget;
                       fsr.stats[string("revisits")] = 0.0;
                       fsr.stats[string("walkers")] = initial_nodes;
                       fsr.stats[string("bpw")] = (samplebudget - independence_cost*initial_nodes)/(1.0*initial_nodes);

                       gsl_rng_env_setup();
                       gsl_rng *rnd = gsl_rng_alloc(gsl_rng_default);
                       gsl_rng_set(rnd, run);
                       pG->setRNG(rnd);

                       vector<vertex *> samplednodes;
                       for (k =0; k < initial_nodes; k++) {
                           r = pG->randomvertex();
                           while( pG->degree(r) == 0 )
                               r = pG->randomvertex();
                           samplednodes.push_back(r);
                       }

                           // std::cout << "hello " << run << std::endl;
                           // std::this_thread::sleep_for(std::chrono::seconds(20));
                           // std::cout << "world " << run << std::endl;

                           sampled_e = NULL;
                           sampled = 0;
                           total_steps = 0;
                           while (sampled < (samplebudget - independence_cost*initial_nodes) && njumps < njumps_limit) {
                               // get one walker and set node as visited
                               r = samplednodes[sampled%initial_nodes];
                               r->flag = 1;

                               // ger random edge
                               e = pG->randomedge(r);

                               // add sample and update walker position
                               fsr.add(e,(_FLOAT)pG->degree(e->v));
                               samplednodes[sampled%initial_nodes] = pG->getvertex(e->v);

                               // increment counters
                               if(pG->flag(e->v) != 1)
                               {
                                   sampled++;
                               } else {
                                   sampled += COST_REPEATED_VISIT;
                                   fsr.stats[string("revisits")] += 1.0;
                               }
                               total_steps++;

                               //if(total_steps % 100 == 0)
                               //cout << "run " << run << ", total_steps " << total_steps << endl;
                               if(total_steps >= 2*samplebudget)
                                break;


                           }
                           fsr.stats[string("sampled")] = sampled;
                           fsr.stats[string("revisits")] /= samplebudget;


                       //delete pG;
                       return fsr.consolidate(trueSum, rnd);
                       //return fsr.est_pdf;
                   })
                );
            }
            time(&current_time);
            for(auto && result: results) {
                auto ret = result.get();
                Fsum->consolidate(ret.first, ret.second);
            };
            cout << "Time spent: " << current_time-start_time << endl;


            #ifdef TESTDIST
            cerr << "Dist = " << endl;
            t =0;
            for (k = 0; k < 14; k++) {
                t+=degdist[k];
            }
            for (k = 0; k < 14; k++) {
                cerr << "   d[" << k << "] = " << (_FLOAT)degdist[k]/t << endl;
            }
            #endif
    
    
    //if (VarDist) {
    //    G->edge_iterator_start();
    //    err = 0.0;
    //    while ((e = G->next_edge()) != NULL) {
    //    //assuming the graph is undirected (i.e. there is an inverse edge) ... otherwise it will explode
    //        q += e->p;
    //        err = ABSMAX(err, (abs(2.0/G->no_edges() - (e->p+G->reverseEdge(e)->p) )));
    //    }
    //    cerr << "# 1-Sum of all sampling probabilities = "<<1.0-q<<endl;
    //}
    return 0.0;
  }
};

class RandomWalkSamplingWithJumps {
    private:
        Graph *G;
        
    public:
        RandomWalkSamplingWithJumps(Graph *pG) {
            G = pG;
        }
        
        _FLOAT runsampled(int64_t samplebudget, int32_t independentWalkers, int64_t semiruns, bool VarDist, _FLOAT jump_weight, F *Fsum, int64_t *) {
            map<LabelPi,_FLOAT> trueSum = Fsum->reallyTrueSum();
            if (Fsum != NULL) Fsum->reset();
            
            if (samplebudget % independentWalkers != 0) {
                cerr << "Warning in multiple random walks: Sample budget must be a multiple of the random walkers";
            }

            ThreadPool tpool(nthreads);
            map<thread::id, Graph *> gpool;

            std::vector< std::future<int> > return_vals;
            for (thread::id &id: tpool.get_ids()) {
                return_vals.emplace_back(tpool.enqueue([&gpool,id]() {
                    gpool[id] = new Graph(filename);
                    cout << "Done " << id << " by " << this_thread::get_id() << endl;
                    return 0;
                }));
            }
            for(auto && return_val: return_vals)
                return_val.get();

            time_t start_time, current_time;
            time(&start_time);
            std::vector< std::future< pair<PDF,PDF> > > results;
            for(int run = 0; run < semiruns; ++run) {
                results.emplace_back(
                   tpool.enqueue([&,run] {

                       // get graph from pool
                       Graph *pG = gpool[this_thread::get_id()];
                       pG->zeroflags();
                       pG->reset_undirected_graph();

                       vertex *r, *jumpvertex;
                       edge *e, *rev_e, *sampled_e;
                       edge virtualedge;
                       int64_t sampled;
                       int32_t i;
                       FSResults fsr(pG);
                       fsr.stats[string("budget")] = samplebudget;

                       gsl_rng_env_setup();
                       gsl_rng *rnd = gsl_rng_alloc(gsl_rng_default);
                       gsl_rng_set(rnd, run);
                       pG->setRNG(rnd);

            
                        for (int walker = 0; walker < independentWalkers; walker++) {
                            /*
                            // STEADY STATE START
                            e = pG->randomedge();
                            r = pG->getvertex(e->v);
                            */
                            r = pG->randomvertex();
                            #ifdef DEBUGDURW
                            cout << "initial sample: " << r->id << endl;
                            #endif

                            sampled_e = NULL;
                            sampled = 0;
                            while (sampled < (samplebudget/independentWalkers) - independence_cost) {
                                // Mark vertex "r" as a visited vertex
                                r->flag = 1;

                                // cerr << "(" << r->id << ") => ";
                            
                                //debug_deg = r->undirected_deg;
                            
                                // Mark all edges (r,v) going out of "r" as real edges if "v" was not already visited
                                r->undirected_deg = 0;
                                for (i =0; i < pG->degree(r); i++) {
                                    e = pG->getedge(r,i);
                                    // cerr << " bu";
                                    if ( pG->flag(e->v) == 0 ) {
                                        bool check_rev = false;
                                        if( e->p == NO_UNDIRECTED_GRAPH_EDGE_FLAG ) {
                                            e->p = UNDIRECTED_GRAPH_EDGE_FLAG;
                                            check_rev = true;
                                        } else if( see_incoming_edges && e->p == INVALID_EDGE_FLAG ) {
                                            e->p = ADDED_EDGE_FLAG;
                                            check_rev = true;
                                        }
                                        if( check_rev ) {
                                            rev_e = pG->reverseEdge(e);
                                            if (rev_e == NULL) {
                                                cerr << " reverse edge does not exist (" << e->v << "," << e->u << ")!!"<<endl;
                                                exit(0);
                                            }
                                            if( rev_e->p == NO_UNDIRECTED_GRAPH_EDGE_FLAG ) {
                                                rev_e->p = UNDIRECTED_GRAPH_EDGE_FLAG;
                                            } else if( rev_e->p == INVALID_EDGE_FLAG ) {
                                                rev_e->p = ADDED_EDGE_FLAG;
                                            }
                                        }
                                    }
                                    if( e->p == UNDIRECTED_GRAPH_EDGE_FLAG || e->p == ADDED_EDGE_FLAG ) {
                                        //#ifdef DEBUGDURW
                                        //cout << "\tis connected to " << e->v << endl;
                                        //#endif
                                        r->undirected_deg++;
                                    }
                                }
                                // Adds the previous edge to the sample set. 
                                // This is done like this because we need the "undirected degree of the vertex e->v and this
                                //  we only get after we visit e->v.
                                // WARNING: e cannot be stored as virtualedge may change in the next run
                                //          if (sampled > 0) // avoids sampling the first edge because it may be too biased
                                if ((sampled_e != NULL) && (Fsum != NULL) && (sampled > 0)) { 
                                    // cerr << " [sampled] ";
                                    //Fsum->add(sampled_e,sampled,RWSAMPLING,((_FLOAT)r->undirected_deg+jump_weight));
                               #ifdef DEBUGDURW
                               cout << "Adding " << r->id << ", weight" << r->undirected_deg+jump_weight << endl;
                               #endif
                                    fsr.add(sampled_e,((_FLOAT)r->undirected_deg+jump_weight));
                                }
                                
                                //if (r->undirected_deg != debug_deg) cerr << "# deg("<<r->id<<")="<<r->undirected_deg<<" , old = "<<debug_deg << endl;
                                
                                
                                // Now select next vertex
                                if (pG->myUnifrand() <= (jump_weight/(jump_weight + r->undirected_deg))) {
                                    jumpvertex = pG->randomvertex();
                                    virtualedge.u = r->id;
                                    virtualedge.v = jumpvertex->id;
                                    virtualedge.p = INVALID_EDGE_FLAG;
                                    sampled_e = &virtualedge;
                                    if(jumpvertex->flag !=1)
                                    {
                                        sampled+=independence_cost;
                                    } else {
                                        sampled += COST_REPEATED_VISIT*independence_cost;
                                    }
                                    fsr.stats[string("pjump")] += 1.0/(samplebudget-independence_cost);
                                } else {
                                    do {
                                        e = pG->randomedge(r);
                                    } while ((e->p == NO_UNDIRECTED_GRAPH_EDGE_FLAG) || (e->p == INVALID_EDGE_FLAG));
                                    sampled_e = e;
                                    //cerr << " selected edge (" << e->u << "," << e->v << ")"<<endl;
                                    if(pG->flag(e->v) != 1)
                                    {
                                        sampled++;
                                    } else {
                                        sampled += COST_REPEATED_VISIT;
                                    }
                                }
                                r = pG->getvertex(sampled_e->v);
                               #ifdef DEBUGDURW
                               cout << "Selected edge: " << sampled_e->u << " -> " << sampled_e->v << endl; 
                               #endif
                                //cerr << endl;
                            }
                                
                        }
                        return fsr.consolidate(trueSum, rnd);
                   })
                );
            }

            time(&current_time);
            for(auto && result: results) {
                auto ret = result.get();
                Fsum->consolidate(ret.first, ret.second);
            }
            cout << "Time spent: " << current_time-start_time << endl;

            return -1;
        }
};

class RandomVertexSampling {
 private:
  Graph *G;
 
 public:
  RandomVertexSampling(Graph *pG) {
    G = pG;
  }

  _FLOAT runsampled(int64_t samples, int32_t runs, bool VarDist, F *Fsum) {
    int64_t i, k;
    _FLOAT err = 0.0;
    vertex *r;
    
    //G->zeroedgep(); // what is this?
    if (Fsum != NULL) Fsum->reset();
    
    for (k=0; k < runs; k++) {
        G->zeroflags();
        if (Fsum != NULL) Fsum->prepare_nextrun();
        
        for (i=0; i*independence_cost < samples; i++) {
            r = G->randomvertex();
            if( r->flag == 0 ) {
                r->flag = 1;
            }
            //cout << r->outdeg << endl;
            if (Fsum != NULL) {
                Fsum->add(r,i,VERTEXSAMPLING,1.0); // compare with the optimal + indepedent
            }
        }
        
        if (Fsum != NULL) Fsum->consolidate();
    }
    
    return err;
  }  
};


class RandomEdgeSampling {
 private:
  Graph *G;
 
 public:
  RandomEdgeSampling(Graph *pG) {
    G = pG;
  }

  _FLOAT runsampled(int64_t samples, int32_t runs, bool VarDist, F *Fsum) {
    int64_t i, k;
    edge *e,*ie;
    _FLOAT err = 0.0;
    _FLOAT q;
    
    G->zeroedgep();
    if (Fsum != NULL) Fsum->reset();
    
    for (k=0; k < runs; k++) {
        if (Fsum != NULL) Fsum->prepare_nextrun();
        
        for (i=0; i*independence_cost < samples; i++) {
            e = G->randomedge();
            e->p += ((_PTYPE)1.0);
            if (Fsum != NULL) Fsum->add(e,i,EDGESAMPLING,0); // compare with the optimal + indepedent
        }
        
        if (Fsum != NULL) Fsum->consolidate();
    }
    if (VarDist) {
        G->edge_iterator_start();
        err = 0;
        q = 0;
        while ((e = G->next_edge()) != NULL) {
        //assuming the graph is undirected (i.e. there is an inverse edge) ... otherwise it will explode
            ie = G->reverseEdge(e);
            if (e->p > 1)  e->p /= (_PTYPE)samples;
            if (ie->p > 1)  ie->p /= (_PTYPE)samples;
            err = ABSMAX(err, (abs(2.0/G->no_edges() - (e->p+ie->p))));
            q += e->p;
        }
        cerr << "# 1-Sum of all sampling probabilities = "<<1.0-q<<endl;
    }
    
    return err;
  }  
};

class GraphStatistics {
    private:
        Graph *G;
        int32_t flags[256];
        
    public:
        GraphStatistics(Graph *pG) {
            G = pG;
            memset(flags,0,sizeof(flags));
        }
        
        void flood_flag(vertex *u,_FLAGTYPE f) {
            int32_t i;
            vertex *v;
            stack<vertex *> s;
            
            s.push(u);
            while (!s.empty()) { 
                u = s.top();
                s.pop();
                G->setflag(u,f);
                for (i=0; i < u->deg; i++) {
                    v = G->getneighbor(u,i);
                    if (G->flag(v) == 0) {
                        s.push(v);
                    }
                    else {
                        if (G->flag(v) != f) {
                            cerr << "Hit flag "<<(int)G->flag(v)<< "  using flag "<<(int)f<<endl;
                        }
                    }
                }
            }
        }
        
        void flag_connected_components() {
            vertex *v = NULL;
            _FLAGTYPE flag;
            
            G->zeroflags();
            G->vertex_iterator_start();
            flag = 1;
            while ((v = G->next_vertex()) != NULL) {
                if (G->flag(v) == 0) {
                    flood_flag(v,flag);
                    if (flag < 255) flag++;
                }
            }
        }
        
        void print_connected() {
            int i;
            vertex *v = NULL;
            
            memset(flags,0,sizeof(flags));
            G->vertex_iterator_start();
            while ((v = G->next_vertex()) != NULL) {
                flags[G->flag(v)]++;
            }
            for (i=0; i < 256; i++) {
                if (flags[i] != 0) {
                    cerr << "flag["<<i<<"]="<<(int)flags[i]<<endl;
                }
            }
        }
        
        inline int32_t count_flag(_FLAGTYPE flag) {
            return flags[flag];
        }
        
        _FLAGTYPE largest_componentflag() {
            vertex *v = NULL;
            _FLAGTYPE flag=0;
            int32_t maxcount = 0;
            memset(flags,0,sizeof(flags));
            G->vertex_iterator_start();
            while ((v = G->next_vertex()) != NULL) {
                flags[G->flag(v)]++;
                if (G->flag(v) < 255) {
                    if (flags[G->flag(v)] > maxcount) {
                        maxcount = flags[G->flag(v)];
                        flag = G->flag(v);
                    }
                }
            }
            return flag;    
        }
        
        void graph_F(F *Fsum) {
            edge *e;
            
            if (Fsum != NULL) Fsum->reset();            
            if (Fsum != NULL) Fsum->prepare_nextrun();
            
            G->edge_iterator_start();
            while ((e = G->next_edge()) != NULL) {
                if (Fsum != NULL) Fsum->add_all(e, EDGESAMPLING, G->degree(e->u), G->degree(e->v));
            }
            if (Fsum != NULL) Fsum->consolidate_TrueSum();
        }

        void graph_mu(F *Fsum) {
          vertex *v;
          map<LabelPi,int> counts;
          //_FLOAT TwoM = 0.0;

          G->vertex_iterator_start();
          while ((v = G->next_vertex()) != NULL) {
              // get undirected degree pi
              // _FLOAT pi = 0;
           //  for (int i =0; i < G->degree(v); i++) {
           //    edge *e = G->getedge(v,i);
           //    if ((e->p == INVALID_EDGE_FLAG) || (e->p == NO_UNDIRECTED_GRAPH_EDGE_FLAG)) {
           //      pi += 1.0;
           //    }
           //  }
           //  TwoM += pi;

            // generate key
            stringstream s;
            if( marginal == GET_OUTDEGREE )
                s << G->outdegree(v);
            else if( marginal == GET_INDEGREE )
                s << G->indegree(v);
            else if( marginal == GET_BOTH )
                s << G->indegree(v) << '\t' << G->outdegree(v);
            mu_star[s.str()] += 1.0;
          }

          // // compute un-normalized mu_star
          // for(auto it = counts.begin(); it != counts.end(); it++) {
          //     string label = it->first.first;
          //     _FLOAT pi = it->first.second;
          //     mu_star[label] += it->second/pi;
          // }

          // // normalize mu_star
          // auto norm = 0.09;
          // cout << "TwoM: " << TwoM << ", no_edges: " << G->no_edges() << endl;
          // for(auto it = mu_star.begin(); it != mu_star.end(); it++)
          //     mu_star[it->first] *= norm;

          // cout << "Normalized mu:" << mu << endl;

        }
};

class RunningStat
{
    public:
        RunningStat() : m_n(0) {}

        void Clear()
        {
            m_n = 0;
        }

        void Push(double x)
        {
            m_n++;

            // See Knuth TAOCP vol 2, 3rd edition, page 232
            if (m_n == 1)
            {
                m_oldM = m_newM = x;
                m_oldS = 0.0;
            }
            else
            {
                m_newM = m_oldM + (x - m_oldM)/m_n;
                m_newS = m_oldS + (x - m_oldM)*(x - m_newM);
    
                // set up for next iteration
                m_oldM = m_newM; 
                m_oldS = m_newS;
            }
        }

        int NumDataValues() const
        {
            return m_n;
        }

        double Mean() const
        {
            return (m_n > 0) ? m_newM : 0.0;
        }

        double Variance() const
        {
            return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
        }

        double StandardDeviation() const
        {
            return sqrt( Variance() );
        }

    private:
        int m_n;
        double m_oldM, m_newM, m_oldS, m_newS;
};


class Fmap : public F {
    protected:
        
        map<LabelPi,_FLOAT> val_e;
        map<string,_FLOAT> val_v;
        map<string,_FLOAT> X, VARTX;
        map<string,RunningStat> stats;
        int32_t sizeofVprime, sizeofEprime;
        int32_t samplepoints;
        F* TrueSum;
        bool tru_ccdf_computed;

    public:

        
        Fmap(Graph *pG, F* pTrueSum) : F(pG) {
            sizeofVprime = 0;
            sizeofEprime = 0;
            samplepoints = 0;
            TrueSum = pTrueSum;
            tru_ccdf_computed = false;
        }
        
        virtual ~Fmap() {
        }
        
        virtual inline _FLOAT g(edge *e) {
            return (((_FLOAT)1.0)/G()->degree(e->v));
        }
        
        virtual void reset() {
            prepare_nextrun();
            samplepoints = 0;
        }
        
        virtual void prepare_nextrun() {
            val_v.clear();
            val_e.clear();
            sizeofEprime = 0;
            sizeofVprime = 0;
        }
        
        /*virtual void add(edge *e, int32_t ) {
            val_e[G()->degree(e->v)] += g_v(e);
            sizeofEprime++;
        }
        */
        virtual void add(edge *e, int32_t , int8_t flag, _FLOAT ) {
        }
        virtual void add(vertex *e, int32_t , int8_t flag, _FLOAT ) {
        }
        
        virtual void consolidate_TrueSum() {
        }
        
        virtual  map<LabelPi,_FLOAT> trueSum() {
            return val_e;
        }

        virtual  map<LabelPi,_FLOAT> reallyTrueSum() { // workaround
            return TrueSum->trueSum();
        }

        /*virtual void setMu(_FLOAT _mu) { mu = _mu; }
        virtual _FLOAT getMu() { return mu; }*/
        virtual void consolidate(PDF est_pdf, PDF new_stats) {
            _FLOAT TrueSumg = normalize_TrueSum();
            map<LabelPi,_FLOAT> fvalue = TrueSum->trueSum();

            long double norm = 0.0;
            for( map<string,_FLOAT>::iterator it = est_pdf.begin(); it != est_pdf.end(); it++ )
              norm += it->second;

            for( map<LabelPi,_FLOAT>::iterator it = fvalue.begin(); it != fvalue.end(); it++ ) {
                string label = it->first.first;
                double est = est_pdf[label] / norm;
                // cout << "B " << label << "\t" << est_pdf[label] << endl;
                //cout << "est_pdf[" << degree << "]=" << est_pdf[degree] << endl;
                X[label] += est;
                VARTX[label] += (est -  it->second/TrueSumg) * (est -  it->second/TrueSumg); 
            }


            for(auto const stat:new_stats) {
                cout << "Stat " << stat.first << ": " << stat.second << endl;
                stats[stat.first].Push(stat.second);
            }

            samplepoints++;
        }
        
        virtual void consolidate() {
            map<string,_FLOAT> est_pdf;
            _FLOAT TrueSumg = normalize_TrueSum();
            map<LabelPi,_FLOAT> fvalue = TrueSum->trueSum();

            map<string,_FLOAT> mu; // normalization constant
            map<string,_FLOAT> m;
            for(map<LabelPi,_FLOAT>::iterator it = val_e.begin(); it != val_e.end(); it++) {
                m[it->first.first] += it->second;
                mu[it->first.first] += it->second/it->first.second; // mu[i] = \sum_j m_{i,j}/j
            }
            for(auto it=mu_star.begin(); it != mu_star.end(); it++)
                fprintf(stdout, "%s\t%Lf\t%Lf\n", it->first.c_str(), mu[it->first], it->second);

            if(reduce_var)
                  for( map<LabelPi,_FLOAT>::iterator it = fvalue.begin(); it != fvalue.end(); it++ ) {
                      string label = it->first.first;
                      if(mu.find(label) == mu.end()) {
                          sizeofVprime -= val_v[label];
                          val_v.erase(label);
                      }
                  }

            if(avg_given) {

              _FLOAT mu0 = 0.0;
              for( map<LabelPi,_FLOAT>::iterator it = fvalue.begin(); it != fvalue.end(); it++ ) {
                  string label = it->first.first;
                  if(mu.find(label) != mu.end()) {
                      //fprintf(stdout, "%s: %Lf * %Lf / %Lf\n", label.c_str(), it->second/TrueSumg, m[label], mu[label]);
                      mu0 += it->second/TrueSumg * m[label] / mu[label];
                  }
              }



              // compute estimator assuming average is known
              for( map<LabelPi,_FLOAT>::iterator it = fvalue.begin(); it != fvalue.end(); it++ ) {
                  string label = it->first.first;
                  if(mu.find(label) == mu.end()) {
                      est_pdf[label] = val_v[label]/sizeofVprime;
                  } else {
                      est_pdf[label] = (val_v[label]+m[label])/
                          (sizeofVprime + sizeofEprime/mu0 *m[label]/mu[label] );
                  }
                  //fprintf(stdout, "%s: (%Lf+%Lf)/(%d + %d/%Lf * %Le/%Le)\n",
                  //        label.c_str(), val_v[label], m[label], sizeofVprime, sizeofEprime, mu0, m[label], mu[label]);
              }


            }
            else {
              // TODO: remove this
              // Compute correct normalizing constant
              _FLOAT mu0 = 0.0;
              for( map<LabelPi,_FLOAT>::iterator it = fvalue.begin(); it != fvalue.end(); it++ ) {
                  string label = it->first.first;
                  if(mu.find(label) != mu.end()) {
                      //fprintf(stdout, "%s: %Lf * %Lf / %Lf\n", label.c_str(), it->second/TrueSumg, m[label], mu[label]);
                      mu0 += it->second/TrueSumg * m[label] / mu[label];
                  }
              }

              // call gradient descent method to compute estimates
              est_pdf = dufs::maxLikelihoodFS(val_v, m, mu);

              // TODO: REMOVE THIS
              _FLOAT mup = 0.0;
              for( map<string,_FLOAT>::iterator it = est_pdf.begin(); it != est_pdf.end(); it++ ) {
                  string label = it->first;
                  if(mu.find(label) != mu.end()) {
                      //fprintf(stdout, "%s: %Lf * %Lf / %Lf\n", label.c_str(), it->second/TrueSumg, m[label], mu[label]);
                      mup += it->second * m[label] / mu[label];
                  }
              }
              fprintf(stdout, "mu0: %Lf, mup: %Lf\n", mu0, mup);
            }

            long double norm = 0.0;
            for( map<string,_FLOAT>::iterator it = est_pdf.begin(); it != est_pdf.end(); it++ )
            {
              // cout << "A " << it->first << "\t" << it->second << endl;
              norm += it->second;
            }
            cout << "Norm:" << norm << endl;


            for( map<LabelPi,_FLOAT>::iterator it = fvalue.begin(); it != fvalue.end(); it++ ) {
                string label = it->first.first;
                double est = est_pdf[label] / norm;
                // cout << "B " << label << "\t" << est_pdf[label] << endl;
                //cout << "est_pdf[" << degree << "]=" << est_pdf[degree] << endl;
                X[label] += est;
                VARTX[label] += (est -  it->second/TrueSumg) * (est -  it->second/TrueSumg); 
            } 


            samplepoints++;
        }
        
        virtual _FLOAT normalize_val(map<LabelPi,_FLOAT> &val_e) {
            cerr << "Fmap::normalize_val()" << endl;
            return ((_FLOAT)sizeofEprime);
        }
        
        virtual _FLOAT normalize_TrueSum() {
            map<LabelPi,_FLOAT> tmp;
            cerr << "Fmap::normalize_TrueSum()" << endl;
            return TrueSum->normalize_val(tmp);
        }    

        virtual void print(const char *preamble) {
            _FLOAT TrueSumg;
            
            
            TrueSumg = normalize_TrueSum();
            cerr << "#" << preamble << "\tE[\\hat{X}]=\\sum_{k=1}^K\\hat{X}_k/K"<< endl;

            string marginal_str;
            if( GET_OUTDEGREE )
                marginal_str = "outdeg";
            else if( GET_INDEGREE )
                marginal_str = "indeg";
            else
                marginal_str = "indeg\toutdeg";

            cerr << "#";
            cerr.precision( 3 );
            for(auto const stat:stats) {
                cerr << "\t" << stat.first << " " << stat.second.Mean();
                double sd = stat.second.StandardDeviation();
                if(sd > 0.0)
                    cerr << "±" << sd;
            }
            cerr << endl;

            cerr << "#" << preamble << "\t" << marginal_str << "\tE[\\hat{X}]\tVar[\\hat{X}]\tE[X]\t1-E[\\hat{X}]/E[X]\tsqrt(E[(\\hat{X}-E[X])^2])"<< endl;
            cerr.precision( 15 );
            cerr.setf(ios::scientific,ios::floatfield);

            #ifdef PRINT_PMF
            map<LabelPi,_FLOAT> fvalue = TrueSum->trueSum();
            for( map<LabelPi,_FLOAT>::iterator it = fvalue.begin(); it != fvalue.end(); it++ ) {
                string label = it->first.first;
                cerr << "PMF " << "\t" << label << "\t" << X[label]/samplepoints << "\t" <<
                   VARTX[label]/(samplepoints*samplepoints) << "\t" << it->second/TrueSumg << "\t" <<
                   1.0-(X[label]/samplepoints)/(it->second/TrueSumg)<< "\t" << sqrt(VARTX[label]/samplepoints) << endl;
            }
            #endif
            cerr << endl;
          }
          virtual void printrealvalue(const char *preamble) {
            //_FLOAT TrueSumg;    
            //TrueSumg = normalize_TrueSum();

           // cerr << "#" << preamble << "\tE[CCDF\\hat{X}]=\\sum_{k=1}^K\\hat{X}_k/K"<< endl;
           // cerr << "#" << preamble << "\ti\tE[CCDF\\hat{X}]\t\tE[CCDF_X]\t\t1-E[CCDF\\hat{X}]/E[CCDF_X]\tsqrt(E[(CCDF\\hat{X}-E[CCDF_X])^2])"<< endl;
           // for (i =0; i <= G()->MAXDEGREE(); i++) {
           //         if (X[i] > 0)
           //                 cerr << preamble << "\t" << i << "\t" << X[i]/samplepoints << "\t" << (TrueSum->trueSum()[i]) << "\t" << 1.0-(X[i]/samplepoints)/(TrueSum->trueSum()[i])<< "\t" << sqrt(VARTX[i]/samplepoints) << endl;
           // }
           //  cerr << preamble << "\t" << (TrueSum->trueSum()[1]/TrueSumg)<< "\t" << (TrueSum->trueSum()[5]/TrueSumg)<< "\t" << (TrueSum->trueSum()[10]/TrueSumg)<< "\t" << (TrueSum->trueSum()[100]/TrueSumg)<< "\t" << (TrueSum->trueSum()[500]/TrueSumg)<< "\t" << (TrueSum->trueSum()[1000]/TrueSumg)<< endl;
            cerr.flush();
          }
        
};



class FBiasedDegDist : public Fmap {
    
  public:        
    FBiasedDegDist(Graph *pG, F* pTrueSum) : Fmap(pG,pTrueSum) {
    }

    virtual void add(edge *e, int32_t , int8_t flag, _FLOAT pi) {
        if (pi > 0.0) {
            stringstream s;
            if( marginal == GET_OUTDEGREE )
                s << G()->outdegree(e->v);
            else if( marginal == GET_INDEGREE )
                s << G()->indegree(e->v);
            else if( marginal == GET_BOTH )
                s << G()->indegree(e->v) << '\t' << G()->outdegree(e->v);

            auto tmp = make_pair(s.str(),pi);
            val_e[tmp] += 1.0;
            sizeofEprime++;

            //cout << s.str() << ' ';
        }
    }

    virtual void add(vertex *v, int32_t , int8_t flag, _FLOAT pi) {
        if (pi > 0.0) {
            stringstream s;
            if( marginal == GET_OUTDEGREE )
                s << v->outdeg;
            else if( marginal == GET_INDEGREE )
                s << v->indeg;
            else if( marginal == GET_BOTH )
                s << v->indeg << '\t' << v->outdeg;

            auto tmp = make_pair(s.str(),pi);
            val_v[s.str()] += 1.0;
            sizeofVprime++;

            //cout << s.str() << ' ';

        }
    }
    
    virtual void add_all(edge *e, int8_t flag, _FLOAT pi_u, _FLOAT pi_v) {
        //cerr << "pi("<<e->u<<") = "<<pi_u << " pi("<<e->v<<") = "<<pi_v<< endl;

        stringstream s;
        if( marginal == GET_OUTDEGREE )
            s << G()->outdegree(e->u);
        else if( marginal == GET_INDEGREE )
            s << G()->indegree(e->u);
        else if( marginal == GET_BOTH )
            s << G()->indegree(e->u) << '\t' << G()->outdegree(e->u);
        //cerr << s.str() << endl;
        val_e[make_pair(s.str(),-1.0)] += 1.0/pi_u;

        s.clear();
        s.str(std::string());
        if( marginal == GET_OUTDEGREE )
            s << G()->outdegree(e->v);
        else if( marginal == GET_INDEGREE )
            s << G()->indegree(e->v);
        else if( marginal == GET_BOTH )
            s << G()->indegree(e->v) << '\t' << G()->outdegree(e->v);
        //cerr << s.str() << endl;
        val_e[make_pair(s.str(),-1.0)] += 1.0/pi_v;

        sizeofEprime+=2;
    }
    
    virtual _FLOAT normalize_TrueSum() {
        _FLOAT TrueSumg;
                
        //cerr << "FBiasedDegDist::normalize_TrueSum()" << endl;
        TrueSumg = 0.0;
        map<LabelPi,_FLOAT> tmp = TrueSum->trueSum();
        for(map<LabelPi,_FLOAT>::iterator it = tmp.begin(); it != tmp.end(); it++ )
            TrueSumg += it->second;
        return TrueSumg;
    }

    virtual _FLOAT normalize_val(map<LabelPi,_FLOAT> &fvalue) {
        _FLOAT C;
        
        //cerr << "FBiasedDegDist::normalize_val()" << endl;
        C = 0.0;
        for(map<LabelPi,_FLOAT>::iterator it = fvalue.begin(); it != fvalue.end(); it++ )
            C += it->second;
        return C;
    }
    
};



#define VARIATIONDISTANCE false


int main(int argc, char *argv[]) {
    signal(SIGSEGV, handler);


    Graph *G;
    GraphStatistics *GS;
    //RandomEdgeSampling *RE;
    FrontierSampling *FS;
    RandomWalkSampling *RW;
    RandomWalkSamplingWithJumps *RWJ;
    //Fdegfreq *Fsum;
    FBiasedDegDist *Fsum;
    FBiasedDegDist *FTrueSum;
    
    // local variables
    float budget;
    int64_t initial_nodes;
    int32_t semiruns;
    int64_t walkers;
    _FLOAT RW_jump_weight;

    // global variables
    //   filename
    //   independence_cost
    //   marginal
    //   nthreads
    //   see_incoming_edges
    

    int bflag, cflag, dflag, fflag, lflag, nflag, rflag, tflag, vflag, mflag, wflag;
    bflag = cflag = dflag = fflag = lflag = nflag = rflag = tflag = vflag = mflag = wflag = 0;
    while (1)
    {
      static struct option long_options[] =
        {
          /* These options set a flag. */
          {"avg_given",     no_argument,       &avg_given,  1},
          {"not_recursive", no_argument,       &not_recursive,  1},
          {"reduce_var",    no_argument,       &reduce_var, 1},
          {"use_hybrid",    no_argument,       &use_hybrid, 1},
          {"use_durw",      no_argument,       &use_durw,   1},
          {"use_mrw",       no_argument,       &use_mrw,   1},
          {"revisit_nocost", no_argument,      &COST_REPEATED_VISIT,   0},
          {"propto_degree", no_argument,       &propto_degree,   1},
          {"inv_propto_degree", no_argument,   &inv_propto_degree,   1},
          /* These options don’t set a flag.
             We distinguish them by their indices. */
          {"budget",       required_argument, 0, 'b'},
          {"indep_cost",   required_argument, 0, 'c'},
          {"distribution", required_argument, 0, 'd'},
          {"filename",     required_argument, 0, 'f'},
          {"njumps_limit", required_argument, 0, 'l'},
          {"nodes",        required_argument, 0, 'n'},
          {"runs",         required_argument, 0, 'r'},
          {"threads",      required_argument, 0, 't'},
          {"visible_in",   required_argument, 0, 'v'},
          {"walkers",      required_argument, 0, 'm'},
          {"jump_weight",  required_argument, 0, 'w'},
          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;
      int c = getopt_long (argc, argv, "b:c:d:f:l:n:r:t:v:m:w:",
                       long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1)
        break;

      switch (c)
        {
        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0)
            break;
          printf ("option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;

        case 'b':
          budget = atof(optarg);
          bflag = 1;
          break;

        case 'c':
          independence_cost = atoi(optarg);
          cflag = 1;
          break;

        case 'd':
          if( strcmp(optarg,"out") == 0 )
              marginal = GET_OUTDEGREE;
          else if( strcmp(optarg,"in") == 0 )
              marginal = GET_INDEGREE;
          else if( strcmp(optarg,"jnt") == 0 )
              marginal = GET_BOTH;
          else {
              cerr << "Marginal distribution " << optarg << "is invalid." << endl;
              exit(1);
          }
          dflag = 1;
          break;

        case 'f':
          filename = new char[strlen(optarg)+1];
          strcpy(filename,optarg);
          fflag = 1;
          break;

        case 'l':
          njumps_limit = atoi(optarg);
          lflag = 1;
          break;

        case 'n':
          initial_nodes = atoi(optarg);
          nflag = 1;
          break;

        case 'r':
          semiruns = atoi(optarg);
          rflag = 1;
          break;

        case 't':
          nthreads = atoi(optarg);
          tflag = 1;
          break;

        case 'v':
          see_incoming_edges = atoi(optarg);
          vflag = 1;
          break;

        case 'm':
          walkers = atoi(optarg);
          mflag = 1;
          break;

        case 'w':
          RW_jump_weight = atof(optarg);
          wflag = 1;
          break;

        case '?':
          /* getopt_long already printed an error message. */
          break;

        default:
          abort ();
        }
    }

    // mandatory flags
    if(bflag+dflag+fflag+rflag+vflag+mflag != 6) {
        if(!bflag) printf("--budget mandatory\n");
        if(!dflag) printf("--distribution mandatory\n");
        if(!fflag) printf("--filename mandatory\n");
        if(!rflag) printf("--runs mandatory\n");
        if(!vflag) printf("--visible_in mandatory\n");
        if(!mflag) printf("--walkers mandatory\n");
        abort();
    }
    if(use_durw && use_mrw) {
       printf("incompatible flags: use_durw and use_mrw\n");
    }
    if(propto_degree && inv_propto_degree) {
       printf("incompatible flags: propto_degree and inv_propto_degree\n");
    }

    // optional flags
    if(!cflag) independence_cost = 1;
    if(!nflag) initial_nodes = walkers;
    if(!tflag) nthreads = 1;
    if(!wflag) RW_jump_weight = 0.0;
    if(!lflag) njumps_limit = numeric_limits<int>::max();

    //if (argc < 8 || argc > 12) {
    //    cerr << "Missing parameters:\n" <<
    //    "   <filename>: trace file name\n" <<
    //    "   <steps>: number of total samples ((steps-walkers) in random walks)\n" <<
    //    "   <walkers>: number of walkers in random walks)\n" <<
    //    "   <nodes>: number of initial nodes)\n" <<
    //    "   <runs>: number of runs\n" <<
    //    "   <independence cost>: cost of sampling independently\n" <<
    //    "   <nthreads>: number of threads\n" <<
    //    "   [<see incoming edges> = 1 [out | in] ]: 0 (false) or 1 (true, default)\n" <<
    //    "   [<average degree is given> = 0]: 0 (false, default) or 1 (true)" <<
    //    "   [<reduce variance trick> = 0]: 0 (false, default) or 1 (true)" << endl;
    //    exit(1);
    //}


    //cerr.precision( 15 );
    //cerr.setf(ios::scientific,ios::floatfield);
    
    
    cerr << "#Starting reading graph file...\n";
    G = new Graph(filename);
    cerr << "#Finished reading graph file...\n";
    cerr << "#Number of edges = " << G->no_edges() << endl;
    cerr << "#Number of vertices = " << G->size() << endl;
    
    cerr << "#budget = " << budget << endl;
    int64_t steps_samples =int(G->size() * budget);
    cerr << "#steps_samples = " << steps_samples << endl;

    if(walkers > steps_samples) {
        cout << "Warning: nwalkers > budget. Setting nwalkers = budget." << endl;
        walkers = steps_samples;
    }

   // cerr << "maxOutdegree:"<<maxOutdegree<<endl;
    cerr << "#Filename: "<<filename<<", <steps> = "<<steps_samples<<", <walkers> = "<<walkers<<", <initial nodes> = " << initial_nodes << ", <independence cost> = " << independence_cost << ", runs = "<<semiruns<<endl;
    GS = new GraphStatistics(G);
    
   // GS->flag_connected_components();
   // flag = GS->largest_componentflag();
   // cerr << "#Largest component has flag = "<< (int)flag <<" and "<< GS->count_flag(flag) << " nodes. The whole graph has " << G->size() << " nodes"<< endl;
   // if (largcomp == 1) {
   //     G->trimGraphDifferent(flag);
   //     cerr << "# Trimming graph "<<endl;
   // }
    
    //  estimates the distribution from random sampling
    // RE = new RandomEdgeSampling(G);
    //RW = new RandomWalkSampling(G);
    
    // Fsum = new Fdegfreq(G);
    FTrueSum = new FBiasedDegDist(G,NULL);
    Fsum = new FBiasedDegDist(G,FTrueSum);
    GS->graph_F(FTrueSum);
    GS->graph_mu(FTrueSum);

    _FLOAT p;
    int64_t samplededges;

    if(use_durw) {
        RWJ = new RandomWalkSamplingWithJumps(G);

        cerr << "#RWSW_withJumps: computing f(G) with ("<<semiruns<<" runs,"<<steps_samples<<" steps, 1 walker)"<<endl;
        p = RWJ->runsampled(steps_samples,walkers,semiruns,false,RW_jump_weight,Fsum,&samplededges);
        Fsum->print("RWSW_wJ ");
        cerr.flush();

        delete RWJ;
    } else if(use_mrw) {
        RW = new RandomWalkSampling(G);
        p = RW->runsampled(steps_samples,walkers,semiruns,false,Fsum,&samplededges);
        Fsum->print("RWMW ");
        cerr.flush();

        delete RW;
    } else {
        FS = new FrontierSampling(G);

        cerr << "#FS: computing f(G) with ("<<semiruns<<" runs,"<<steps_samples<<" steps,"<< walkers <<" walkers)"<<endl;    
        if(!wflag) RW_jump_weight = numeric_limits<double>::min();
        p = FS->runsampled(steps_samples,walkers,initial_nodes,semiruns,false,RW_jump_weight,Fsum,&samplededges);
        Fsum->print("FS ");
        cerr.flush();

        delete FS;
    }


    /*cerr << "RWMW: computing f(G) with ("<<semiruns<<" runs,"<<steps_samples<<" steps,"<< walkers <<" walkers)"<<endl;    
    p = RW->runsampled(steps_samples,walkers,semiruns,false,Fsum,&samplededges);
    Fsum->print("RWMW ");
    cerr.flush();*/


    /*cerr << "#RWSW: computing f(G) with ("<<semiruns<<" runs,"<<steps_samples<<" steps, 1 walker)"<<endl;    
    p = RW->runsampled(steps_samples,1,semiruns,false,Fsum,&samplededges);
    Fsum->print("RWSW ");
    cerr.flush();*/

    /*p = RE->runsampled(steps_samples,semiruns,false,Fsum);
    Fsum->print("RE ");
    cerr.flush();*/
        
    delete Fsum;
    delete FTrueSum;
    
    // delete RE;
    // delete FS;
    // delete RW;
    if (GS != NULL) delete GS;
    delete G;
    
    return 0;
}

