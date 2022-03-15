#define USE_SRAND true
#define USE_LOG_SAMPLE true

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <time.h>
#include <numeric>

using namespace std;


template <class T>
T get_sum(T* al, int size){
  T sum = 0;
  for (int i = 0; i < size; i++){
    sum += al[i];
  }
  return sum;
}

template <class T>
T get_sum_vec(vector<T>  in){
  typedef typename vector<T>::iterator myiterator;
  T sum = 0;
  for(myiterator it =in.begin(); it != in.end() ; ++it){
    sum += *it;
  } 
  return sum;
}

template <class T>
T get_max(T* array, int N){
  T max = 0;
  for (int n = 0; n < N; n ++){
    if(array[n] > max){
      max = array[n];
    }
  }
  return max;
}

template <typename T>
T get_max_vec(vector<T> in){
  typedef typename vector<T>::iterator myiterator;
  T max = 0;
  for(myiterator it =in.begin(); it != in.end() ; ++it){
    if (*it > max){
        max = *it;
    }
  } 
  return max;
}

int sample_from_mult_log(vector< double > pval, int SEED){
  vector<double> pval_cul_log;
  vector<double> pval_cul;
  pval_cul_log.push_back(log(pval[0]));
  pval_cul.push_back(pval[0]);
  int K = pval.size();
  for (int k = 1; k < K; k ++){
    pval_cul.push_back(pval[k] + pval_cul[k-1]); 
    pval_cul_log.push_back( log(pval[k] + pval_cul[k]) );
  }

  double u_log = log( ((double) random() /RAND_MAX) );
  unsigned int picked;
  for(picked = 0; picked < pval_cul_log.size() ; ++picked){
      if(pval_cul_log[picked] >= u_log){
        break;
      }
  }
  return picked;
}

int sample_from_mult(vector< double > pval, int SEED){
  vector<double> pval_cul;
  pval_cul.push_back(pval[0]);
  int K = pval.size();
  for (int k = 1; k < K; k ++){
    pval_cul.push_back(pval[k] + pval_cul[k-1]); 
  }
  double u = ((double) random() /RAND_MAX)*pval_cul[K-1];
  int picked;
  for(picked = 0; picked < K ; ++picked){
      if(pval_cul[picked] >= u){
        break;
      }
  }
  return picked;
}

template <class T>
vector< vector< T> > convert_flat_to_2d(T* in,  int ROW, int COL){
  vector< vector< T> > out;
  for (int r = 0; r < ROW; r ++){
    vector<T> out_t;
    for (int c = 0; c < COL; c ++){
      out_t.push_back(in[r*COL + c]);
    }
    out.push_back(out_t);
  }
  return out;
}

template <typename T>
void convert_2d_to_flat(vector< vector<T> > in, T* out){
  for(unsigned int i = 0;i < in.size(); i++){
    for(unsigned int j = 0;j < in[i].size(); j++){
      out[i*in[i].size() + j] = in[i][j];
    }
  }
}

template <class T>
void convert_3d_to_flat(vector< vector< vector<T> > > in, T* out){
  for(unsigned int i = 0;i < in.size(); i++){
    for(unsigned int j = 0;j < in[i].size(); j++){
      for(unsigned int k = 0;k < in[i][j].size(); k++){
        out[i*in[i].size()*in[i][j].size() + j*in[i][j].size()+k ] = in[i][j][k];
      }
    }
  }
}


template <class T>
vector< vector< vector< T> > >  convert_flat_to_3d(T* in, int ROW, int COL, int K){
  vector< vector< vector< T> > > out;
  for (int r = 0; r < ROW; r ++){
    vector < vector<T> > out_t;
    for (int c = 0; c < COL; c ++){
      vector<T>  out_tt;
      for (int k = 0; k < K; k ++){
        // if array has negative, that means it is
        // irregular array, no elements in there.
        out_tt.push_back(in[r*COL*K + c*K + k]);
      }
      out_t.push_back(out_tt);
    }
    out.push_back(out_t);
  }
  return out;
}

template <typename T>
vector<T> normalize_vec(vector<T> in){
  typedef typename vector<T>::iterator myiterator;
  T sum = 0;
  vector<T> out;
  for(myiterator it =in.begin(); it != in.end() ; ++it){
    sum += *it;
  } 
  for(myiterator it =in.begin(); it != in.end() ; ++it){
    out.push_back( (*it)/sum );
  } 
  return out;
}

vector<double> getAlVal_vec(int psf, int wsf, double lamb, double PEAKSIZE, int len_fg){
  vector<double> out;
  for(int fg=0; fg<len_fg; fg++){
      out.push_back(lamb);
  };
  if (wsf == 1){
    out[psf] = lamb*PEAKSIZE;
  }
  out = normalize_vec(out);

  for(int fg=0; fg<len_fg; fg++){
      out[fg] = out[fg]*lamb*len_fg;
  };
  return out;
};


// now normalize this param using log trick.
vector<double> normalizeLog(vector<double> log_multParam){
    vector<double> exp_log_multParam_minus_c;
    double log_sum_exp_log_multParam;
    double c = get_max_vec(log_multParam);
    exp_log_multParam_minus_c = log_multParam;
    vector<double> multParam;
    int S = log_multParam.size();
    for(int s=0; s < S; ++s){
      exp_log_multParam_minus_c[s] = exp( exp_log_multParam_minus_c[s] - c );
    }
    log_sum_exp_log_multParam = log( get_sum_vec( exp_log_multParam_minus_c ));
    for(int s=0; s < S; ++s){
      multParam.push_back( exp( log_multParam[s] - c - log_sum_exp_log_multParam));
    }
    return multParam;
}

// now normalize this param using log trick.
vector<double> normalizeLog_output_log(vector<double> log_multParam){

    vector<double> exp_log_multParam_minus_c;
    double log_sum_exp_log_multParam;
    double c = get_max_vec(log_multParam);
    exp_log_multParam_minus_c = log_multParam;
    vector<double> multParam_log_normalized;
    int S = log_multParam.size();
    for(int s=0; s < S; ++s){
      exp_log_multParam_minus_c[s] = exp( exp_log_multParam_minus_c[s] - c );
    }
    log_sum_exp_log_multParam = log( get_sum_vec( exp_log_multParam_minus_c ));
    for(int s=0; s < S; ++s){
      multParam_log_normalized.push_back( log_multParam[s] - c - log_sum_exp_log_multParam);
    }
    return multParam_log_normalized;
}

extern "C" {
      void normalize(double* testin, int size, double* out){
        double sum = 0;
        for(int i = 0; i < size; i++){
          sum += testin[i];
        }
        for(int i = 0; i < size; i++){
          out[i] = testin[i]/sum;
        }
      };


      void get_al_val(int psf, int wsf, double lamb, double PEAKSIZE, int len_fg, double* out){
        for(int fg=0; fg<len_fg; fg++){
            out[fg] = lamb;
        };
        if (wsf == 1){
          out[psf] = lamb*PEAKSIZE;
        }
        normalize(out, len_fg, out);

        for(int fg=0; fg<len_fg; fg++){
            out[fg] = out[fg]*lamb*len_fg;
        };
      };

      void log_Beta(double* al, int al_size, double* out){
        double sum_al = get_sum(al, al_size);
        double sum_lgamma_al = 0;
        for(int i = 0; i < al_size; i++ ){
             sum_lgamma_al += lgamma(al[i]);
        }
        out[0] = -lgamma(sum_al) + sum_lgamma_al;
      };

      double log_Beta_vec(vector<double> al){
        double sum_al = get_sum_vec(al);
        double sum_lgamma_al = 0;
        for(unsigned int i = 0; i < al.size(); i++ ){
             sum_lgamma_al += lgamma(al[i]);
        }
        return -lgamma(sum_al) + sum_lgamma_al;
      };

      /*
       * nd: N x S [n][s]
       * nfg: F x S x FG  [f][s][ fg]
       * z : N x F [n][f]
       * x_by_sf : S x F X N (max)
       */
      
      void update_counting_variables(int* nd, int* nfg, int N, int S, int F, int* z, int* x, int* FG){
        // get histrogram of z, store in nd         
        int maxFG = get_max(FG,F);
        vector< vector< vector<int> > > x_by_sf; // s by f
        vector< vector<int> > ND = convert_flat_to_2d(nd, N, S);
        vector< vector< vector<int> > > NFG = convert_flat_to_3d(nfg, F, S, maxFG);
        vector< vector<int> > Z = convert_flat_to_2d(z, N, F);
        vector< vector<int> > X = convert_flat_to_2d(x, N, F);

        for (int n = 0; n < N; n ++){
            for (int s = 0; s < S; s++){
              ND[n][s] = count(Z[n].begin() , Z[n].end(), s );
              
            }
        };
        convert_2d_to_flat(ND, nd);
        for (int s = 0; s < S; s++){
          vector <vector<int> > x_by_sf_s;
          for (int f = 0; f < F; f++){
            vector<int> x_by_sf_sf;
            for (int n = 0; n < N; n ++){
             if (Z[n][f] == s){
              x_by_sf_sf.push_back( X[n][f] );
             }
            }
            x_by_sf_s.push_back(x_by_sf_sf);
          }
          x_by_sf.push_back( x_by_sf_s );
        }
        for (int f = 0; f < F; f++){
          for (int s = 0; s < S; s++){
            for (int fg = 0; fg < maxFG; fg++){
              if (fg < FG[f]){
                NFG[f][s][fg] = count(x_by_sf[s][f].begin(), x_by_sf[s][f].end(), fg);
              }
              else {
                NFG[f][s][fg] = -1;
              }
            }
          }
        }
        convert_3d_to_flat(NFG, nfg);
        // arrnage x by s and f.
        // update nfg according to x by s and f
      };


      double get_log_multParam4d(double* al,vector<int> NFG_F_S,int len_fg,double q){

        double al_plus_nfg_f_s_sum[len_fg];
        double log_beta_first_term [1];
        double log_beta_second_term [1];

        for(int fg = 0; fg < len_fg; fg++){
          al_plus_nfg_f_s_sum[fg] = al[fg] + NFG_F_S[fg];
        }
        log_Beta(al_plus_nfg_f_s_sum, len_fg, log_beta_first_term );
        log_Beta(al, len_fg, log_beta_second_term );
        double out = log_beta_first_term[0] - log_beta_second_term[0] + log(q);
        return out ;
      }


      void sampleW_subroutine(int N, int S, int F, int*  nd, double alpha, int*  nfg, int* p, int* w, double lamb, int PEAKSIZE, int*  FG, int* x, int* z, double* q, double* wScores_out, int SEED)
      {
        int maxFG = get_max(FG,F);
        vector< vector<int> > ND = convert_flat_to_2d(nd, N, S);
        vector< vector< vector<int> > > NFG = convert_flat_to_3d(nfg, F, S, maxFG);
        vector< vector<int> > P =  convert_flat_to_2d(p, S, F);
        vector< vector<int> > W = convert_flat_to_2d(w, S, F);
        vector< vector<double> > Q = convert_flat_to_2d(q, S, F);
        vector< vector<int> > Z = convert_flat_to_2d(z, N, F);
        vector< vector<int> > X = convert_flat_to_2d(x, N, F);
        vector< vector< vector<double> > > WSCORES;

        vector< double > log_multParam4w;
        vector< double > multParam4w;
        for(int s = 0; s < S; s++){
          vector< vector<double> > WSCORES_f_level;
          for(int f = 0; f < F; f++){
            log_multParam4w.clear();
            multParam4w.clear();
            double al_on[FG[f]];
            double al_off[FG[f]];
            get_al_val(P[s][f], 0, lamb, PEAKSIZE, FG[f], al_off );
            get_al_val(P[s][f], 1, lamb, PEAKSIZE, FG[f], al_on);
            log_multParam4w.push_back( get_log_multParam4d(al_off, NFG[f][s], FG[f], 1 - Q[s][f]));
            log_multParam4w.push_back( get_log_multParam4d(al_on, NFG[f][s], FG[f], Q[s][f]));
            for (unsigned int aa= 0; aa < log_multParam4w.size(); ++aa){
              multParam4w.push_back( exp(log_multParam4w[aa]) );
            }
            multParam4w = normalize_vec( multParam4w );
            if (USE_LOG_SAMPLE) {
            W[s][f] = sample_from_mult_log( multParam4w, SEED);
            }
            else{
            W[s][f] = sample_from_mult( multParam4w, SEED);
            }
            WSCORES_f_level.push_back(multParam4w);
          } // end of F
          WSCORES.push_back(WSCORES_f_level);
        }// end of S
      convert_2d_to_flat(W, w);
      convert_3d_to_flat(WSCORES, wScores_out);
      }


      void sampleW_subroutine_selective(int N, int S, int F, int*  nd, double alpha, int*  nfg, int* p, int* w, double lamb, int PEAKSIZE, int*  FG, int* x, int* z, double* q, double* wScores_out, int SEED, int* subspaces_to_sample, int num_subspaces_to_sample)
      {
        int maxFG = get_max(FG,F);
        vector< vector<int> > ND = convert_flat_to_2d(nd, N, S);
        vector< vector< vector<int> > > NFG = convert_flat_to_3d(nfg, F, S, maxFG);
        vector< vector<int> > P =  convert_flat_to_2d(p, S, F);
        vector< vector<int> > W = convert_flat_to_2d(w, S, F);
        vector< vector<double> > Q = convert_flat_to_2d(q, S, F);
        vector< vector<int> > Z = convert_flat_to_2d(z, N, F);
        vector< vector<int> > X = convert_flat_to_2d(x, N, F);
        vector< vector< vector<double> > > WSCORES;

        vector< double > log_multParam4w;
        vector< double > multParam4w;
        int s = 0;
        for(int ii = 0; ii < num_subspaces_to_sample; ii++){
          s = subspaces_to_sample[ii];
          vector< vector<double> > WSCORES_f_level;
          for(int f = 0; f < F; f++){
            log_multParam4w.clear();
            multParam4w.clear();
            double al_on[FG[f]];
            double al_off[FG[f]];
            get_al_val(P[s][f], 0, lamb, PEAKSIZE, FG[f], al_off );
            get_al_val(P[s][f], 1, lamb, PEAKSIZE, FG[f], al_on);
            log_multParam4w.push_back( get_log_multParam4d(al_off, NFG[f][s], FG[f], 1 - Q[s][f]));
            log_multParam4w.push_back( get_log_multParam4d(al_on, NFG[f][s], FG[f], Q[s][f]));
            for (unsigned int aa= 0; aa < log_multParam4w.size(); ++aa){
              multParam4w.push_back( exp(log_multParam4w[aa]) );
            }
            multParam4w = normalize_vec( multParam4w );

            if (USE_LOG_SAMPLE) {
            W[s][f] = sample_from_mult_log( multParam4w, SEED);
            }
            else{
            W[s][f] = sample_from_mult( multParam4w, SEED);
            }
            WSCORES_f_level.push_back(multParam4w);
          } // end of F
          WSCORES.push_back(WSCORES_f_level);
        }// end of S
      convert_2d_to_flat(W, w);
      convert_3d_to_flat(WSCORES, wScores_out);
      }

      void sampleZ_subroutine(int N, int S, int F, int* nd, double alpha, int* nfg, int*  p ,int*  w, double lamb, int PEAKSIZE, int* FG, int* x, int* z, double* zScore_out, int SEED){

        int maxFG = get_max(FG,F);
        vector< vector<int> > ND = convert_flat_to_2d(nd, N, S);
        vector< vector< vector<int> > > NFG = convert_flat_to_3d(nfg, F, S, maxFG);
        vector< vector<int> > P =  convert_flat_to_2d(p, S, F);
        vector< vector<int> > W = convert_flat_to_2d(w, S, F);
        vector< vector<int> > Z = convert_flat_to_2d(z, N, F);
        vector< vector<int> > X = convert_flat_to_2d(x, N, F);
        vector< vector< vector<double> > > multParam4z_all;

        vector<double> log_multParam4z;
        vector<double> multParam4z;
        //double log_multParam4z_max = -10000;
        //int log_multParam4z_max_idx;
        for(int n=0; n < N; ++n){
          vector< vector<double> > multParam4z_all_n_level;
          for(int f=0; f < F; ++f){
            ND[n][ Z[n][f] ] -= 1;
            NFG[f][ Z[n][f]][ X[n][f] ] -= 1;
            log_multParam4z.clear();
            multParam4z.clear();
            //exp_log_multParam4z_minus_c.clear();
            for(int s=0; s < S; ++s){
              double al[FG[f]];
              get_al_val(P[s][f], W[s][f], lamb, PEAKSIZE, FG[f], al);
              log_multParam4z.push_back( log(alpha/S + ND[n][s] )
                                        + log(NFG[f][s][X[n][f]] + al[ X[n][f] ])
                                        - log( get_sum_vec(NFG[f][s] ) + get_sum(al, FG[f]) )   );
            }
            multParam4z = normalizeLog(log_multParam4z);
            multParam4z_all_n_level.push_back(multParam4z);
            for(int s=0; s < S; ++s){
              //cout << "paramz[" << n << "][" << s << "]" << multParam4z[s] << endl;
            }
            if (USE_LOG_SAMPLE) {
            Z[n][f] = sample_from_mult_log( multParam4z , SEED);
            }
            else{
            Z[n][f] = sample_from_mult( multParam4z , SEED);
            }
            ND[n][Z[n][f]] += 1;
            NFG[f][ Z[n][f]][ X[n][f] ] += 1;
        }
        multParam4z_all.push_back(multParam4z_all_n_level);
      }
    convert_2d_to_flat(Z, z);
    convert_2d_to_flat(ND, nd);
    convert_3d_to_flat(NFG, nfg);
    convert_3d_to_flat(multParam4z_all, zScore_out);
  }

  void sampleP_subroutine(int N, int S, int F, int* nd, double alpha, int* nfg, int*  p ,int* p_index, int*  w, double lamb, int PEAKSIZE, int* FG, int* x, int* z, double* G0, double* pScores_out,int SEED){
    int maxFG = get_max(FG,F);
    vector< vector<int> > ND = convert_flat_to_2d(nd, N, S);
    vector< vector< vector<int> > > NFG = convert_flat_to_3d(nfg, F, S, maxFG);
    vector< vector<int> > P =  convert_flat_to_2d(p, S, F);
    vector< vector<int> > W = convert_flat_to_2d(w, S, F);
    vector< vector<int> > Z = convert_flat_to_2d(z, N, F);
    vector< vector<int> > X = convert_flat_to_2d(x, N, F);
    vector< vector<double> > PSCORES; // = convert_flat_to_2d(pScores_out, S, N);
    // temp buffers
    vector<double> log_multParam4p;
    vector<double> multParam4p;
    vector<double> log_multParam4p_f;
    vector<double> al_plus_nfg_f_s_sum;
    vector<double> al;
    vector<int> p_candidate;
    double logBeta_first_term;
    double logBeta_second_term;

    for(int s=0; s < S; ++s){
      log_multParam4p.clear();
      for(int n=0; n < N; ++n){
        p_candidate.clear();
        p_candidate = X[n];
        log_multParam4p_f.clear();
        for(int f=0; f < F; ++f){
          al.clear();
          al.resize(FG[f]);
          al_plus_nfg_f_s_sum.clear();
          al_plus_nfg_f_s_sum.resize(FG[f]);
          al = getAlVal_vec(p_candidate[f], W[s][f], lamb, PEAKSIZE, FG[f]);
          for(int fg = 0; fg < FG[f]; fg++){
            al_plus_nfg_f_s_sum[fg] = al[fg] + NFG[f][s][fg];
          }
          logBeta_first_term = log_Beta_vec(al_plus_nfg_f_s_sum);
          logBeta_second_term = log_Beta_vec(al);
          log_multParam4p_f.push_back(logBeta_first_term - logBeta_second_term);
        }// end of F
        log_multParam4p.push_back( log( G0[n]) + get_sum_vec(log_multParam4p_f));
      }// end of N
      // normalize using log-sum-exp trick. 
      multParam4p = normalizeLog(log_multParam4p);
      if (USE_LOG_SAMPLE) {
      p_index[s] = sample_from_mult_log( multParam4p , SEED) ;
      }
      else{
      p_index[s] = sample_from_mult( multParam4p , SEED) ;
      }
      P[s] = X[ p_index[s] ];
      PSCORES.push_back(multParam4p);

    }// end of S
    convert_2d_to_flat(P, p);
    convert_2d_to_flat(PSCORES, pScores_out);
  }// end of sampleP_subroutine

  /*
   * Given p, w and z combinations, find the scores.
   * pScores: SxN
   * wScores: SxFx2
   * zScores: NxFxS 
   * simply some parts of sample*_subroutines.
   * */
  void get_scores(int N, int S, int F, int* nd, double alpha, int* nfg, int*  p ,int* p_index, int*  w, double lamb, int PEAKSIZE, int* FG, int* x, int* z, double* G0, double* q, double* zScores_out, double* pScores_out, double* wScores_out){
 /* p sampling
  */ 
    int maxFG = get_max(FG,F);
    vector< vector<int> > ND = convert_flat_to_2d(nd, N, S);
    vector< vector< vector<int> > > NFG = convert_flat_to_3d(nfg, F, S, maxFG);
    vector< vector<int> > P =  convert_flat_to_2d(p, S, F);
    vector< vector<int> > W = convert_flat_to_2d(w, S, F);
    vector< vector<int> > Z = convert_flat_to_2d(z, N, F);
    vector< vector<int> > X = convert_flat_to_2d(x, N, F);
    vector< vector<double> > Q = convert_flat_to_2d(q, S, F);
    vector< vector<double> > PSCORES; // = convert_flat_to_2d(pScores_out, S, N);
    vector< vector< vector<double> > > WSCORES;
    vector< vector< vector<double> > > ZSCORES;

    // temp buffers
    vector<double> log_multParam4p;
    vector<double> multParam4p;
    vector<double> log_multParam4p_f;
    vector<double> al_plus_nfg_f_s_sum;
    vector<double> al;
    vector<int> p_candidate;
    double logBeta_first_term;
    double logBeta_second_term;

    for(int s=0; s < S; ++s){
      log_multParam4p.clear();
      for(int n=0; n < N; ++n){
        p_candidate.clear();
        p_candidate = X[n];
        log_multParam4p_f.clear();
        for(int f=0; f < F; ++f){
          al.clear();
          al.resize(FG[f]);
          al_plus_nfg_f_s_sum.clear();
          al_plus_nfg_f_s_sum.resize(FG[f]);
          al = getAlVal_vec(p_candidate[f], W[s][f], lamb, PEAKSIZE, FG[f]);
          for(int fg = 0; fg < FG[f]; fg++){
            al_plus_nfg_f_s_sum[fg] = al[fg] + NFG[f][s][fg];
          }
          logBeta_first_term = log_Beta_vec(al_plus_nfg_f_s_sum);
          logBeta_second_term = log_Beta_vec(al);
          log_multParam4p_f.push_back(logBeta_first_term - logBeta_second_term);
        }// end of F
        log_multParam4p.push_back( log( G0[n]) + get_sum_vec(log_multParam4p_f));
      }// end of N
      // normalize using log-sum-exp trick. 
      multParam4p = normalizeLog(log_multParam4p);
      PSCORES.push_back(multParam4p);

    }// end of S

  /*
   * 
   */
    vector< double > log_multParam4w;
    vector< double > multParam4w;
    for(int s = 0; s < S; s++){
      vector< vector<double> > WSCORES_f_level;
      for(int f = 0; f < F; f++){
        log_multParam4w.clear();
        multParam4w.clear();
        double al_on[FG[f]];
        double al_off[FG[f]];
        get_al_val(P[s][f], 0, lamb, PEAKSIZE, FG[f], al_off );
        get_al_val(P[s][f], 1, lamb, PEAKSIZE, FG[f], al_on);
        log_multParam4w.push_back( get_log_multParam4d(al_off, NFG[f][s], FG[f], 1 - Q[s][f]));
        log_multParam4w.push_back( get_log_multParam4d(al_on, NFG[f][s], FG[f], Q[s][f]));
        for (unsigned int aa= 0; aa < log_multParam4w.size(); ++aa){
          multParam4w.push_back( exp(log_multParam4w[aa]) );
        }
        multParam4w = normalize_vec( multParam4w );

        WSCORES_f_level.push_back(multParam4w);
      } // end of F
      WSCORES.push_back(WSCORES_f_level);
    }// end of S

  /*
   * z
   */
    vector<double> log_multParam4z;
    vector< vector< vector<double> > > multParam4z_all;
    vector<double> multParam4z;
    //double log_multParam4z_max = -10000;
    //int log_multParam4z_max_idx;
    for(int n=0; n < N; ++n){
      vector< vector<double> > multParam4z_all_n_level;
      for(int f=0; f < F; ++f){
        //ND[n][ Z[n][f] ] -= 1;
        //NFG[f][ Z[n][f]][ X[n][f] ] -= 1;
        log_multParam4z.clear();
        multParam4z.clear();
        for(int s=0; s < S; ++s){
          double al[FG[f]];
          get_al_val(P[s][f], W[s][f], lamb, PEAKSIZE, FG[f], al);
          log_multParam4z.push_back( log(alpha/S + ND[n][s] )
                                    + log(NFG[f][s][X[n][f]] + al[ X[n][f] ])
                                    - log( get_sum_vec(NFG[f][s] ) + get_sum(al, FG[f]) )   );
        }
        multParam4z = normalizeLog(log_multParam4z);
        multParam4z_all_n_level.push_back(multParam4z);
        //ND[n][Z[n][f]] += 1;
        //NFG[f][ Z[n][f]][ X[n][f] ] += 1;
    }
    multParam4z_all.push_back(multParam4z_all_n_level);
  }
  convert_2d_to_flat(PSCORES, pScores_out);
  convert_3d_to_flat(WSCORES, wScores_out);
  convert_3d_to_flat(multParam4z_all, zScores_out);
  }// end of get_scores

} // end of extern C




vector< vector< vector<double>  > > compute_phi(double lamb, int PEAKSIZE, vector< vector<int> >  P,vector< vector<int> > W, int S, int F, int* FG, vector< vector< vector<int> > > NFG)
{
    double sums_al_nfg;

    vector< vector< vector<double> > > PHI;
    for(int s=0; s< S; ++s){
      vector< vector<double> > phi_s;
      for(int f=0; f< F; ++f){
        vector<double> phi_s_f;
        double al[FG[f]];
        get_al_val(P[s][f], W[s][f], lamb, PEAKSIZE, FG[f], al);
        sums_al_nfg = get_sum(al, FG[f]) + get_sum_vec(NFG[f][s]);
        for(int fg=0; fg< FG[f]; ++fg){
          phi_s_f.push_back( (NFG[f][s][fg] + al[fg] )/sums_al_nfg);
        }
        phi_s.push_back(phi_s_f);
      }
      PHI.push_back(phi_s);
    }
    return PHI;
}

vector< vector<double> >  compute_pi(int N, int S, vector< vector<int> > Z)
{
    vector< vector<double> >  PI;
    vector<double>  counts;
    for(int n=0; n < N; ++n){
      counts.clear();
      for(int s=0; s< S; ++s){
      counts.push_back( count(Z[n].begin() , Z[n].end(), s ) );
      }
      PI.push_back( normalize_vec(counts));
    }
    return PI;
}


extern "C" {
  void set_random_seed(int SEED){
    if (USE_SRAND){
      cout << "setting random seed from bcm.cpp" << endl;
      srand(SEED);
    }
  }
  double calculate_likelihood(int N, int S, int F, int* nfg, int*  p ,int*  w, double lamb, int PEAKSIZE, int* FG, double * q,  int* z, int *x){
    
    int maxFG = get_max(FG,F);
    vector< vector< vector<int> > > NFG = convert_flat_to_3d(nfg, F, S, maxFG);
    vector< vector<int> > P =  convert_flat_to_2d(p, S, F);
    vector< vector<int> > W = convert_flat_to_2d(w, S, F);
    vector< vector<double> > Q = convert_flat_to_2d(q, S, F);

    vector< vector<int> > Z = convert_flat_to_2d(z, N, F);
    vector< vector<int> > X = convert_flat_to_2d(x, N, F);

    vector< vector< vector<double> > > PHI = compute_phi(lamb, PEAKSIZE, P, W, S, F, FG, NFG);
    vector < vector<double> >  PI = compute_pi(N, S, Z);

    //vector<double> al_plus_nfg_f_s_sum;
    //vector<double> al;
    //double logBeta_first_term;
    //double logBeta_second_term;
    double logLikelihood = 0;

    for(int s=0; s < S; ++s){
      for(int f=0; f < F; ++f){
        logLikelihood += -log(N) + log(Q[s][f]);
      }
    }
    for(int n=0; n < N; ++n){
      for(int f=0; f < F; ++f){
        logLikelihood += log(PHI[ Z[n][f] ][f][X[n][f]])  + PI[n][Z[n][f]];
      }
    }
    return logLikelihood;
  }

  double sample_all(int N, int S, int F, int* nd, double alpha, int* nfg, int*  p ,int* p_index, int*  w, double lamb, int PEAKSIZE, int* FG, int* x, int* z, double* G0, double* q, double* zScores_out, double* pScores_out, double* wScores_out, int SEED){
      sampleP_subroutine(N, S, F, nd, alpha, nfg, p ,p_index, w, lamb, PEAKSIZE, FG, x, z, G0, pScores_out, SEED);
      sampleW_subroutine(N, S, F, nd, alpha,  nfg, p, w, lamb, PEAKSIZE, FG, x, z, q, wScores_out, SEED);
      sampleZ_subroutine(N, S, F, nd, alpha, nfg, p ,w, lamb, PEAKSIZE, FG, x, z, zScores_out, SEED);
      return calculate_likelihood(N, S, F, nfg, p , w, lamb, PEAKSIZE, FG, q, z, x);
  }

  double sample_all_iterations(int N, int S, int F, int* nd, double alpha, int* nfg, int*  p ,int* p_index, int*  w, double lamb, int PEAKSIZE, int* FG, int* x, int* z, double* G0, double* q, int ITERATION, double* zScores_out, double* pScores_out, double* wScores_out, int SEED){
    for(int iter = 0; iter < ITERATION ; ++iter){
      sampleP_subroutine(N, S, F, nd, alpha, nfg, p ,p_index, w, lamb, PEAKSIZE, FG, x, z, G0, pScores_out, SEED);
      sampleW_subroutine(N, S, F, nd, alpha,  nfg, p, w, lamb, PEAKSIZE, FG, x, z, q, wScores_out, SEED);
      sampleZ_subroutine(N, S, F, nd, alpha, nfg, p ,w, lamb, PEAKSIZE, FG, x, z, zScores_out, SEED);
    }
    return calculate_likelihood(N, S, F, nfg, p , w, lamb, PEAKSIZE, FG, q, z, x);
  }

} // end of extern C




















