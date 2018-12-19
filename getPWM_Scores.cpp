#include <Rcpp.h>
using namespace Rcpp;

#include <stdlib.h>
#include <stdio.h>

#include <map>
#include <vector>

using namespace std;

double initMap(NumericMatrix m, map<char, vector<double> > &mpw);
double getScoreAtPos(string &s, int pos, map<char, vector<double> > &mpw, int lg, double minScore);
int printMap(map<char, vector<double> >& m);


// [[Rcpp::export]]
void test(NumericMatrix m, string s){

  map<char, vector<double> > mpw;
  
  //////////////////////////////////////////////////////////////
  // map initialization
  double minScore = initMap(m, mpw);

  int lg = mpw['A'].size();

  for(int pos=0 ; pos<s.size() ; pos++){
    double score = getScoreAtPos(s, pos, mpw, lg, minScore);
    fprintf(stdout, "%d -> %.4lf\n", pos, score);
  }
  
  //printMap(mpw);
  
  //return 0;
}

// [[Rcpp::export]]
double getScoreAtPos(NumericMatrix m, string &s, int pos){
  map<char, vector<double> > mpw;
  double minScore = initMap(m, mpw);
  int lg = mpw['A'].size();
  double score = getScoreAtPos(s, pos, mpw, lg, minScore);
  return score;
}


// [[Rcpp::export]]
NumericVector getScores(NumericMatrix m, string &s){

  map<char , vector<double> > mpw;
  double minScore = initMap(m, mpw);

  //printMap(mpw);

  int lgMap = mpw['A'].size();
  
  int slg = s.size();
  NumericVector vres(slg);
  
  for(int pos=0 ; pos<slg ; pos++){
    double score = getScoreAtPos(s, pos, mpw, lgMap, minScore);
    vres[pos] = score;
  }
  return vres;
}

double getScoreAtPos(string &s, int pos, map<char, vector<double> > &mpw, int lg, double minScore){
  double score = 0.0;
  int slg = s.size();
  int i;
  int j;
  for(i=pos, j=0 ; i<slg && j<lg ; i++,j++){
    char base = toupper(s[i]);
    if( base=='A' || base=='T' || base=='C' || base=='G' ){
      double baseScore = mpw[base][j];
      //fprintf(stderr, "base:%c (->%.4lf)\n", base, baseScore);
      score += baseScore;
    } else {
      score += minScore;
    }
  }
  return score;
}

int printMap(map< char , vector<double> > &m){
  
  map< char , vector<double> >::iterator itm;
  vector<char> vbases;
  
  for( itm = m.begin() ; itm!=m.end() ; itm++ ){
    vbases.push_back(itm->first);
  }
  int nbBases = vbases.size();

  Rcpp::Rcout<<"pos";
  for(int j=0 ; j<nbBases ; j++){
    Rcpp::Rcout<<"\t"<<vbases[j];
  }
  Rcpp::Rcout<<"\n";
  
  int lg = m[vbases[0]].size();
  for(int i=0 ; i<lg ; i++){
    Rcpp::Rcout<<i;
    for(int j=0 ; j<nbBases ; j++){
      char base = vbases[j];
      double val = m[base][i];
      Rcpp::Rcout<<"\t"<<val;
    }
    Rcpp::Rcout<<"\n";
  }

  return 0;
}

double initMap(NumericMatrix m, map<char , vector<double> > &mpw){
  CharacterVector cnames = colnames(m);
  int nbc = cnames.size();
  double minScore = 0.0;
  
  vector<char> vbases;

  for(int ibase=0 ; ibase<nbc ; ibase++){

    const char *basep = CHAR(cnames(ibase));
    char base = (*basep);

    NumericVector nv = m( _, ibase);
    vector<double> vtmp = as< vector<double> >(nv);
    //mpw[base] = vtmp;
    int lg = vtmp.size();
    for(int i=0 ; i<lg ; i++){
      //fprintf(stdout, "%.2f\n", vtmp[i]);
      double score = vtmp[i];
      mpw[base].push_back(score);
      if(score<minScore){
	minScore = score;
      }
    }
  }
  return minScore;
}
