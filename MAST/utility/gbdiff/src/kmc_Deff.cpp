#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib> 
#include <cmath>  
#include <cstring>

// input parameters that can be changed by users

// ---------------------------------------------

const double TEMP = 1373.0;     // unit K
const double ACTENGL = 2.70;    // unit eV
const double ACTENGH = 2.04;    // unit eV
const double DEFF0L = 1.0e-9;   // unit m^2/s
const double DEFF0H = 1.0e-9;   // unit m^2/s
const double FRACTION = 0.0;    // between 0 and 1
const int NLOOP = 10000;        // default value

// ---------------------------------------------

const int LOWDIFF = 0;
const int HIGHDIFF = 1;
const int NATOMS = 1;
const int NGRAINS = 4000;
const int NT = 100000000;    
const int NT2 = 100;       
const int NT3 = NT/NT2;     
const int NOUT = 50;
const double KB = 8.617e-5;
const double FACTOR = 2.15;
const double XMIN = 0.0;
const double XMAX = 1.0;
const double SMALL = 1e-6;
const double ALPHA = 2.0;
const double DELTA = 1.0e-4;
  
using namespace std;

double temperature(const double & time);

class GBnetwork {

private:

  int ngb_;
  int nxmin_, nxmax_;
  int ncs_;                             
  vector<double> acteng_;
  vector<double> deff0_;
  vector<int> gbtype_;                  
  vector< vector<double> > node_pos_;   
  vector< vector<int> > elem_;          
  vector< vector<int> > gb_elem_;       
  vector< vector<int> > surf_gb_;       
  vector< vector<int> > surf_node_;     
  vector< vector<int> > surf_node_gb_;  
  vector< vector<int> > node_con_;      
  vector< vector<double> > node_rate_ref_; 
  vector< vector< vector<double> > > node_diff_len_; 
  vector< vector<int> > node_numstep_;  
  vector<int> iatoms_;
  vector<double> events_;                
  vector<double> eventsa_;              
  vector<double> timea_;
  vector<double> msda_, msdxa_, msdya_, msdza_;
  vector<double> Deff_;
  double aveD_;
  vector< vector<double> > disD_;
  double stddev_;

public:

  GBnetwork();
  void input();
  void getSurfNodes();
  double calDiff(const int n1, const int n2, const int gbtype, 
                 const int boundary, vector<double> & length,  
                 const double temp, const double dhdl = -1.0) const;
  void createRateList(const double & fh, const double & dhdl = - 1.0);
  void oneStep(vector<double> & mdis, double & dt, 
               const double temp, const double dhdl);
  void allSteps(const double & dhdl = - 1.0);
  void statAna();
  double aveDeff() { return aveD_; }
  void outputData1(const double temp, const double fh, const double dhdl = -1.0) const;
};


class KMC : public GBnetwork {

  public:

    void performKMC();
    void outputData2(const double temp, const double dhdl, 
                     const vector<double> & fh,
                     const vector<double> & aveD) const;
};


int main()
{
  KMC s;

  s.performKMC();
 
  return 0;
}


GBnetwork::GBnetwork() {

  acteng_.push_back(ACTENGL);
  acteng_.push_back(ACTENGH);
  deff0_.push_back(DEFF0L);
  deff0_.push_back(DEFF0H);

  iatoms_.resize(NATOMS);

  ncs_ = (int) ((ALPHA - 1.0) * log(1.0/DELTA) / DELTA + 1.0);

  srand(time(NULL));
}


void
GBnetwork::input() {

  char tessfile[100];
  char meshfile[100];
  char ng[100];
  int i, j, k;
  int igb, nnodes, count;
  int k1, k2, k3;
  double x, y, z;
  vector<int> gb_list;
  vector<int> one_elem;
  vector<double> pos;
  vector<int> newindex;
  vector< vector<double> > node_pos_ini;

  strcpy(tessfile, "n");
  sprintf(ng, "%d", NGRAINS);
  strcat(tessfile, ng);
  strcat(tessfile, ".tess");

  strcpy(meshfile, "n");
  sprintf(ng, "%d", NGRAINS);
  strcat(meshfile, ng);
  strcat(meshfile, ".msh");

  ifstream input;

  input.open(tessfile);

  if(!input) { 
    cerr << "Error: The input file can not be opened!" << endl;
    exit(1);
  }

  input >> ngb_;

  surf_gb_.clear();

  for (i = 0; i < 6; ++i) {  

    input >> k;

    gb_list.clear();

    for (j = 0; j < k; ++j) {

      input >> igb;
      gb_list.push_back(igb-1);
    }

    surf_gb_.push_back(gb_list);
  }

  input.close();

  nxmin_ = nxmax_ = 0;

  input.open(meshfile);

  if(!input) { 
    cerr << "Error: The input file can not be opened!" << endl;
    exit(1);
  }

  input >> nnodes;

  node_pos_ini.clear();

  for (i = 0; i < nnodes; ++i) {

    input >> j >> x >> y >> z;

    pos.clear();  
    pos.push_back(x); pos.push_back(y); pos.push_back(z);

    node_pos_ini.push_back(pos);

    if (fabs(x - XMIN) < SMALL) ++nxmin_;
    if (fabs(x - XMAX) < SMALL) ++nxmax_;
  }

  elem_.clear();
 
  gb_elem_.resize(ngb_);

  count = -1;

  while (true) {

    input >> k;

    if (input.eof()) break;

    ++count;

    input >> k >> k >> k >> igb >> k >> i >> j >> k;

    one_elem.clear();

    one_elem.push_back(igb-1); one_elem.push_back(i-1);
    one_elem.push_back(j-1);   one_elem.push_back(k-1); 

    elem_.push_back(one_elem);

    gb_elem_[igb-1].push_back(count);

  }  

  input.close(); 

  k1 = k2 = k3 = -1;
 
  newindex.resize(nnodes);

  for (i = 0; i < nnodes; ++i) {

    if (fabs(node_pos_ini[i][0] - XMIN) < SMALL) {

      ++k1;
      newindex[i] = k1;
    }
    else if (fabs(node_pos_ini[i][0] - XMAX) < SMALL) {

      ++k2;
      newindex[i] = nnodes - nxmax_ + k2;
    }
    else {

      ++k3;
      newindex[i] = nxmin_ + k3;
    }
  }

  node_pos_.clear();  node_pos_.resize(nnodes);

  for (i = 0; i < nnodes; ++i) {

    node_pos_[newindex[i]] = node_pos_ini[i];

  }

  for (i = 0; i < elem_.size(); ++i) {

    for (j = 1; j <=3; ++j) {

      elem_[i][j] = newindex[elem_[i][j]];

    }
  }

  for (i = 0; i < 6; ++i) {

    for (j = 0; j < surf_gb_[i].size(); ++j) {

      for (k = 0; k < gb_elem_[surf_gb_[i][j]].size(); ++k) {

        if (elem_[gb_elem_[surf_gb_[i][j]][k]][0] != surf_gb_[i][j]) {

          cerr << "Error: data input is wrong!" << endl;
          exit(1);
        }
      }
    }
  }

  getSurfNodes();

}


void
GBnetwork::getSurfNodes() {

  int i, j, k, m, n, p, in;
  vector<int> node_list;
  vector<int> gb_list;

  surf_node_.clear();
  surf_node_gb_.clear();

  for (i = 0; i < 6; ++i) {

    node_list.clear();
    gb_list.clear();

    for (j = 0; j < surf_gb_[i].size(); ++j) {    

      for (k = 0; k < gb_elem_[surf_gb_[i][j]].size(); ++k) {

        for (m = 1; m <=3; ++m) {

          n = elem_[gb_elem_[surf_gb_[i][j]][k]][m];

          in = 1;

          for (p = 0; p < node_list.size(); ++p) {

            if (n == node_list[p]) { in = 0; break; }
          }

          if (in == 1) {

            node_list.push_back(n);
            gb_list.push_back(surf_gb_[i][j]);

          }
        }   
      }
    }

    surf_node_.push_back(node_list);
    surf_node_gb_.push_back(gb_list);
  }

}


double
GBnetwork::calDiff(const int n1, const int n2, const int gbtype, 
                   const int boundary, vector<double> & length,  
                   const double temp, const double dhdl) const {

  double rate, dis2, D;
  double maxlen2, minlen2, medlen2;
  double dx, dy, dz;

  maxlen2 = 1.0 / pow(1.0 * NGRAINS, 2.0/3.0);  
  minlen2 = maxlen2 * 0.01;
  medlen2 = maxlen2 * 0.16;

  length.clear();

  if (boundary == 0) { 

    dis2 = (node_pos_[n1][0] - node_pos_[n2][0]) * (node_pos_[n1][0] - node_pos_[n2][0]) +
           (node_pos_[n1][1] - node_pos_[n2][1]) * (node_pos_[n1][1] - node_pos_[n2][1]) +
           (node_pos_[n1][2] - node_pos_[n2][2]) * (node_pos_[n1][2] - node_pos_[n2][2]);

    if (dis2 < minlen2) {

      dx = (node_pos_[n2][0] - node_pos_[n1][0]) * minlen2 / dis2;
      dy = (node_pos_[n2][1] - node_pos_[n1][1]) * minlen2 / dis2;
      dz = (node_pos_[n2][2] - node_pos_[n1][2]) * minlen2 / dis2;

      dis2 = minlen2;
    }
    else if (dis2 > maxlen2) {

      dx = (node_pos_[n2][0] - node_pos_[n1][0]) * maxlen2 / dis2;
      dy = (node_pos_[n2][1] - node_pos_[n1][1]) * maxlen2 / dis2;
      dz = (node_pos_[n2][2] - node_pos_[n1][2]) * maxlen2 / dis2;

      dis2 = maxlen2;
    }
    else {

      dx = node_pos_[n2][0] - node_pos_[n1][0];
      dy = node_pos_[n2][1] - node_pos_[n1][1];
      dz = node_pos_[n2][2] - node_pos_[n1][2];
    }

    length.push_back(dx); length.push_back(dy); length.push_back(dz); 
  }
  else {

    dis2 = (node_pos_[n1][0] - node_pos_[n2][0]) * (node_pos_[n1][0] - node_pos_[n2][0]) +
           (node_pos_[n1][1] - node_pos_[n2][1]) * (node_pos_[n1][1] - node_pos_[n2][1]) +
           (node_pos_[n1][2] - node_pos_[n2][2]) * (node_pos_[n1][2] - node_pos_[n2][2]);

    dx = (node_pos_[n1][0] - node_pos_[n2][0]) * medlen2 / dis2;
    dy = (node_pos_[n1][1] - node_pos_[n2][1]) * medlen2 / dis2;
    dz = (node_pos_[n1][2] - node_pos_[n2][2]) * medlen2 / dis2;

    dis2 = medlen2;

    length.push_back(dx); length.push_back(dy); length.push_back(dz);
  }

  if (dhdl < 0.0) { 

    D = deff0_[gbtype] * exp( -acteng_[gbtype] / (KB * temp));
  }
  else { 

    if (gbtype == LOWDIFF) {

      D = deff0_[gbtype] * exp( -acteng_[gbtype] / (KB * temp));
    }
    else {

      D = dhdl * deff0_[LOWDIFF] * exp( -acteng_[LOWDIFF] / (KB * temp));
    }
  }  

  rate = D / dis2;

  return (rate);
}


void
GBnetwork::createRateList(const double & fh, const double & dhdl) {

  int i, j, k, m;
  int index, n1, n2;
  double diff;
  vector<double> length;
  vector<int> steps;

  node_con_.clear();        node_con_.resize(node_pos_.size());
  node_rate_ref_.clear();   node_rate_ref_.resize(node_pos_.size());
  node_diff_len_.clear();   node_diff_len_.resize(node_pos_.size());

  gbtype_.resize(ngb_);

  for (i = 0; i < gbtype_.size(); ++i) {

    if (rand() / (double) RAND_MAX < fh) { gbtype_[i] = HIGHDIFF; }
    else { gbtype_[i] = LOWDIFF; }
  }

  for (i = 0; i < elem_.size(); ++i) {

    for (n1 = 1; n1 < 3; ++n1) {
      for (n2 = n1+1; n2 <= 3; ++n2) {

        index = 1; 

        for (j = 0; j < node_con_[elem_[i][n1]].size(); ++j) {

          if (elem_[i][n2] == node_con_[elem_[i][n1]][j]) {index = 0; break;}
        }

        if (index == 1) {

          diff = calDiff(elem_[i][n1], elem_[i][n2], gbtype_[elem_[i][0]],
                         0, length, temperature(0.0), dhdl);

          node_con_[elem_[i][n1]].push_back(elem_[i][n2]);
          node_con_[elem_[i][n2]].push_back(elem_[i][n1]);

          node_rate_ref_[elem_[i][n1]].push_back(diff);
          node_rate_ref_[elem_[i][n2]].push_back(diff);

          node_diff_len_[elem_[i][n1]].push_back(length);
          for (int k = 0; k < length.size(); ++k) { length[k] *= -1.0; }
          node_diff_len_[elem_[i][n2]].push_back(length);
        }
      }  
    }  
  }  

  node_numstep_.clear();

  for (i = 0; i < node_con_.size(); ++i) {

    steps.clear();

    steps.resize(node_con_[i].size(), 0);

    node_numstep_.push_back(steps);
  }

}


void
GBnetwork::oneStep(vector<double> & mdis, double & dt, 
                   const double temp, const double dhdl) {

  int i, k;
  int oldatom, newatom;
  double rate, totrates, rr;

  events_.clear();
  eventsa_.clear();

  events_ = node_rate_ref_[iatoms_[0]];

  for (i = 0; i < events_.size(); ++i) {

    if (i == 0) { eventsa_.push_back(events_[i]); }
    else { eventsa_.push_back(eventsa_[i-1] + events_[i]); }
  }

  totrates = eventsa_[eventsa_.size()-1];

  rr = rand()/ (double) RAND_MAX;

  if ( rr > 0.0) {dt = -log(rr) / totrates;}
  else { dt = 0.0; }

  rate = rand() / (double) RAND_MAX * totrates;

  for (i = 0; i < eventsa_.size(); ++i) {

    if (rate <= eventsa_[i]) {k = i; break;}
  }

  mdis.clear();

  mdis = node_diff_len_[iatoms_[0]][k];

  oldatom = iatoms_[0];

  newatom = node_con_[iatoms_[0]][k];

  node_numstep_[oldatom][k] = node_numstep_[oldatom][k] + 1;

  if (node_numstep_[oldatom][k] > ncs_) {

    node_numstep_[oldatom][k] = 0;

    node_rate_ref_[oldatom][k] = node_rate_ref_[oldatom][k] / ALPHA;

    for (i = 0; i < node_con_[newatom].size(); ++i) {

      if (node_con_[newatom][i] == oldatom) {

        node_numstep_[newatom][i] = 0;

        node_rate_ref_[newatom][i] = node_rate_ref_[newatom][i] / ALPHA;

        break;
      }
    }
  }

  iatoms_[0] = node_con_[iatoms_[0]][k];
}


void
GBnetwork::allSteps(const double & dhdl) {

  int k, kk;
  double tottime, totdx, totdy, totdz;
  double temp, dt;
  vector<double> mdis;

  timea_.resize(NT2, 0.0);  msda_.resize(NT2, 0.0);
  msdxa_.resize(NT2, 0.0);  msdya_.resize(NT2, 0.0);  msdza_.resize(NT2, 0.0);
  Deff_.resize(NLOOP, 0.0);

  for (kk = 0; kk < NLOOP; ++kk) {

    tottime = totdx = totdy = totdz = 0.0;

    temp = temperature(0.0);

    iatoms_[0] = (int) (0.999999 * rand() * nxmin_ / RAND_MAX);

    for (k = 0; k < NT; ++k) {

      oneStep(mdis, dt, temp, dhdl);

      tottime += dt;
      totdx += mdis[0];  totdy += mdis[1];  totdz += mdis[2];

      temp = temperature(tottime);

      if ((k+1)%NT3 == 0) {

        msdxa_[(k+1)/NT3-1] += totdx * totdx / (1.0 * NLOOP);
        msdya_[(k+1)/NT3-1] += totdy * totdy / (1.0 * NLOOP);
        msdza_[(k+1)/NT3-1] += totdz * totdz / (1.0 * NLOOP);

        msda_[(k+1)/NT3-1] += (totdx * totdx + totdy * totdy + totdz * totdz) / (1.0 * NLOOP);

        timea_[(k+1)/NT3-1] += tottime / (1.0 * NLOOP);
      }

      if (iatoms_[0] >= node_pos_.size() - nxmax_) break;
    }

    Deff_[kk] = (XMAX - XMIN) * (XMAX - XMIN) / (2.0 * tottime * FACTOR);
  }
}


void
GBnetwork::statAna() {

  int i, j;
  double maxD, minD, invD;

  aveD_ = 0.0; maxD = -1.0; minD = 1000000.0;

  for (i = 0; i < NLOOP; ++i) {

    aveD_ += Deff_[i];
  
    if (maxD < Deff_[i]) maxD = Deff_[i];
    if (minD > Deff_[i]) minD = Deff_[i];
  }

  aveD_ /= 1.0 * NLOOP;

  maxD = maxD * 0.8;

  invD = (maxD - minD) / (1.0 * NOUT);

  disD_.clear();   disD_.resize(NOUT);

  for (i = 0; i < NOUT; ++i) {

    disD_[i].push_back(minD + invD * (i+0.5));
    disD_[i].push_back(0.0);
  }

  for (i = 0; i < NOUT; ++i) {
    for (j = 0; j < NLOOP; ++j) {

      if (Deff_[j] >= (minD + invD * i) && Deff_[j] <= (minD + invD * (i+1)) ) {
         disD_[i][1] += 1.0;
      }
    }

    disD_[i][1] /= 1.0 * NLOOP;
  }

  stddev_ = 0.0;

  for (i = 0; i < NLOOP; ++i)
    stddev_ += (aveD_ - Deff_[i]) * (aveD_ - Deff_[i]) / (1.0 * NLOOP);

  stddev_ = sqrt(stddev_);
}

void
GBnetwork::outputData1(const double temp, const double fh, const double dhdl) const {

  int i, j;

  ofstream output;

  output.open("Deff.txt");

  output << aveD_ << endl;

  output.close();

  output.open("Deff_dis.txt");

  for (i = 0; i < NOUT; ++i) {
    for (j = 0; j < disD_[i].size(); ++j) {
      output << disD_[i][j] << " ";
    }
    output << endl;
  }

  output.close();

  output.open("Deff_dev.txt");

  output << stddev_ << endl;

  output.close();

}


void
KMC::performKMC() {

  int i, j;
  double temp;
  vector<double> dhdl, fh;
  vector<double> aveD;

  temp = temperature(0.0);

  fh.push_back(FRACTION);

  GBnetwork a;

  a.input();

  aveD.clear();

  for (j = 0; j < fh.size(); ++j) {

    a.createRateList(fh[j]);

    a.allSteps();

    a.statAna();

    aveD.push_back(a.aveDeff());

    a.outputData1(temp, fh[j]);
  }
}

double temperature(const double & time) {

  return (TEMP);
}