#include "kalman.h"
#include <iostream>
#include <vector>

void print_pair(kalman::algorithm::MatrixCMPair pair) {
  pair.first.print();
  pair.second.print();
}

int main(int argc, char **argv) {

  std::vector<kalman::algorithm::MatrixCMPair> pairs(2, kalman::algorithm::pair_generator());

  int imax = 100;
  for (int i=0; i<imax; ++i) {
    double u [6] = { 0, 0, 0, 1./30., 0, 0.9 };
    double w [6] = { 0.01, 0.01, 0.01, 0, 10, 0 };
    kalman::algorithm::predict_collection(pairs.begin(), pairs.end(), MatrixCM(6,1,u), MatrixCM(6,1,w));

    double z [2] = { 100 + 100*i, 0 };
    double v [2] = { 30., 30. };
    kalman::algorithm::update_collection(pairs.begin(), pairs.end(), MatrixCM(2,1,z), MatrixCM(2,1,v));

    double z2 [2] = { 0, 100 + 100*i };
    double v2 [2] = { 30., 30. };
    kalman::algorithm::update_collection(pairs.begin(), pairs.end(), MatrixCM(2,1,z2), MatrixCM(2,1,v2));
  }

  std::for_each(pairs.begin(), pairs.end(), print_pair);

  return 0;
}
