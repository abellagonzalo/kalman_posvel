#include "kalman.h"
#include <iostream>

int main(int argc, char **argv) {

  float x_init [4] = { 100, 0, 0, 0 };
  float P_init [16] = { 100,   0,    0,    0,
                          0, 100,    0,    0,
                          0,   0, 1000,    0,
                          0,   0,    0, 1000 };

/*
  kalman::Kalman mykalman(x_init, P_init);

  int imax = 100;
  for (int i=0; i<imax; ++i) {
    float u [6] = { 0, 0, 0, 30, 0, 0.9 };
    float w [6] = { 0.01, 0.01, 0.01, 0, 10, 0 };
    mykalman.predict(u, w);

    float z [2] = { 100 + 100*i, 0 };
    float v [2] = { .1, .1 };
    mykalman.update(z, v);
  }


  std::cerr << "*** x ***";
  mykalman.x.print();
  std::cerr << "*** P ***";
  mykalman.P.print();
*/
  return 0;
}