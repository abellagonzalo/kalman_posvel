#include "classic_matrices.h"

namespace kalman {

  class Kalman {
  public:
    Kalman() {
      init();
    }

    Kalman(float* x_init, float* P_init)
    : x(4, 1, x_init)
    , P(4, 4, P_init)
    { 
      init(); 
    }

    void reset(float* x_init, float* P_init) {
      x = MatrixCM(4, 1, x_init);
      P = MatrixCM(4, 4, P_init);
    }

    void predict(float* u_init, float* w_init) {
      using namespace matrices;
      MatrixCM u(6, 1, u_init);
      MatrixCM w(6, 1, w_init);
      x = f(x, u, zero6);
      P = A(x, u, w) * P * A(x, u, w).transpose() + W(x, u, w) * Q(w) * W(x, u, w).transpose();      
    }

    void update(float* z_init, float* v_init) {
      using namespace matrices;
      MatrixCM z(2, 1, z_init);
      MatrixCM v(2, 1, v_init);
      MatrixCM aux = H(x, v) * P * H(x, v).transpose() + V(x, v) * R(v) * V(x, v).transpose();
      MatrixCM gain = P * H(x, v).transpose() * aux.inverse();
      x = x + gain * ( z - h(x, zero4) );
      P = ( iden4 - gain * H(x,v) ) * P;      
    }

    MatrixCM x;
    MatrixCM P;

  private:

    void init() {
      float zero4_init [4] = { 0., 0., 0., 0.};
      float zero6_init [6] = { 0., 0., 0., 0., 0., 0.};
      zero4 = MatrixCM(4, 1, zero4_init);
      zero6 = MatrixCM(6, 1, zero6_init);
      iden4 = MatrixCM(4);
    }

    MatrixCM zero4;
    MatrixCM zero6;
    MatrixCM iden4;
  };

}

int main(int argc, char **argv) {

  using namespace kalman::matrices;

  float x_init [4] = { 100, 0, 0, 0 };
  float P_init [16] = { 100,   0,    0,    0,
                          0, 100,    0,    0,
                          0,   0, 1000,    0,
                          0,   0,    0, 1000 };

  kalman::Kalman kalman(x_init, P_init);

  int imax = 100;
  for (int i=0; i<imax; ++i) {
    float u [6] = { 0, 0, 0, 30, 0, 0.9 };
    float w [6] = { 0.01, 0.01, 0.01, 0, 10, 0 };
    kalman.predict(u, w);

    float z [2] = { 100 + 100*i, 0 };
    float v [2] = { .1, .1 };
    kalman.update(z, v);
  }


  std::cerr << "*** x ***";
  kalman.x.print();
  std::cerr << "*** P ***";
  kalman.P.print();

  return 0;
}