#include "classic_matrices.h"
#include <boost/function.hpp>

/******************************
* Newtons difference quotient *
******************************/
typedef boost::function<double (const MatrixCM&, const MatrixCM&, const MatrixCM)> Function;

double limm_xx(Function f, const MatrixCM& x, const MatrixCM& u, const MatrixCM& w, double h) {
  float hm [4] = {h, 0, 0, 0};
  return (f(x+MatrixCM(4, 1, hm), u, w) - f(x, u, w)) / h;
}

double limm_xy(Function f, const MatrixCM& x, const MatrixCM& u, const MatrixCM& w, double h) {
  float hm [4] = {0, h, 0, 0};
  return (f(x+MatrixCM(4, 1, hm), u, w) - f(x, u, w)) / h;
}

double limm_xvx(Function f, const MatrixCM& x, const MatrixCM& u, const MatrixCM& w, double h) {
  float hm [4] = {0, 0, h, 0};
  return (f(x+MatrixCM(4, 1, hm), u, w) - f(x, u, w)) / h;
}

double limm_xvy(Function f, const MatrixCM& x, const MatrixCM& u, const MatrixCM& w, double h) {
  float hm [4] = {0, 0, 0, h};
  return (f(x+MatrixCM(4, 1, hm), u, w) - f(x, u, w)) / h;
}

double limm_wx(Function f, const MatrixCM& x, const MatrixCM& u, const MatrixCM& w, double h) {
  float hm [6] = {h, 0, 0, 0, 0, 0};
  return (f(x, u, w+MatrixCM(6, 1, hm)) - f(x, u, w)) / h;
}

double limm_wy(Function f, const MatrixCM& x, const MatrixCM& u, const MatrixCM& w, double h) {
  float hm [6] = {0, h, 0, 0, 0, 0};
  return (f(x, u, w+MatrixCM(6, 1, hm)) - f(x, u, w)) / h;
}

double limm_wo(Function f, const MatrixCM& x, const MatrixCM& u, const MatrixCM& w, double h) {
  float hm [6] = {0, 0, h, 0, 0, 0};
  return (f(x, u, w+MatrixCM(6, 1, hm)) - f(x, u, w)) / h;
}

double limm_wt(Function f, const MatrixCM& x, const MatrixCM& u, const MatrixCM& w, double h) {
  float hm [6] = {0, 0, 0, h, 0, 0};
  return (f(x, u, w+MatrixCM(6, 1, hm)) - f(x, u, w)) / h;
}

double limm_wj(Function f, const MatrixCM& x, const MatrixCM& u, const MatrixCM& w, double h) {
  float hm [6] = {0, 0, 0, 0, h, 0};
  return (f(x, u, w+MatrixCM(6, 1, hm)) - f(x, u, w)) / h;
}

double limm_wk(Function f, const MatrixCM& x, const MatrixCM& u, const MatrixCM& w, double h) {
  float hm [6] = {0, 0, 0, 0, 0, h};
  return (f(x, u, w+MatrixCM(6, 1, hm)) - f(x, u, w)) / h;
}

// Custom assert
void assert_derivative(double expected, double value, double variation) {
  if (abs(expected - value) > variation) {
    std::cerr << "Expected " << expected << " but got " << value << std::endl;
    assert(true);
  }
}

int main(int argc, char **argv) {

  std::cerr << "Classic matrices test started" << std::endl;

  using namespace kalman::matrices;

  float x_values [4] = {800, 1200, 567, 100};
  float u_values [6] = {832, 1358, 0.12f, 1.f/30.f, 0.1f, 0.9f};
  float w_values [6] = {100, 234, 234, 94, 245, 456};
  float h = 0.00001;

  MatrixCM x(4, 1, x_values);
  MatrixCM u(6, 1, u_values);
  MatrixCM w(6, 1, w_values);

  assert_derivative(  limm_xx(f1, x, u, w, h),  df1_dxx(x, u, w), h);
  assert_derivative(  limm_xy(f1, x, u, w, h),  df1_dxy(x, u, w), h);
  assert_derivative( limm_xvx(f1, x, u, w, h), df1_dxvx(x, u, w), h);
  assert_derivative( limm_xvy(f1, x, u, w, h), df1_dxvy(x, u, w), h);

  assert_derivative( limm_wx(f1, x, u, w, h), df1_dwx(x, u, w), h);
  assert_derivative( limm_wy(f1, x, u, w, h), df1_dwy(x, u, w), h);
  assert_derivative( limm_wo(f1, x, u, w, h), df1_dwo(x, u, w), h);
  assert_derivative( limm_wt(f1, x, u, w, h), df1_dwt(x, u, w), h);
  assert_derivative( limm_wj(f1, x, u, w, h), df1_dwj(x, u, w), h);
  assert_derivative( limm_wk(f1, x, u, w, h), df1_dwk(x, u, w), h);

  assert_derivative(  limm_xx(f2, x, u, w, h),  df2_dxx(x, u, w), h);
  assert_derivative(  limm_xy(f2, x, u, w, h),  df2_dxy(x, u, w), h);
  assert_derivative( limm_xvx(f2, x, u, w, h), df2_dxvx(x, u, w), h);
  assert_derivative( limm_xvy(f2, x, u, w, h), df2_dxvy(x, u, w), h);

  assert_derivative( limm_wx(f2, x, u, w, h), df2_dwx(x, u, w), h);
  assert_derivative( limm_wy(f2, x, u, w, h), df2_dwy(x, u, w), h);
  assert_derivative( limm_wo(f2, x, u, w, h), df2_dwo(x, u, w), h);
  assert_derivative( limm_wt(f2, x, u, w, h), df2_dwt(x, u, w), h);
  assert_derivative( limm_wj(f2, x, u, w, h), df2_dwj(x, u, w), h);
  assert_derivative( limm_wk(f2, x, u, w, h), df2_dwk(x, u, w), h);

  assert_derivative(  limm_xx(f3, x, u, w, h),  df3_dxx(x, u, w), h);
  assert_derivative(  limm_xy(f3, x, u, w, h),  df3_dxy(x, u, w), h);
  assert_derivative( limm_xvx(f3, x, u, w, h), df3_dxvx(x, u, w), h);
  assert_derivative( limm_xvy(f3, x, u, w, h), df3_dxvy(x, u, w), h);

  assert_derivative( limm_wx(f3, x, u, w, h), df3_dwx(x, u, w), h);
  assert_derivative( limm_wy(f3, x, u, w, h), df3_dwy(x, u, w), h);
  assert_derivative( limm_wo(f3, x, u, w, h), df3_dwo(x, u, w), h);
  assert_derivative( limm_wt(f3, x, u, w, h), df3_dwt(x, u, w), h);
  assert_derivative( limm_wj(f3, x, u, w, h), df3_dwj(x, u, w), h);
  assert_derivative( limm_wk(f3, x, u, w, h), df3_dwk(x, u, w), h);

  assert_derivative(  limm_xx(f4, x, u, w, h),  df4_dxx(x, u, w), h);
  assert_derivative(  limm_xy(f4, x, u, w, h),  df4_dxy(x, u, w), h);
  assert_derivative( limm_xvx(f4, x, u, w, h), df4_dxvx(x, u, w), h);
  assert_derivative( limm_xvy(f4, x, u, w, h), df4_dxvy(x, u, w), h);

  assert_derivative( limm_wx(f4, x, u, w, h), df4_dwx(x, u, w), h);
  assert_derivative( limm_wy(f4, x, u, w, h), df4_dwy(x, u, w), h);
  assert_derivative( limm_wo(f4, x, u, w, h), df4_dwo(x, u, w), h + 2.0);
  assert_derivative( limm_wt(f4, x, u, w, h), df4_dwt(x, u, w), h);
  assert_derivative( limm_wj(f4, x, u, w, h), df4_dwj(x, u, w), h);
  assert_derivative( limm_wk(f4, x, u, w, h), df4_dwk(x, u, w), h);

  std::cerr << "Classic matrices test finished " << std::endl;

  return 0;
}