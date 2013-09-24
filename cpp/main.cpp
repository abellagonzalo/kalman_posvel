#include "Matrix.h"
#include <cmath>

MatrixCM kalman_predict_state(const MatrixCM& mov_model, const MatrixCM& state,
                              const MatrixCM& rot_model, const MatrixCM& tra_model) {
  return rot_model*(mov_model*state) + tra_model;
}

MatrixCM kalman_predict_estimation(const MatrixCM& estimation,
                                   const MatrixCM& mov_model,
                                   const MatrixCM& pred_noise) {
  MatrixCM transpose(mov_model);
  transpose.transpose();
  return mov_model * estimation * transpose + pred_noise;
}

MatrixCM kalman_gain(const MatrixCM& estimation,
                     const MatrixCM& transition,
                     const MatrixCM& obs_noise) {
  MatrixCM transpose(transition);
  transpose.transpose();
  MatrixCM inverse(transition * estimation * transpose + obs_noise);
  inverse.inverse();
  return transition * transpose * inverse;
}

MatrixCM kalman_update_state(const MatrixCM& state,
                             const MatrixCM& gain,
                             const MatrixCM& observation,
                             const MatrixCM& transition) {
  return state + gain * (observation - transition * state);
}

MatrixCM kalman_update_estimation(const MatrixCM& estimation,
                                  const MatrixCM& gain,
                                  const MatrixCM& transition) {
  return (MatrixCM(4) - gain * transition) * estimation;
}

MatrixCM create_movement_model(int last_ts, int ts) {
  float t = ts-last_ts;
  float data[4*4] = {1., 0., t,  0.,
                     0., 1., 0., t,
                     0., 0., 1., 0.,
                     0., 0., 0., 1.};
  return MatrixCM(4,4,data);
}

MatrixCM create_rotation_model(float alpha) {
  float c = cos(alpha);
  float s = sin(alpha);
  float data[4*4] = {c,  s,  0., 0.,
                     -s, c,  0., 0.,
                     0., 0., c,  s,
                     0., 0., -s, c};
  return MatrixCM(4,4,data);
}

MatrixCM create_translation_model(double x, double y) {
  float data[4] = {-x, -y, 0., 0.};
  return MatrixCM(4,1,data);
}

class BicaKalman {

public:
  void timestamp(long timestamp) {

  }

  void odometry(int x, int y, double alpha) {

  }

protected:
  MatrixCM movementModel() {
    return MatrixCM();
  }

  MatrixCM rotationModel() {
    return MatrixCM();
  }

  MatrixCM translationModel() {
    return MatrixCM();
  }

  MatrixCM predictionNoise() {
    return MatrixCM();
  }

  MatrixCM transitionObs2State() {
    return MatrixCM();
  }

  MatrixCM ObservationsNoise() {
    return MatrixCM();
  }

  MatrixCM observations() {
    return MatrixCM();
  }
};

template <class T>
class GenericKalman : public T 
{
public:
  GenericKalman() { }
  ~GenericKalman() { }

  void predict() {
    state = kalman_predict_state(T::movementModel(), state, T::rotationModel(), T::translationModel());
    estimation = kalman_predict_estimation(estimation, T::movementModel(), T::predictionNoise());
  }

  void update() {
    MatrixCM gain = kalman_gain(estimation, T::transitionObs2State(), T::ObservationsNoise());
    state = kalman_update_state(state, gain, T::observations(), T::transitionObs2State());
    estimation = kalman_update_estimation(estimation, gain, T::transitionObs2State());
  }

  void predict(long time, int movx, int movy, int alpha) {

  }

  void update(int x, int y, int vx, int vy) {
    
  }

  double x() const {
    return state.e(0,0);
  }

  double y() const {
    return state.e(1,0);
  }

  double vx() const {
    return state.e(2,0);
  }

  double vy() const {
    return state.e(3,0);
  }

  double estimationX() {
    return estimation.e(0,0);
  }
  double estimationY() {
    return estimation.e(1,1);
  }

  double estimationVX() {
    return estimation.e(2,2);
  }

  double estimationVY() {
    return estimation.e(3,3);
  }

private:
  MatrixCM state;
  MatrixCM estimation;

  MatrixCM mov_model_;
  MatrixCM rot_model_;
  MatrixCM trans_model;

  MatrixCM transition_;

  MatrixCM pred_noise;
  MatrixCM obs_noise;  
};

class Tracker {
public:
  std::pair<double, double> position() const {
    return std::make_pair(0.,0.);
  }

  std::pair<double, double> velocity() const {
    return std::make_pair(0.,0.);
  }

  void update(std::pair<double, double> newPosition, long timestamp) {
    velocity_ = make_pair(calc_vel(newPosition.first, position_.first, timestamp, timestamp_),
                          calc_vel(newPosition.second, position_.second, timestamp, timestamp_));

  };

private:
  double calc_vel(double newPos, double lastPos, long newTime, long lasttime) {
    return (newPos - lastPos) / double(newTime - lasttime);
  }

  long timestamp_;
  std::pair<double, double> position_;
  std::pair<double, double> velocity_;
};

int main(int argc, char **argv) {

  GenericKalman<BicaKalman> kalman;

  float data[4] = {100., 100., 100., 100.};
  MatrixCM state(4,1,data);

  MatrixCM estimation(4);

  float prediction_noise_data[4*4] = {.2, .2, .2, .2,
                                      .2, .2, .2, .2,
                                      .2, .2, .2, .2,
                                      .2, .2, .2, .2};
  MatrixCM predictionNoise(4,4,prediction_noise_data);

  MatrixCM newState = kalman_predict_state(create_movement_model(0,1), 
                                           state,
                                           create_rotation_model(M_PI/2), 
                                           create_translation_model(0, 0));
  newState.print();

  return 0;
}