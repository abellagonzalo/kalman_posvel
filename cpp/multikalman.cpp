#include <iostream>
#include <assert.h>
#include <algorithm>
#include <vector>
#include <boost/bind.hpp>
#include "Matrix.h"
#include "classic_matrices.h"


namespace kalman {

    struct Kalman {
        bool operator==(const Kalman& rhs) const {
            return x == rhs.x && P == rhs.P;
        }
        bool operator!=(const Kalman& rhs) const {
            return !(*this == rhs);
        }
        MatrixCM x;
        MatrixCM P;
    };

    Kalman kalman_generator() {
        double x [4] = { 0, 0, 0, 0 };
        double P [16] = { 1000, 0, 0, 0,
                         0, 1000, 0, 0,
                         0, 0, 1000, 0,
                         0, 0, 0, 1000 };
        Kalman kalman = { MatrixCM(4, 1, x), MatrixCM(4, 4, P)};
        return kalman;
    }

    void print_kalman(const Kalman& k) {
        std::cerr << "*** x ***";
        k.x.print();
        std::cerr << "*** P ***";
        k.P.print();
    }

    double zero4_init [4] = { 0., 0., 0., 0.};
    double zero6_init [6] = { 0., 0., 0., 0., 0., 0.};
    MatrixCM zero4(4, 1, zero4_init);
    MatrixCM zero6(6, 1, zero6_init);
    MatrixCM iden4(4);

    Kalman predict_f(const Kalman& k, double* u_init, double* w_init) {
        using namespace matrices;
        MatrixCM u(6, 1, u_init);
        MatrixCM w(6, 1, w_init);
        Kalman out = { f(k.x, u, zero6),
                       A(k.x, u, w) * k.P * A(k.x, u, w).transpose() + W(k.x, u, w) * Q(w) * W(k.x, u, w).transpose() };
        return out;
    }

    template <typename T>
    MatrixCM diff(const T& t1, const T& t2) {
        return t1 - t2;
    }

    double mahalanobis_mat(const MatrixCM& x1, const MatrixCM& x2, MatrixCM cov) {
        return sqrt( (diff(x1, x2).transpose() * cov.inverse() * diff(x1, x2)).determinant() );
    }

    double mahalanobis(const Kalman& k, double* z) {
        double x1  [2] = { k.x.e(0,0), k.x.e(1,0) };
        double cov [4] = { k.P.e(0,0), k.P.e(0,1),
                           k.P.e(1,0), k.P.e(1,1)};
        return mahalanobis_mat(MatrixCM(2,1,x1), MatrixCM(2,1,z), MatrixCM(2,2,cov));
    }

    bool min_mahalanobis(const Kalman& k1, const Kalman& k2, double* z) {
        std::cerr << mahalanobis(k1, z) << " - " << mahalanobis(k2, z) << std::endl;
        return mahalanobis(k1, z) < mahalanobis(k2, z);
    }

    Kalman update_f(const Kalman& k, double* z_init, double* v_init) {
        using namespace matrices;
        MatrixCM z(2, 1, z_init);
        MatrixCM v(2, 1, v_init);
        MatrixCM aux = H(k.x, v) * k.P * H(k.x, v).transpose() + V(k.x, v) * R(v) * V(k.x, v).transpose();
        MatrixCM gain = k.P * H(k.x, v).transpose() * aux.inverse();
        Kalman out = { k.x + gain * ( z - h(k.x, zero4) ),
                     ( iden4 - gain * H(k.x,v) ) * k.P };
        return out;
    }

    class MultiKalman {
    public:
        typedef std::vector<Kalman>::iterator iterator;
        typedef std::vector<Kalman>::const_iterator const_iterator;

        MultiKalman(int min, int max)
        : min_(min), max_(max), kalmans_(min, kalman_generator())
        {
            kalmans_.reserve(max);
        }

        int min() const
        { 
            return min_; 
        }

        int max() const
        {
            return max_; 
        }

        int size() const
        {
            return kalmans_.size();
        }

        iterator begin() {
            return kalmans_.begin();
        }

        const_iterator begin() const {
            return kalmans_.begin();
        }

        iterator end() {
            return kalmans_.end();
        }

        const_iterator end() const {
            return kalmans_.end();
        }

        void predict(double* u, double* w)
        {
            std::transform(kalmans_.begin(), kalmans_.end(), kalmans_.begin(), boost::bind(predict_f, _1, u, w));
        }

        void update(double* z, double* v)
        {
            if (kalmans_.empty()) return;
            std::vector<Kalman>::iterator it = std::min_element(kalmans_.begin(), kalmans_.end(),
                                                                boost::bind(min_mahalanobis, _1, _2, z));
            *it = update_f(*it, z, v);
        }

    private:
        int min_;
        int max_;
        std::vector<Kalman> kalmans_;
    };

}

bool is_equal_to_kalman_generator(const kalman::Kalman& k) {
    return k == kalman::kalman_generator();
}

int main(int argc, char** argv) {

    // Create kalman filters
    kalman::MultiKalman multi(1,3);

    // Check min && max
    assert( multi.min() == 1 );
    assert( multi.max() == 3 );

    // Check is initialized correctly
    assert( std::count_if(multi.begin(), multi.end(), is_equal_to_kalman_generator) == 1 );

    int imax = 100;
    for (int i=0; i<imax; ++i) {
        double u [6] = { 0, 0, 0, 1./30., 0, 0.9 };
        double w [6] = { 0.01, 0.01, 0.01, 0, 10, 0 };
        multi.predict(u, w);

        double z [2] = { 100 + 100*i, 0 };
        double v [2] = { .1, .1 };
        multi.update(z, v);

        double z2 [2] = { 0, 100 + 100*i };
        double v2 [2] = { .1, .1 };
        // multi.update(z2, v2);
    }
    std::for_each(multi.begin(), multi.end(), kalman::print_kalman);

    std::cout << "Finish!!" << std::endl;

    return 0;
}