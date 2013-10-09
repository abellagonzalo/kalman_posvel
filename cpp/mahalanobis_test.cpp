#include <vector>
#include "kalman.h"

namespace kalman {

    namespace selector {

        namespace {

            template <typename T>
            MatrixCM diff(const T& t1, const T& t2) {
                return t1 - t2;
            }

        } // anonymous

        double mahalanobis(const MatrixCM& x1, const MatrixCM& x2, MatrixCM cov) {
            return sqrt( (diff(x1, x2).transpose() * cov.inverse() * diff(x1, x2)).determinant() );
        }

    } // selector

    class MultiKalman {
    public:
        MultiKalman(int n)
        : n_(n) 
        {

        }

    private:
        int n_;
        std::vector<algorithm::MatrixCMPair> pairs;
    };

} // kalman


int main(int argc, char** argv) {    



    return 0;
}