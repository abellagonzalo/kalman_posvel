namespace kalman {

    double mahalanobis_distance() {
        return 0.0;
    }

    class MultiKalman {
    public:
        MultiKalman(int min, int max)
        : min_(min)
        , max_(max)
        { }

    private:
        int min_;
        int max_;
    };

}

int main(int argc, char** argv) {

    kalman::MultiKalman multiKalman(1,3);

    return 0;
}