
#ifndef NL_BEND_PROCESSOR_H
#define NL_BEND_PROCESSOR_H

#include <cmath>
#include <vector>
#include <tuple>

#include <Eigen/Dense>
#include <string_view>

#include "EigenUtility.h"

enum NLMODE {
    LINEAR,
    MODEWISE,
    SUM
    // NEAREST
};

template <class T>
class NlBendProcessor {
    private:
        using Vec = Eigen::Array<T, -1, 1>;
        // Numerical epsilon value
        constexpr static T NUM_EPS{1e-12};
        // Number of modes
        int Nmodes{1};

        // Linear part: system "matrices" (diagonal)
        Vec M, K, R;
        // Higer level modal parameters
        Vec Amps, Omega, Decays;

        // Nonlinear part: function parametrization: 
        // How can we do that?
        NLMODE nlMode = LINEAR;


        // Time-scheme parameters
        float sr;
        T dt;
        bool controlTerm{true};
        T lambda0{0};

        // System state (modal coordinates)
        Vec qlast, qnow, qnext;
        T r;

        // Intermediate vectors
        Vec RHS, LHS;

        // Displacement (spatial coordinates)
        Vec qspat;

        // Nonlinear function evaluation
        Vec g, dqV;
        T V;

        // Drift variable
        T epsilon{0}, maxV{0};

        // Input and output dimensions
        int Nins{1}, Nouts{1};
        
        void setModalMatrices();

    public:
        NlBendProcessor(float sampleRate, int Nmodes = 1);

        void ReinitDsp(float sampleRate);

        void computeVAndVprime();

        void computeV();

        void process(Eigen::Ref<const Vec> input, Eigen::Ref<Vec> out, T &epsilonOut);

        // Higher level modal parameters
        void setLinearParameters (Vec Amps, Vec Omega, Vec Decays);
        void setAmps (Vec Amps);
        void setFreqs (Vec Freqs);
        void setDecays (Vec Decays);
    
        void setNlMode (NLMODE mode) {this->nlMode = mode;};

        T getOmega(int i){
            return this->K(i);
        }
        
        void setLambda0(T lambda0){this->lambda0 = std::clamp(lambda0, T(0), T(10000));};
        // Discretization parameters
        int getNmodes() {return Nmodes;};

        static constexpr int getNins(int Nmodes){return Nmodes;}
        static constexpr int getNouts(int Nmodes) {return Nmodes;}
};

#endif
