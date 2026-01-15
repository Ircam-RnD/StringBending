#include "NlBendProcessor.h"
#include "iostream"

template <class T>
NlBendProcessor<T>::NlBendProcessor(float sampleRate, int Nmodes){
  
  this->Nmodes = Nmodes;
  this->Nins = Nmodes;
  this->Nouts = Nmodes;


  ReinitDsp(sampleRate);
};

template <class T>
void NlBendProcessor<T>::ReinitDsp(float sampleRate){
  sr = sampleRate;
  dt = 1 / (float(sr));

  // Reinit state
  qnow = Eigen::Array<T, -1, 1>::Zero(Nmodes);
  qnext = qnow;
  qlast = qnow;
  qspat = qnow;

  r = 0;

  RHS = LHS = qnow;

  // Reinit nonlinear variables
  g = qnow;
  dqV = qnow;
  V = 0;

  // Reinit system matrices
  Eigen::Array<T, -1, 1> Mcopy, Kcopy, Rcopy;
  Mcopy = M;
  Kcopy = K;
  Rcopy = R;

  M = Eigen::Array<T, -1, 1>::Ones(Nmodes);
  K = M;
  R = M;

  Amps = M;
  Omega = 2 * M_PI * 100 * Amps;
  Decays = M;

  if (Mcopy.size() != 0){
    int maxsize = std::max(Mcopy.size(), M.size());
    M.head(maxsize) = Mcopy.head(maxsize);
    K.head(maxsize) = Kcopy.head(maxsize);
    R.head(maxsize) = Rcopy.head(maxsize);
  }
};

template <class T>
void NlBendProcessor<T>::setModalMatrices(){
  // For force-velocity transfer function
  R = Amps.cwiseInverse();
  M = 1/(2*6.9) * Decays.cwiseProduct(R);
  K = M.cwiseProduct(Omega).cwiseProduct(Omega);

  // For force-displacement transfer function
  // R = Omega.cwiseInverse();
  // M = 1/(2*6.9) * Decays.cwiseProduct(R);
  // K = M.cwiseProduct(Omega).cwiseProduct(Omega);
};

template <class T>
void NlBendProcessor<T>::setLinearParameters(Eigen::Array<T, -1, 1> Amps, Eigen::Array<T, -1, 1> Freqs, Eigen::Array<T, -1, 1> Decays){
  SafeSetEigen(this->Amps, Amps);
  // Omega is clamped to ensure stability of the scheme
  SafeSetEigen(this->Omega, ClipEigen(2 * M_PI * Freqs, T(0), T(2.0 * sr)).eval());
  SafeSetEigen(this->Decays, Decays);
  setModalMatrices();
};

template <class T>
void NlBendProcessor<T>::setAmps(Eigen::Array<T, -1, 1> Amps){
  SafeSetEigen(this->Amps, Amps);
  setModalMatrices();
};


template <class T>
void NlBendProcessor<T>::setFreqs(Eigen::Array<T, -1, 1> Freqs){
  SafeSetEigen(this->Omega, ClipEigen(2 * M_PI * Freqs, T(0), T(2.0 * sr)).eval());
  setModalMatrices();
};

template <class T>
void NlBendProcessor<T>::setDecays(Eigen::Array<T, -1, 1> Decays){
  SafeSetEigen(this->Decays, Decays);
  setModalMatrices();
};

template <class T>
void NlBendProcessor<T>::computeVAndVprime(){
  switch (nlMode){
    case LINEAR:
      break;
    case MODEWISE:
      V = (K * qnow.pow(4)).sum() / 4;
      dqV = (K * qnow.pow(3));
      break;
    case SUM:
      V = (qnow*qnow).sum() * (qnow * K * qnow).sum() / 4;
      dqV = 0.5 * (
        qnow * (qnow * K * qnow).sum()
        + K * qnow * (qnow*qnow).sum()
      );
      break;
  };
};

template <class T>
void NlBendProcessor<T>::computeV(){
  switch (nlMode){
    case LINEAR:
      break;
    case MODEWISE:
      V = (K * ((qnow+qlast)/2).pow(4)).sum() / 4;
      break;
    case SUM:
      V = (((qnow+qlast)/2)*((qnow+qlast)/2)).sum() * (((qnow+qlast)/2) * K * ((qnow+qlast)/2)).sum() / 4;
      break;
  };
};

template <class T>
void NlBendProcessor<T>::process(Eigen::Ref<const Eigen::Array<T, -1, 1>> input, Eigen::Ref<Eigen::Array<T, -1, 1>> out, T &epsilonOut){
  // Nonlinear part
  computeVAndVprime();
  g = dqV / (sqrt(2*V) + NUM_EPS);
  if (V>maxV){
    maxV = V;
  }

  if (controlTerm){
    computeV();
    epsilon = r - sqrt(2*V);
    g -= lambda0 * epsilon * dt * ((qnow-qlast)>0).select(Eigen::Array<T, -1, 1>::Ones(Nmodes), -Eigen::Array<T, -1, 1>::Ones(Nmodes)) / ((qnow-qlast).matrix().template lpNorm<1>() + NUM_EPS);
  }

  // Linear part
  RHS = (- K * dt * dt + 2 * M) * qnow
    - (M - R * dt / 2) * qlast
    + input * dt * dt;

  // Nonlinear part
  RHS += pow(dt, 2) * (0.25 * g * (g * qlast).sum() - g * r);

  // Solve using Shermann-Morrisson
  auto inter = (M + dt * R/2).cwiseInverse(); // TODO: Precompute (and allocate) vector
  qnext = inter * RHS - 0.25 * pow(dt, 2) *(inter * g * (inter*g*RHS).sum()) 
    / (1 + 0.25 * pow(dt, 2) * (g * inter * g).sum());
  
  r = r + 0.5 * (g*(qnext - qlast)).sum();

  qlast = qnow;
  qnow = qnext;

  out = (qnow-qlast) / dt;
  // out = qnow;
  epsilonOut = epsilon / (sqrt(2*maxV)+NUM_EPS);
};

template class NlBendProcessor<double>;
template class NlBendProcessor<float>;
