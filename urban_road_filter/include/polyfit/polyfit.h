/*
 * Copyright 2012 Roland Kwitt
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License.  You may obtain a copy
 * of the License at
 *
 *    http:www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  See the
 * License for the specific language governing permissions and limitations
 * under the License.
 */

#ifndef POLYFIT_H
#define POLYFIT_H

// #include "polyfit/logging.h"

#include <Eigen/Dense>

#include <set>
#include <vector>

// Random number generation
#include <boost/random/mersenne_twister.hpp>

#if BOOST_VERSION >= 104700
	#include <boost/random/uniform_int_distribution.hpp>
#else
	#include <boost/random/uniform_int.hpp>
#endif


// extern logging logger;

template <typename T>
class PolyFit
{
  public:

    enum SolverType 
    {
      EIGEN_COL_PIV_HOUSEHOLDER,
      EIGEN_FULL_PIV_HOUSEHOLDER,   
      EIGEN_JACOBI_SVD
    };

    struct Options
    {
      Options() : 
        solver(EIGEN_JACOBI_SVD), 
        polyDeg(2),
        maxTrial(100),
        maxMPts(5),
        tolerance(0.1) {}

      SolverType solver; // Default: SVD solver
      int polyDeg; // Default: 2-nd deg. poly.
      int maxTrial; // Max. RANAC iterations
      int maxMPts; // Max. #pts to fit polynomial
      double tolerance; // Pt. to model thr. for inliers
    };

  private:

#if BOOST_VERSION >= 104700
    // RNG
    boost::random::mt19937 gen;
    typedef boost::random::uniform_int_distribution<> DistributionType;
#else
    // RNG
    boost::mt19937 gen;
    typedef boost::uniform_int<int> DistributionType;
#endif

    // Define own matrix type, based on template parameter T
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatType;

		// Define own vector type, based on template parameter T
		typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VecType;
    
    // Polynomial coefficients
    VecType pCoef;

    std::vector<T> &x; // Sample points
    std::vector<T> &y; // Signal values

    // Sample + Signal in Eigen format
    VecType xM;
    VecType yM;

    // Vandermonde matrix
    MatType V;

    // Fitting options
    Options &options;

  public:

    // Public accessor methods
    int getDegree() const { return options.polyDeg; }

    SolverType getSolverType() const { return options.solver; }

    const VecType getCoefficients() const { return pCoef; }

    // Constructor
    PolyFit(std::vector<T> &x, std::vector<T> &y, Options &options): 
      x(x), y(y), options(options) 
    {
      if (!x.size() || !y.size() || (x.size() != y.size()))
      {
        // logger.log(ERROR, "[%s]: Size mismatch (%d,%d)\n",
        //     __func__, x.size(), y.size());
        throw std::exception();
      }

      if (options.polyDeg < 1)
      {
        // logger.log(ERROR, "[%s]: Polynomial degree < 1\n",
        //     __func__);
        throw std::exception();
      }
      dataToEigen();
    }

    // Evaluate polynomial at positions in 'at'
    std::vector<T> evalPoly(const std::vector<T> &at)
    { 
      VecType x;
      VecType y;

      x.resize(at.size());
      y.resize(at.size());

      for (int i=0; i<at.size(); ++i)
        x(i,0) = at[i];

      // Horner method
      if (pCoef.rows() > 0)
        y.setConstant(pCoef(0,0));

      int nc = pCoef.rows();
      for (int i=1; i<nc; ++i)
        y = (x.cwiseProduct(y)).array() + pCoef(i,0);

      // Copy data over to output vector
      std::vector<T> out(y.rows());
      for (int i=0;i<y.rows();++i)
        out[i] = y(i,0);

      return out;
    }


    // Points where signal is sampled
    void setSamplePoints(std::vector<T> &v)
    {
      if (!v.size())
      {
        // logger.log(ERROR, "[%s]: No new x-values!\n",
        //     __func__);
        throw std::exception();
      }

      if (x.size())
      {
        // logger.log(WARNING, "[%s]: Resetting x-val.\n",
        //     __func__);
      }
      x = v;
      dataToEigen();
    }


    // Signal value at sample points
    void setSignalPoints(std::vector<T> &v)
    {
      if (!v.size())
      {
        // logger.log(ERROR, "[%s]: No new y-values!\n",
        //     __func__);
        throw std::exception();
      }

      if (y.size())
      {
        // logger.log(WARNING, "[%s]: Resetting y-val.\n",
        //     __func__);
      }
      y = v;
      dataToEigen();
    }


    // Least-Squares (LS) 
    T solveLS(void)
    {
      // logger.log(DEBUG, "[%s]: LS-solving a (x,y)=(%d,%d) problem!\n",
      //     __func__, xM.rows(), yM.rows());

      // Size check    
      if (xM.rows() != yM.rows())
      {
        // logger.log(ERROR, "[%s]: Size mismatch!\n",
        //     __func__);
        throw std::exception();
      }

      buildV();

      switch (options.solver)
      {
        case EIGEN_COL_PIV_HOUSEHOLDER:
          pCoef = V.colPivHouseholderQr().solve(yM.col(0));  
          break;
        case EIGEN_FULL_PIV_HOUSEHOLDER:
          pCoef = V.fullPivHouseholderQr().solve(yM.col(0));
          break;
        case EIGEN_JACOBI_SVD:
          pCoef = V.jacobiSvd(Eigen::ComputeThinU | 
              Eigen::ComputeThinV).solve(yM.col(0));
          break;
      }
      return getMSE();
    }


    // Mean-Square Error (MSE)
    T getMSE(void)
    {
      assert(V.cols() == pCoef.rows() && 
          yM.rows() == V.rows());

      T mse = 0.0;
      Eigen::VectorXd diff = ((V*pCoef).array() - yM.col(0).array());
      for (int i=0; i<diff.rows(); ++i)
        mse += pow(diff(i),2);
      return mse;
    }


    // Sample s points from [1...N]
    void randsample(int N, int s, std::set<int> &ind)
    {
      assert(N > 0 && s > 0 && (N>s));

      DistributionType dist(1, N);
      int count = 0;
      while(static_cast<int>(ind.size()) < s)
      { 
        ind.insert(dist(gen));      
        count++;
        if(count > s){
          break;
        }
      }
    }


    // RANSAC-based solution
    void solveRLS(void)
    {
      std::vector<int> bestSubset;
      Eigen::Matrix<T, Eigen::Dynamic, 1> bestCoef;
      
      // Probability that we can at least choose
      // one sample with options.maxMPts free of
      // outliers (sugg. by Peter Kovesi's code)
      double p = 0.99;

      int trialCount = 0; int N = 1;
      while (N > trialCount)
      {
        // logger.log(DEBUG, "[%s]: Trial %d\n",
        //     __func__, trialCount);         
        std::set<int> indexSet;
        randsample(static_cast<int>(y.size()-1),
            options.maxMPts, indexSet);
         
        assert(static_cast<int>(indexSet.size()) == options.maxMPts);
        
        // Build (sample,signal) for current subset
        dataToEigenSubset(indexSet);
        
        buildV();
        
        
        // Solve the LS-fitting problem
        double mse = solveLS();
      
        
        //  logger.log(DEBUG, "[%s]: MSE=%.8f in %d-th RLS trial\n",
        //     __func__, mse, trialCount);

        // Compute inlier indices
        std::vector<int> inliers;
        int nInliers = getInliers(inliers);  
        // logger.log(DEBUG, "[%s]: #Inliers = %d\n",
        //     __func__, nInliers);

        if (inliers.size() > bestSubset.size())
        {
          // Store settings
          bestSubset = inliers;
          bestCoef = pCoef;

          // Refine estimate of N (based on Peter Kovesi's MATLAB impl. of RANSAC)
          double fracInliers = static_cast<double>(inliers.size())/double(yM.rows());
          double pNoOutliers = 1 - pow(fracInliers,options.maxMPts);
          pNoOutliers = std::max(1e-9, pNoOutliers);  
          pNoOutliers = std::min(1-1e-9, pNoOutliers);
          N = log(1-p)/log(pNoOutliers);

          // logger.log(DEBUG, "[%s]: Inlier fraction=%.2f, P(Choose no outlier)=%.2f --> N=%d\n",
          //     __func__, fracInliers, pNoOutliers, N);
        }

        ++trialCount;
   
        if (trialCount > options.maxTrial)
        {
          // logger.log(WARNING, "[%s]: Maximum nr. of trials (%d) reached!\n",
          //     __func__, options.maxTrial);
          break;
        }
      }
      pCoef = bestCoef;
    }

  private:

    void dataToEigenSubset(std::set<int> &indexSet)
    {
      int nVal = static_cast<int>(indexSet.size());
      xM.resize(nVal,1);
      yM.resize(nVal,1);
      
      int cnt = 0;  
      std::set<int>::iterator it = indexSet.begin();
      while (it != indexSet.end())
      {
        xM(cnt,0) = x[*it];
        yM(cnt,0) = y[*it];
        ++cnt; ++it;
      }
    }

    // Determine model inliers, based on max. allowed
    // distance to fitted polynomial
    int getInliers(std::vector<int> &inliers)
    {
      // Full data to Eigen
      dataToEigen();
      buildV();

      Eigen::VectorXd diff = ((V*pCoef).array() - yM.col(0).array());
      for (int i=0; i<diff.rows(); ++i)
      {
        if (pow(diff(i),2) < options.tolerance)
          inliers.push_back(i);
      }
      return static_cast<int>(inliers.size());
    }
    

    // All (sample,signal) data to Eigen format
    void dataToEigen(void)
    {
      xM.resize(static_cast<int>(x.size()),1);
      yM.resize(static_cast<int>(y.size()),1);

      for (int i=0;i<static_cast<int>(x.size()); ++i)
        xM(i,0) = x[i];

      for (int i=0;i<static_cast<int>(y.size()); ++i)
        yM(i,0) = y[i];
    }


    // Vandermonde matrix
    void buildV(void)
    {
      V.resize(xM.rows(), options.polyDeg + 1);

      for (int i=0; i<xM.rows(); ++i)
        V(i,options.polyDeg) = 1;

      for (int d=0; d<options.polyDeg; ++d) 
        for (int i=0; i<xM.rows(); ++i) 
          V(i,d) = pow(xM(i,0), options.polyDeg - d);
    }
};


// For convenience
template <typename T>
std::ostream& operator<< (std::ostream &out, PolyFit<T> &f)
{
   out << "Polynomial coefficients for " << 
   f.getDegree() << "-nd order:"  << std::endl;     
   out << f.getCoefficients() << std::endl;
   out << "Solver type: " << f.getSolverType();
   return out;
}

#endif
