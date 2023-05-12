/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Kinetx/RandomNetworkFactory.h"
#include "Kinetx/Network.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <random>

namespace Scine {
namespace Kinetx {

Network RandomNetworkFactory::random() {
  // All values in SI units:
  unsigned int nSamples = 26;           // B+1 in the original paper
  unsigned int maxEdges = 125;          // L_max in the original paper
  double avgCovarianceDeviation = 2.09; // sigma_A in the original paper
  double muBi = 2;                      // mu_exp(N_bi) in the original paper
  double muUni = 5;                     // mu_exp(N_uni) in the original paper
  double pNew = 0.5;                    // p(V_new) in the original paper
  double vertexMean = -25.0;            // mu(A^--A^+)  in the original paper
  double vertexStdDev = 50.0;           // sigma(A^--A^+) in the original paper
  double temperature = 298.15;          // T in the original paper
  double minEdgeVal = 0.0;              // min(A^--A^+) in the original paper
  double maxEdgeVal = 100.0;            // max(A^--A^+) in the original paper
  double nInitialVertices = 2;
  unsigned int seed = 123456;

  unsigned int nVertices = 0; // N in the original paper
  unsigned int nEdges = 0;    // L in the original paper

  // Initialize vertices
  std::vector<double> vertexMasses;
  vertexMasses.reserve(maxEdges);
  for (unsigned int i = 1; i < nInitialVertices + 1; i++) {
    // Add initial vertices
    vertexMasses.emplace_back(double(i));
    nVertices++;
  }

  // Initialize Edge data
  using ESMi = Eigen::SparseMatrix<int>;
  ESMi sf(maxEdges, maxEdges);
  ESMi sb(maxEdges, maxEdges);
  sf.reserve(2 * maxEdges);
  sb.reserve(2 * maxEdges);

  /*
   * Generate new Edges until the requested amount is reached
   */
  std::default_random_engine generator{seed};
  std::exponential_distribution<double> exprndUni(1.0 / muUni);
  std::exponential_distribution<double> exprndBi(1.0 / muBi);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  unsigned int lastNewVertex = 0;
  while (nEdges < maxEdges) {
    // Keep track of new and old vertices
    const unsigned int nOldVertices = lastNewVertex;
    lastNewVertex = nVertices;
    // Uni molecular reactions
    for (unsigned int v = nOldVertices; v < lastNewVertex && nEdges < maxEdges; v++) {
      // `1 +` taken from source code the paper does not read that way
      const unsigned int nUni = floor(exprndUni(generator));
      // Find all vertces with the same mass
      std::vector<int> sameMass;
      if (nUni > 0) {
        for (unsigned int j = 0; j < v; j++) {
          if (vertexMasses[j] == vertexMasses[v]) {
            sameMass.push_back(j);
          }
        }
      }
      for (unsigned int i = 0; i < nUni && nEdges < maxEdges; i++) {
        bool generateNew{uniform(generator) > pNew};
        if (generateNew) {
          // Generate new vertex with same mass
          vertexMasses.emplace_back(vertexMasses[v]);
          nVertices++;
          // Add edge
          sf.insert(nEdges, v) = 1;
          sb.insert(nEdges, nVertices - 1) = 1;
          nEdges++;
        }
        else {
          // Exclusively back connect to old vertices with the same mass
          if (sameMass.size() == 0) {
            continue;
          }
          // Pick old element
          std::uniform_int_distribution<int> dist(0, sameMass.size() - 1);
          const unsigned int old = dist(generator);
          // Add edge
          sf.insert(nEdges, v) = 1;
          sb.insert(nEdges, sameMass[old]) = 1;
          sameMass.erase(sameMass.begin() + old);
          nEdges++;
        }
      }
    }

    // Bi molecular reactions
    for (unsigned int v1 = nOldVertices; v1 < lastNewVertex && nEdges < maxEdges; v1++) {
      for (unsigned int v2 = 0; v2 < lastNewVertex && nEdges < maxEdges; v2++) {
        // `1 +` taken from source code the paper does not read that way
        const unsigned int nBi = floor(exprndUni(generator));
        // Find all vertices with the same mass as v1+v2
        std::vector<int> sameMass;
        const double mass = vertexMasses[v1] + vertexMasses[v2];
        if (nBi > 0) {
          for (unsigned int j = 0; j < vertexMasses.size(); j++) {
            if (vertexMasses[j] == mass) {
              sameMass.push_back(j);
            }
          }
        }
        // Find all old vertice pairs with the same mass as v1+v2
        std::vector<std::pair<int, int>> sameMassPairs;
        if (nBi > 0) {
          for (unsigned int j = 0; j < nOldVertices; j++) {
            for (unsigned int k = 0; k <= j; k++) {
              if ((vertexMasses[j] + vertexMasses[k]) == mass) {
                sameMassPairs.push_back({int(j), int(k)});
              }
            }
          }
        }
        for (unsigned int i = 0; i < nBi && nEdges < maxEdges; i++) {
          bool generateNew{uniform(generator) > pNew};
          if (generateNew) {
            // Generate new vertex with same mass as v1+v2
            vertexMasses.emplace_back(vertexMasses[v1] + vertexMasses[v2]);
            nVertices++;
            // Add edge
            if (v1 != v2) {
              sf.insert(nEdges, v1) = 1;
              sf.insert(nEdges, v2) = 1;
            }
            else {
              sf.insert(nEdges, v2) = 2;
            }
            sb.insert(nEdges, nVertices - 1) = 1;
            nEdges++;
          }
          else {
            // Exclusively back connect to old vertices with the same mass
            if (sameMass.size() == 0 && sameMassPairs.size() == 0) {
              continue;
            }
            // Pick old element
            std::uniform_int_distribution<int> dist(0, sameMass.size() + sameMassPairs.size() - 1);
            const unsigned int old = dist(generator);
            if (old < sameMass.size()) {
              // Add edge
              if (v1 != v2) {
                sf.insert(nEdges, v1) = 1;
                sf.insert(nEdges, v2) = 1;
              }
              else {
                sf.insert(nEdges, v2) = 2;
              }
              sb.insert(nEdges, sameMass[old]) = 1;
              sameMass.erase(sameMass.begin() + old);
              nEdges++;
            }
            else {
              const auto pair = sameMassPairs[old - sameMass.size()];
              // Add edge
              if (v1 != v2) {
                sf.insert(nEdges, v1) = 1;
                sf.insert(nEdges, v2) = 1;
              }
              else {
                sf.insert(nEdges, v2) = 2;
              }
              if (pair.first != pair.second) {
                sb.insert(nEdges, pair.first) = 1;
                sb.insert(nEdges, pair.second) = 1;
              }
              else {
                sb.insert(nEdges, pair.first) = 2;
              }
              sameMassPairs.erase(sameMassPairs.begin() + old - sameMass.size());
              nEdges++;
            }
          }
        } /* nBi tries */
      }
    } /* outer bi-mol */
  }
  sf.conservativeResize(nEdges, nVertices);
  sf.data().squeeze();
  sb.conservativeResize(nEdges, nVertices);
  sb.data().squeeze();

  /*
   * Generate vertex free energies (Compounds)
   */
  // Generate random covariance matrix (sigma)
  srand(seed);
  Eigen::MatrixXd p = 0.5 * Eigen::MatrixXd::Random(nEdges, nEdges);
  Eigen::MatrixXd sigma = p * p.transpose();
  Eigen::VectorXd eval = sigma.selfadjointView<Eigen::Lower>().eigenvalues();
  sigma.array() *= (avgCovarianceDeviation * avgCovarianceDeviation / eval.maxCoeff());
  Eigen::MatrixXd normalNoise(nEdges, nSamples);
  std::normal_distribution<double> stdGauss(0.0, 1.0);
  for (unsigned int i = 0; i < nEdges; i++) {
    for (unsigned int j = 0; j < nSamples; j++) {
      normalNoise(i, j) = stdGauss(generator);
    }
  }
  Eigen::MatrixXd vertexNoise = sigma.llt().matrixL() * normalNoise;

  // Generate normal distributed vertex energies
  std::normal_distribution<double> vertexDist(vertexMean, vertexStdDev);
  Eigen::VectorXd vertexGs = Eigen::VectorXd::Zero(nVertices);
  for (unsigned int i = 0; i < nVertices; i++) {
    vertexGs[i] = vertexDist(generator);
  }

  /*
   * Generate edge free energies (Reactions)
   */
  // Calculate total free energies for both sides of the reactions
  Eigen::VectorXd lhs = Eigen::VectorXd::Zero(nEdges);
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(nEdges);
  Eigen::MatrixXd lhsNoise = Eigen::MatrixXd::Zero(nEdges, nSamples);
  Eigen::MatrixXd rhsNoise = Eigen::MatrixXd::Zero(nEdges, nSamples);
  for (int col = 0; col < sf.outerSize(); ++col) {
    for (Eigen::SparseMatrix<int>::InnerIterator it(sf, col); it; ++it) {
      lhs[col] += it.value() * vertexGs[it.row()];
      for (unsigned int j = 0; j < nSamples; j++)
        lhsNoise(col, j) += it.value() * vertexNoise(it.row(), j);
    }
  }
  for (int col = 0; col < sb.outerSize(); ++col) {
    for (Eigen::SparseMatrix<int>::InnerIterator it(sb, col); it; ++it) {
      rhs[col] += it.value() * vertexGs[it.row()];
      for (unsigned int j = 0; j < nSamples; j++)
        rhsNoise(col, j) += it.value() * vertexNoise(it.row(), j);
    }
  }

  /*
   * Generate rate constants
   */
  // Calculate Delta G (dg) for R->TS (forward, f) and TS->P (backward, b)
  //  (activation free energies)
  using ESMd = Eigen::SparseMatrix<double>;
  ESMd dgf(nEdges, nSamples);
  ESMd dgb(nEdges, nSamples);
  std::uniform_real_distribution<double> edgeDist(minEdgeVal, maxEdgeVal);
  std::normal_distribution<double> edgeNoise(0.0, 2.0 * avgCovarianceDeviation);
  for (unsigned int i = 0; i < nEdges; i++) {
    const double gEdge = std::max(lhs[i], rhs[i]);
    const double baseHeight = edgeDist(generator);
    for (unsigned int j = 0; j < nSamples; j++) {
      const double noise = edgeNoise(generator) + (lhsNoise(i, j) + rhsNoise(i, j)) / 2.0;
      double height = baseHeight + noise;
      if (height < 0.0) {
        height *= -1.0;
      }
      dgf.insert(i, j) = gEdge - lhs[i] + height;
      dgb.insert(i, j) = gEdge - rhs[i] + height;
    }
  }

  // Generate rate constants from activation free energies (Eyring equation)
  ESMd rf(nEdges, nSamples);
  ESMd rb(nEdges, nSamples);
  constexpr double gas = 8.314e-03;      // gas constant [kJ / (mol * K)]
  constexpr double plank = 6.626e-34;    // Planck constant [J * s]
  constexpr double boltzman = 1.381e-23; // Boltzmann constant [J / K]
  for (int k = 0; k < dgf.outerSize(); ++k) {
    for (ESMd::InnerIterator it(dgf, k); it; ++it) {
      rf.insert(it.row(), it.col()) = boltzman * temperature / plank * std::exp(-1.0 * it.value() / (gas * temperature));
    }
  }
  for (int k = 0; k < dgb.outerSize(); ++k) {
    for (ESMd::InnerIterator it(dgb, k); it; ++it) {
      rb.insert(it.row(), it.col()) = boltzman * temperature / plank * std::exp(-1.0 * it.value() / (gas * temperature));
    }
  }

  // Build and return Network
  Eigen::VectorXd masses = Eigen::Map<Eigen::VectorXd>(vertexMasses.data(), vertexMasses.size());
  auto rateConstants = std::make_tuple<ESMd, ESMd>(std::move(rf), std::move(rb));
  auto stoichiometry = std::make_tuple<ESMi, ESMi>(std::move(sf), std::move(sb));
  return Network(masses, rateConstants, stoichiometry);
}

} /* namespace Kinetx */
} /* namespace Scine */
