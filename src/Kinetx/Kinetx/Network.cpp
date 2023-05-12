/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Kinetx/Network.h"

namespace Scine {
namespace Kinetx {

Network::Network(Eigen::VectorXd masses, std::tuple<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> rateConstants,
                 std::tuple<Eigen::SparseMatrix<int>, Eigen::SparseMatrix<int>> stoichiometry, std::vector<std::string> labels)
  : nCompounds(std::get<0>(stoichiometry).cols()),
    nReactions(std::get<0>(stoichiometry).rows()),
    masses(masses),
    rateConstants(rateConstants),
    stoichiometry(stoichiometry),
    stoichiometryTransposed({std::get<0>(stoichiometry).transpose().eval(), std::get<1>(stoichiometry).transpose().eval()}),
    totalStoichiometry(std::get<1>(stoichiometry) - std::get<0>(stoichiometry)),
    addedStoichiometry(std::get<1>(stoichiometry) + std::get<0>(stoichiometry)),
    compoundLabels((labels.size() == nCompounds) ? labels : std::vector<std::string>(nCompounds, "")) {
  // Error checks
  if (std::get<0>(rateConstants).rows() != nReactions)
    throw std::runtime_error("Forward rate constants' dimensions did not match: to few or too many compounds.");
  if (std::get<1>(rateConstants).rows() != nReactions)
    throw std::runtime_error("Backward rate constants' dimensions did not match: to few or too many compounds.");
  if (std::get<1>(stoichiometry).cols() != nCompounds)
    throw std::runtime_error("Stoichiometry compound number is not consistent.");
  if (std::get<1>(stoichiometry).rows() != nReactions)
    throw std::runtime_error("Stoichiometry reaction number is not consistent.");
}

// void Network::load(std::string /*path*/) {
//}
// void Network::write(std::string path) {
//  const auto& sf = std::get<0>(this->stoichiometry);
//  const auto& sb = std::get<1>(this->stoichiometry);
//
//  FILE* pFile;
//  std::string edgeFile = path + ".edges.dat";
//  pFile = fopen(edgeFile.c_str(), "w");
//  for (int col = 0; col < sf.outerSize(); ++col) {
//    for (Eigen::SparseMatrix<int>::InnerIterator it(sf, col); it; ++it) {
//      fprintf(pFile, "%5d %5d\n", int(col + sf.rows()), int(it.row()));
//    }
//  }
//  for (int col = 0; col < sb.outerSize(); ++col) {
//    for (Eigen::SparseMatrix<int>::InnerIterator it(sb, col); it; ++it) {
//      fprintf(pFile, "%5d %5d\n", int(it.row()), int(col + sb.rows()));
//    }
//  }
//  fclose(pFile);
//
//  std::string vertFile = path + ".vertices.dat";
//  pFile = fopen(vertFile.c_str(), "w");
//  for (long int i = 0; i < sb.rows(); ++i) {
//    fprintf(pFile, "%5d ts\n", int(i));
//  }
//  for (long int i = 0; i < sb.cols(); ++i) {
//    fprintf(pFile, "%5d reactant\n", int(i + sf.rows()));
//  }
//  fclose(pFile);
//}

} /* namespace Kinetx */
} /* namespace Scine */
