/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Kinetx/Network.h"
#include "Kinetx/NetworkBuilder.h"
#include "gmock/gmock.h"

using namespace testing;
namespace Scine {
namespace Kinetx {
namespace Tests {

TEST(NetworkBuilderTest, Constructor) {
  NetworkBuilder builder();
  ASSERT_TRUE(true);
}

TEST(NetworkBuilderTest, constructNetwork) {
  NetworkBuilder networkBuilder;
  unsigned int nCompounds = 5;
  unsigned int nReactions = 4;
  unsigned int nChannelsPerReaction = 1;
  networkBuilder.reserve(nCompounds - 1, nReactions - 1, nChannelsPerReaction);
  Eigen::VectorXd concentrations = Eigen::VectorXd::Zero(5);
  concentrations << 0.5, 0.4, 0.0, 0.0, 0.0;
  for (unsigned int i = 0; i < nCompounds; ++i) {
    networkBuilder.addCompound(1.0);
  } // for i
  Eigen::VectorXd rateConstantsForward_ref(4);
  rateConstantsForward_ref << 0.1, 0.05, 0.02, 0.02;
  Eigen::VectorXd rateConstantsBackward_ref(4);
  rateConstantsBackward_ref << 0.05, 0.05, 0.001, 0.0000001;
  networkBuilder.addReaction({rateConstantsForward_ref(0)}, {rateConstantsBackward_ref(0)}, {{0, 1}}, {{2, 1}});
  networkBuilder.addReaction({rateConstantsForward_ref(1)}, {rateConstantsBackward_ref(1)}, {{1, 1}}, {{3, 1}});
  networkBuilder.addReaction({rateConstantsForward_ref(2)}, {rateConstantsBackward_ref(2)}, {{2, 1}, {3, 1}},
                             {{0, 1}, {1, 1}});
  networkBuilder.addReaction({rateConstantsForward_ref(3)}, {rateConstantsBackward_ref(3)}, {{0, 1}, {1, 1}}, {{4, 1}});
  Network network = networkBuilder.generate();
  Eigen::MatrixXd rateConstants_forward = Eigen::MatrixXd(std::get<0>(network.rateConstants));
  Eigen::MatrixXd rateConstants_backward = Eigen::MatrixXd(std::get<1>(network.rateConstants));
  ASSERT_EQ(network.nCompounds, nCompounds);
  ASSERT_EQ(network.nReactions, nReactions);
  ASSERT_NEAR((rateConstantsForward_ref - rateConstants_forward).array().abs().sum(), 0.0, 1e-14);
  ASSERT_NEAR((rateConstantsBackward_ref - rateConstants_backward).array().abs().sum(), 0.0, 1e-14);
  ASSERT_THROW(networkBuilder.reserve(nCompounds - 1, nReactions, nChannelsPerReaction), std::runtime_error);
  ASSERT_THROW(networkBuilder.reserve(nCompounds, nReactions - 1, nChannelsPerReaction), std::runtime_error);
  ASSERT_THROW(networkBuilder.reserve(nCompounds, nReactions, 0), std::runtime_error);
  ASSERT_THROW(networkBuilder.addReaction({rateConstantsForward_ref(2), rateConstantsForward_ref(2)},
                                          {rateConstantsBackward_ref(2)}, {{2, 1}, {3, 1}}, {{0, 1}, {1, 1}}),
               std::runtime_error);
  ASSERT_THROW(networkBuilder.addReaction({rateConstantsForward_ref(2)}, {rateConstantsBackward_ref(2)},
                                          {{2, 1}, {3, 1}}, {{nCompounds, 1}, {1, 1}}),
               std::runtime_error);
}

} /* namespace Tests */
} /* namespace Kinetx */
} /* namespace Scine */
