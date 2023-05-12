/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#include "Kinetx/Network.h"
#include "Kinetx/RungeKutta/CashKarp5.h"
#include "Kinetx/RungeKutta/ImplicitEuler.h"
//#include "Kinetx/RandomNetworkFactory.h"
#include "Kinetx/ReferenceNetworks.h"
#include <gmock/gmock.h>

using namespace testing;
namespace Scine {
namespace Kinetx {
namespace Tests {

TEST(NetworkTest, Constructor) {
  const unsigned int nReactions = 4;
  const unsigned int nCompounds = 5;
  Eigen::VectorXd masses = Eigen::VectorXd::Constant(nCompounds, 1);
  Eigen::MatrixXd lhsRates(nReactions, 1);
  lhsRates << 0.1, 0.05, 0.02, 0.02;
  Eigen::MatrixXd rhsRates(nReactions, 1);
  rhsRates << 0.05, 0.05, 0.001, 0.0000001;
  Eigen::MatrixXi stoichiometryForward = Eigen::MatrixXi::Zero(nReactions, nCompounds);
  stoichiometryForward(0, 0) = 1;
  stoichiometryForward(1, 1) = 1;
  stoichiometryForward(2, 2) = 1;
  stoichiometryForward(2, 3) = 1;
  stoichiometryForward(3, 0) = 1;
  stoichiometryForward(3, 1) = 1;
  Eigen::MatrixXi stoichiometryBackward = Eigen::MatrixXi::Zero(nReactions, nCompounds);
  stoichiometryBackward(0, 2) = 1;
  stoichiometryBackward(1, 3) = 1;
  stoichiometryBackward(2, 0) = 1;
  stoichiometryBackward(2, 1) = 1;
  stoichiometryBackward(3, 4) = 1;
  std::vector<std::string> labels = {"c1", "c2", "c3", "c4", "c5"};
  ASSERT_NO_THROW(Network network(masses, {lhsRates.sparseView(), rhsRates.sparseView()},
                                  {stoichiometryForward.sparseView(), stoichiometryBackward.sparseView()}, labels));
  auto wrongStoichiometry = stoichiometryBackward;
  wrongStoichiometry.conservativeResize(nReactions, nCompounds + 1);
  auto wrongRateConstants = lhsRates;
  wrongRateConstants.conservativeResize(nReactions + 1, 1);
  ASSERT_THROW(Network network(masses, {wrongRateConstants.sparseView(), rhsRates.sparseView()},
                               {stoichiometryForward.sparseView(), stoichiometryBackward.sparseView()}, labels),
               std::runtime_error);
  ASSERT_THROW(Network network(masses, {lhsRates.sparseView(), wrongRateConstants.sparseView()},
                               {stoichiometryForward.sparseView(), stoichiometryBackward.sparseView()}, labels),
               std::runtime_error);
  Network network(masses, {lhsRates.sparseView(), rhsRates.sparseView()},
                  {stoichiometryForward.sparseView(), stoichiometryBackward.sparseView()}, labels);

  Eigen::VectorXd referenceConcentrations(labels.size());
  referenceConcentrations << 3.337327e-02, 5.991047e-05, 6.674503e-02, 5.839152e-05, 3.998817e-01;
  Eigen::VectorXd concentrations(labels.size());
  concentrations << 0.5, 0.4, 0.0, 0.0, 0.0;

  ImplicitEuler solver(network);
  Eigen::VectorXd rFlux1, rForwardFlux1, rBackwardFlux1;
  const Eigen::MatrixXd results =
      solver.runIntegration(concentrations, 0.0, 7.75, rFlux1, rForwardFlux1, rBackwardFlux1, 1000, 100, 1e-10);
  const double error = (results.col(0) - referenceConcentrations).array().abs().maxCoeff();
  ASSERT_NEAR(error, 0.0, 1e-8);

  CashKarp5 solver2(network);
  Eigen::VectorXd rFlux2, rForwardFlux2, rBackwardFlux2;
  const Eigen::MatrixXd results2 =
      solver2.runIntegration(concentrations, 0.0, 7.75, rFlux2, rForwardFlux2, rBackwardFlux2, 1000, 100, 1e-10);
  const double error2 = (results2.col(0) - referenceConcentrations).array().abs().maxCoeff();
  ASSERT_NEAR(error2, 0.0, 1e-8);
  Eigen::VectorXd rFluxRef(4);
  rFluxRef << 9.422121e-01, 7.592755e-01, 6.809862e-01, 4.153817e-01;
  double errorFlux = (rFlux1.array() - rFluxRef.array()).abs().maxCoeff();
  ASSERT_NEAR(errorFlux, 0.0, 0.2);
  errorFlux = (rFlux1.array() - rFlux2.array()).abs().maxCoeff();
  ASSERT_NEAR(errorFlux, 0.0, 0.2);
  Eigen::VectorXd rForwardFluxRef(4);
  rForwardFluxRef << 2.487043e+02, 1.206172e+01, 7.018483e-01, 4.181250e-01;
  double errorFluxForward = (rForwardFlux1.array() - rForwardFluxRef.array()).abs().maxCoeff();
  ASSERT_NEAR(errorFluxForward, 0.0, 0.2);
  errorFluxForward = (rForwardFlux1.array() - rForwardFlux2.array()).abs().maxCoeff();
  ASSERT_NEAR(errorFluxForward, 0.0, 130);
  Eigen::VectorXd rBackwardFluxRef(4);
  rBackwardFluxRef << 2.477620e+02, 1.130244e+01, 2.090625e-02, 2.743272e-03;
  double errorFluxBackward = (rBackwardFlux1.array() - rBackwardFluxRef.array()).abs().maxCoeff();
  ASSERT_NEAR(errorFluxBackward, 0.0, 0.2);
  errorFluxBackward = (rBackwardFlux1.array() - rBackwardFlux2.array()).abs().maxCoeff();
  ASSERT_NEAR(errorFluxBackward, 0.0, 130);
  double additivity = ((rForwardFlux1 - rBackwardFlux1).cwiseAbs() - rFlux1).array().abs().maxCoeff();
  ASSERT_NEAR(additivity, 0.0, 1e-2);
  Eigen::VectorXd rFlux3, rForwardFlux3, rBackwardFlux3;
  const Eigen::MatrixXd results3 =
      solver2.runIntegrationByTime(concentrations, 0.0, 7.75, rFlux2, rForwardFlux2, rBackwardFlux2, 1e+8, 1000, 1e-10);
  const double error3 = (results2 - results3).array().abs().sum();
  ASSERT_NEAR(error3, 0.0, 1e-9);
}

// TODO Random network generation appears to be broken.
// TEST(NetworkTest, IntegrateRandomNetwork) {
//  Network randomNetwork = RandomNetworkFactory::random();
//  Eigen::VectorXd concentrations(randomNetwork.nCompounds);
//  concentrations(0) = 0.5;
//  concentrations(1) = 0.4;
//  concentrations(2) = 0.3;
//  CashKarp5 solver(randomNetwork);
//  const Eigen::MatrixXd results = solver.runIntegration(concentrations, 0.0, 1e-6, 1000, 5, 1e-10);
//}

TEST(NetworkTest, IntegrateExampleNetwork) {
  auto exampleNetworkAndConcentrations = ReferenceNetworks::getBrayLiebhafsky();
  const Eigen::VectorXd& concentrations = exampleNetworkAndConcentrations.second;
  CashKarp5 solver(exampleNetworkAndConcentrations.first);
  Eigen::VectorXd rFlux, rForwardFlux, rBackwardFlux;
  const Eigen::MatrixXd results =
      solver.runIntegration(concentrations, 0.0, 1e-6, rFlux, rForwardFlux, rBackwardFlux, 1000, 5, 1e-10);
}

TEST(NetworkTest, IntegrateTiming) {
  const unsigned int nReactions = 135;
  Eigen::VectorXd lhsRateConstants(nReactions);
  lhsRateConstants << 62124359349210.9, 27.7628983943279, 12435.2725043387, 12423.4954446288, 840263563.210643,
      77613857660.7285, 2.49276338695478e-06, 271264074.481528, 39770639163.7788, 0.000831048310776266,
      0.00265721801163604, 294.537834463143, 5.51133284417578e-16, 2.94664670525467, 3.22140548279155, 1186254442.17973,
      3.80509308137194e-06, 2.67127939065399e-06, 3.06151114507626e-07, 3.37980273357383e-07, 22923721013.9432,
      85334099.3280445, 0.126143889839647, 3.6997785602884e-06, 64097633.5738327, 757529.506694868, 0.414418776412466,
      165.954848631939, 37274615609526.5, 338.660933934075, 321685.42528564, 138058.321375682, 1.30094635845809e-08,
      197948256182.167, 24692.2962193414, 4587820.27531595, 1805339.67462677, 122811697317901, 3.46903162381707e-13,
      197.674898199635, 0.000220853350475281, 765.546684273683, 31494544.4707981, 9.19611020975147e-05, 222.669110475415,
      1.46574018297296e-11, 0.156232973573933, 0.0411614090623102, 1.25583722518024e-18, 180160642112711,
      23.7365833916396, 665.627077380831, 8.32681779151497e-06, 210.534025724134, 404.958713852586, 134.568399831855,
      743.736007343142, 777.728191890678, 554889.796776671, 722697.841584704, 402489.630306367, 0.987187432838048,
      33336.9848107515, 1336150729229.63, 669.661503480084, 0.01783217821359, 694269.32955873, 2111.99671150532,
      144.328328148146, 1.82229480833925, 8477.38592778171, 583405.399002643, 33035.3365483423, 6.40942756430428,
      474547.497469992, 3.39319350911653e-06, 65.4562589689802, 87722640941357.5, 0.000293381662741161, 7.07248886143678,
      1.15040175372236e-07, 5.6308077606605e-05, 105785.701850385, 7.92315881639237e-10, 0.0114216454749059,
      19314.39786267, 1433390784.30051, 2.01417179759802, 7.86809939035748e-10, 0.000475407569663483, 0.807880911796242,
      26.0529995866734, 0.0760635666682921, 4.78712543091802e-16, 6.10082256268928e-08, 35205.7476816651,
      15109.9115626323, 7.08715843393486e-11, 1.578812852098e-07, 11811842.0393824, 16.6956132646489,
      4.96268477187418e-06, 200.647171772802, 152639.717877085, 2.45840218255256e-06, 3.30727087971018e-05,
      2.24299628652706e-16, 1982.94541710999, 2351168802.59594, 9.2751874147162e-09, 4.60457751456594e-07,
      1.981087722071e-14, 12.2227184647711, 2.05634012428369e-17, 544357704500.959, 2671.05023664789,
      2.72137663026263e-07, 1282361.12774153, 1.32494122485328e-07, 8.87026454008282e-06, 2672.58380885067,
      9.04391104024767e-06, 1.00034147840176e-07, 29240880313785.8, 0.499768258609831, 7.24950781858782e-09,
      8.49324540989336e-06, 419268962704.965, 3.09817777222851e-12, 1.56412351543756e-09, 0.00211571098471571,
      8.7669456102916e-09, 3.01859243693822e-07, 0.169738543119018, 0.00506794971599136;
  Eigen::VectorXd rhsRateConstants(nReactions);
  rhsRateConstants << 32935919088980, 1.10865903592661, 582.845552861574, 4.74458646040522, 6627944897015.12,
      6212435934921.08, 308.578126024709, 35869098.1788939, 767671.843332567, 6418865507345.88, 1200.975440871,
      4778986.55699372, 6135.85382701471, 8028959566.46065, 4019730213.33723, 21604455579904.2, 34.1051584502082,
      114.903563362347, 137074.211234799, 586419.79512364, 856732442083.612, 1040573.25043198, 159985069.824359,
      12258.5738231834, 6212435934921.08, 5.79666433402401e-07, 6212435934921.08, 6212435934921.08, 6619692128184.21,
      13.5647604927723, 0.120740964925733, 0.0673093774881122, 6212435934921.08, 8632483448081.17, 8.37853769152033,
      0.0490036372957926, 0.0278898085793505, 15237502481075.7, 6212435934921.08, 37196522896.024, 14727736.3805283,
      1944.61832867037, 6557223602.66771, 6212435934921.08, 6212435934921.08, 55.0608599960409, 129101316525.194,
      21622913744.3301, 39.6942723897158, 46720819951313.9, 6310174493093.13, 6.65642571866675e-05, 12409108556664.5,
      695323043833.616, 4.19895208741217e-13, 12356959858188.4, 1.16953497998633e-06, 4.60685188397145e-07,
      0.00432786535627985, 2.90808976381034e-10, 2.01706948421525e-07, 6212435934921.08, 3.40879421861125e-09,
      29991784444600.4, 7.1917889092062e-05, 301625792283.409, 2.64796224969178e-10, 2.07187662103155e-12,
      600257062272.051, 8257810370640.58, 5.24935901030364e-06, 0.00441484336855049, 3.75751728157174e-09,
      1710389143500.62, 2.24559207703816e-07, 13905.5733788499, 6212435934921.08, 13667463799058.9, 65757.7970043802,
      1747139857918.85, 337.614417787251, 17667699.5608621, 6222417871104.33, 3081355791465.47, 1536662774.11615,
      6212435934921.08, 17161935724193.9, 196812209368.595, 22.3103643602754, 1217.20573619071, 266728.665153158,
      0.253851472306758, 24847524.5764019, 1264360716.72874, 158230542954.103, 14683731.3225577, 6212435934921.08,
      13.0929767452944, 0.184651673172134, 17930328859257.5, 70967.0718985688, 0.30403583848242, 4620.36548615971,
      90517798.4878715, 20586.7607521191, 499027277046.417, 176999848841.922, 210617290535.108, 7331302557061.88,
      34786053932.9041, 293874226.904294, 6212435934921.08, 6237614204633.88, 6212435934921.08, 13736151865044.8,
      1312214480698.43, 3014274.84096, 4109205824.26303, 12.2403724349562, 27906.935964984, 267183346419.141,
      26779.4103372219, 13.1077131278882, 12707348976171, 6212435934921.08, 1552.55012900126, 12213011.3952068,
      6722442605724.07, 0.0723199665333495, 5004292246953.64, 1058055752772.09, 183761.395296694, 1222.84555616168,
      412885256.870685, 2445026939004.22;
  const unsigned int nCompounds = 138;
  Eigen::VectorXd startConcentrations = Eigen::VectorXd::Zero(nCompounds);
  startConcentrations(0) = 1.2;
  startConcentrations(1) = 1.0;
  startConcentrations(6) = 0.05;
  Eigen::MatrixXi forward = Eigen::MatrixXi::Zero(nReactions, nCompounds);
  Eigen::MatrixXi backward = Eigen::MatrixXi::Zero(nReactions, nCompounds);
  forward(0, 0) += 1;
  forward(0, 1) += 1;
  backward(0, 2) += 1;
  forward(1, 2) += 1;
  backward(1, 3) += 1;
  forward(2, 2) += 1;
  backward(2, 4) += 1;
  forward(3, 2) += 1;
  backward(3, 5) += 1;
  forward(4, 0) += 1;
  forward(4, 6) += 1;
  backward(4, 7) += 1;
  forward(5, 7) += 1;
  backward(5, 8) += 1;
  forward(6, 7) += 1;
  backward(6, 9) += 1;
  forward(7, 7) += 1;
  backward(7, 10) += 1;
  forward(8, 7) += 1;
  backward(8, 11) += 1;
  forward(9, 7) += 1;
  backward(9, 12) += 1;
  forward(10, 7) += 1;
  backward(10, 13) += 1;
  forward(11, 7) += 1;
  backward(11, 14) += 1;
  forward(12, 7) += 1;
  backward(12, 15) += 1;
  forward(13, 7) += 1;
  backward(13, 6) += 1;
  backward(13, 16) += 1;
  forward(14, 7) += 1;
  backward(14, 6) += 1;
  backward(14, 17) += 1;
  forward(15, 1) += 1;
  forward(15, 6) += 1;
  backward(15, 18) += 1;
  forward(16, 18) += 1;
  backward(16, 19) += 1;
  forward(17, 18) += 1;
  backward(17, 20) += 1;
  forward(18, 18) += 1;
  backward(18, 21) += 1;
  forward(19, 18) += 1;
  backward(19, 22) += 1;
  forward(20, 18) += 1;
  backward(20, 23) += 1;
  forward(21, 18) += 1;
  backward(21, 24) += 1;
  forward(22, 18) += 1;
  backward(22, 25) += 1;
  forward(23, 18) += 1;
  backward(23, 26) += 1;
  forward(24, 18) += 1;
  backward(24, 27) += 1;
  backward(24, 28) += 1;
  forward(25, 18) += 1;
  backward(25, 29) += 1;
  forward(26, 18) += 1;
  backward(26, 30) += 1;
  forward(27, 18) += 1;
  backward(27, 31) += 1;
  forward(28, 0) += 1;
  forward(28, 4) += 1;
  backward(28, 32) += 1;
  forward(29, 32) += 1;
  backward(29, 33) += 1;
  forward(30, 32) += 1;
  backward(30, 34) += 1;
  forward(31, 32) += 1;
  backward(31, 35) += 1;
  forward(32, 32) += 1;
  backward(32, 36) += 1;
  forward(33, 0) += 1;
  forward(33, 5) += 1;
  backward(33, 37) += 1;
  forward(34, 37) += 1;
  backward(34, 38) += 1;
  forward(35, 37) += 1;
  backward(35, 39) += 1;
  forward(36, 37) += 1;
  backward(36, 40) += 1;
  forward(37, 0) += 1;
  forward(37, 11) += 1;
  backward(37, 41) += 1;
  forward(38, 41) += 1;
  backward(38, 42) += 1;
  forward(39, 41) += 1;
  backward(39, 43) += 1;
  forward(40, 41) += 1;
  backward(40, 44) += 1;
  backward(40, 45) += 1;
  forward(41, 41) += 1;
  backward(41, 46) += 1;
  forward(42, 41) += 1;
  backward(42, 47) += 1;
  forward(43, 41) += 1;
  backward(43, 48) += 1;
  forward(44, 41) += 1;
  backward(44, 49) += 1;
  forward(45, 41) += 1;
  backward(45, 50) += 1;
  forward(46, 41) += 1;
  backward(46, 11) += 1;
  backward(46, 16) += 1;
  forward(47, 41) += 1;
  backward(47, 11) += 1;
  backward(47, 17) += 1;
  forward(48, 41) += 1;
  backward(48, 51) += 1;
  forward(49, 1) += 1;
  forward(49, 4) += 1;
  backward(49, 52) += 1;
  forward(50, 52) += 1;
  backward(50, 53) += 1;
  forward(51, 52) += 1;
  backward(51, 54) += 1;
  forward(52, 52) += 1;
  backward(52, 55) += 1;
  forward(53, 52) += 1;
  backward(53, 56) += 1;
  forward(54, 52) += 1;
  backward(54, 57) += 1;
  forward(55, 52) += 1;
  backward(55, 58) += 1;
  forward(56, 52) += 1;
  backward(56, 59) += 1;
  forward(57, 52) += 1;
  backward(57, 60) += 1;
  forward(58, 52) += 1;
  backward(58, 61) += 1;
  forward(59, 52) += 1;
  backward(59, 62) += 1;
  forward(60, 52) += 1;
  backward(60, 63) += 1;
  forward(61, 52) += 1;
  backward(61, 64) += 1;
  forward(62, 52) += 1;
  backward(62, 65) += 1;
  forward(63, 1) += 1;
  forward(63, 5) += 1;
  backward(63, 66) += 1;
  forward(64, 66) += 1;
  backward(64, 67) += 1;
  forward(65, 66) += 1;
  backward(65, 68) += 1;
  forward(66, 66) += 1;
  backward(66, 69) += 1;
  forward(67, 66) += 1;
  backward(67, 70) += 1;
  forward(68, 66) += 1;
  backward(68, 71) += 1;
  forward(69, 66) += 1;
  backward(69, 72) += 1;
  forward(70, 66) += 1;
  backward(70, 73) += 1;
  forward(71, 66) += 1;
  backward(71, 74) += 1;
  forward(72, 66) += 1;
  backward(72, 75) += 1;
  forward(73, 66) += 1;
  backward(73, 76) += 1;
  forward(74, 66) += 1;
  backward(74, 77) += 1;
  forward(75, 66) += 1;
  backward(75, 78) += 1;
  forward(76, 66) += 1;
  backward(76, 79) += 1;
  forward(77, 1) += 1;
  forward(77, 11) += 1;
  backward(77, 80) += 1;
  forward(78, 80) += 1;
  backward(78, 81) += 1;
  forward(79, 80) += 1;
  backward(79, 82) += 1;
  forward(80, 80) += 1;
  backward(80, 83) += 1;
  forward(81, 80) += 1;
  backward(81, 84) += 1;
  forward(82, 80) += 1;
  backward(82, 85) += 1;
  forward(83, 80) += 1;
  backward(83, 86) += 1;
  forward(84, 80) += 1;
  backward(84, 87) += 1;
  forward(85, 80) += 1;
  backward(85, 88) += 1;
  forward(86, 6) += 1;
  forward(86, 4) += 1;
  backward(86, 89) += 1;
  forward(87, 89) += 1;
  backward(87, 90) += 1;
  forward(88, 89) += 1;
  backward(88, 91) += 1;
  forward(89, 89) += 1;
  backward(89, 92) += 1;
  forward(90, 89) += 1;
  backward(90, 93) += 1;
  forward(91, 89) += 1;
  backward(91, 94) += 1;
  forward(92, 89) += 1;
  backward(92, 95) += 1;
  forward(93, 89) += 1;
  backward(93, 96) += 1;
  forward(94, 89) += 1;
  backward(94, 97) += 1;
  forward(95, 89) += 1;
  backward(95, 98) += 1;
  forward(96, 89) += 1;
  backward(96, 99) += 1;
  forward(97, 89) += 1;
  backward(97, 100) += 1;
  forward(98, 89) += 1;
  backward(98, 101) += 1;
  forward(99, 6) += 1;
  forward(99, 5) += 1;
  backward(99, 102) += 1;
  forward(100, 102) += 1;
  backward(100, 103) += 1;
  forward(101, 102) += 1;
  backward(101, 104) += 1;
  forward(102, 102) += 1;
  backward(102, 105) += 1;
  forward(103, 102) += 1;
  backward(103, 106) += 1;
  forward(104, 102) += 1;
  backward(104, 107) += 1;
  forward(105, 102) += 1;
  backward(105, 108) += 1;
  forward(106, 102) += 1;
  backward(106, 109) += 1;
  forward(107, 102) += 1;
  backward(107, 110) += 1;
  forward(108, 6) += 1;
  forward(108, 11) += 1;
  backward(108, 111) += 1;
  forward(109, 111) += 1;
  backward(109, 112) += 1;
  forward(110, 111) += 1;
  backward(110, 113) += 1;
  forward(111, 111) += 1;
  backward(111, 114) += 1;
  backward(111, 45) += 1;
  forward(112, 111) += 1;
  backward(112, 115) += 1;
  forward(113, 111) += 1;
  backward(113, 116) += 1;
  backward(113, 117) += 1;
  forward(114, 4) += 1;
  forward(114, 5) += 1;
  backward(114, 118) += 1;
  forward(115, 118) += 1;
  backward(115, 119) += 1;
  forward(116, 111) += 1;
  backward(116, 6) += 1;
  backward(116, 13) += 1;
  forward(117, 118) += 1;
  backward(117, 120) += 1;
  forward(118, 118) += 1;
  backward(118, 121) += 1;
  forward(119, 118) += 1;
  backward(119, 122) += 1;
  forward(120, 118) += 1;
  backward(120, 123) += 1;
  forward(121, 118) += 1;
  backward(121, 124) += 1;
  forward(122, 118) += 1;
  backward(122, 125) += 1;
  forward(123, 4) += 1;
  forward(123, 11) += 1;
  backward(123, 126) += 1;
  forward(124, 126) += 1;
  backward(124, 127) += 1;
  forward(125, 126) += 1;
  backward(125, 128) += 1;
  forward(126, 126) += 1;
  backward(126, 129) += 1;
  forward(127, 5) += 1;
  forward(127, 11) += 1;
  backward(127, 130) += 1;
  forward(128, 130) += 1;
  backward(128, 131) += 1;
  forward(129, 130) += 1;
  backward(129, 11) += 1;
  backward(129, 132) += 1;
  forward(130, 130) += 1;
  backward(130, 133) += 1;
  forward(131, 130) += 1;
  backward(131, 134) += 1;
  forward(132, 130) += 1;
  backward(132, 135) += 1;
  forward(133, 130) += 1;
  backward(133, 136) += 1;
  forward(134, 130) += 1;
  backward(134, 137) += 1;
  backward(134, 5) += 1;
  std::vector<std::string> labels;
  for (unsigned int i = 0; i < nCompounds; ++i) {
    labels.push_back(std::to_string(i));
  }
  Eigen::VectorXd masses = Eigen::VectorXd::Constant(nCompounds, 1);
  Network network(masses, {lhsRateConstants.sparseView(), rhsRateConstants.sparseView()},
                  {forward.sparseView(), backward.sparseView()}, labels);
  CashKarp5 solver(network);
  Eigen::VectorXd rFlux, rForwardFlux, rBackwardFlux;
  const Eigen::MatrixXd results =
      solver.runIntegration(startConcentrations, 0.0, 1e-8, rFlux, rForwardFlux, rBackwardFlux, 100, 5, 1e-10);
}

TEST(NetworkTest, NumericalStability) {
  const unsigned int nReactions = 24;
  Eigen::VectorXd lhsRateConstants(nReactions);
  lhsRateConstants << 1147444934826.03, 0.000584295479706238, 0.287691621543586, 10833677500.3147, 0.000785238153654118,
      1467146988.41555, 144577781292854, 59371337.2532026, 12906812456.8549, 0.1134011882258, 2526063918.35424,
      41657.7983794888, 0.00170166658943427, 8.45683649867227e-06, 0.000455046144157151, 5.76988572273214e-09,
      811.912974344609, 37900.1997021746, 91327266.0697536, 1513.22433395654, 9179541669387.54, 603382.465629377,
      16064197921428.2, 6884656252040.66;
  Eigen::VectorXd rhsRateConstants(nReactions);
  rhsRateConstants << 2397745253112.45, 3139940077.22654, 2.36001560997978, 1147442708673.44, 28.6622869714074,
      1147442708673.44, 38285343965887.5, 498043.275363807, 56224692724998.7, 4589771522347.86, 35759.2146051663,
      331554130805.803, 113228191.493695, 2880497607532.27, 1147468094790.28, 18838951622.1974, 1147442708673.44,
      162386962329.182, 78122652621.6779, 138.472765594363, 9179541669387.54, 9.13744651937366e-05, 16064197921428.2,
      16287953394.2778;
  const unsigned int nCompounds = 25;
  Eigen::VectorXd startConcentrations = Eigen::VectorXd::Zero(nCompounds);
  startConcentrations(0) = 1.2;
  startConcentrations(1) = 0.05;
  Eigen::MatrixXi forward = Eigen::MatrixXi::Zero(nReactions, nCompounds);
  Eigen::MatrixXi backward = Eigen::MatrixXi::Zero(nReactions, nCompounds);
  forward(0, 0) += 1;
  forward(0, 1) += 1;
  backward(0, 2) += 1;
  forward(1, 2) += 1;
  backward(1, 1) += 1;
  backward(1, 3) += 1;
  forward(2, 2) += 1;
  backward(2, 4) += 1;
  forward(3, 2) += 1;
  backward(3, 5) += 1;
  forward(4, 2) += 1;
  backward(4, 6) += 1;
  forward(5, 2) += 1;
  backward(5, 7) += 1;
  forward(6, 0) += 1;
  forward(6, 5) += 1;
  backward(6, 8) += 1;
  forward(7, 8) += 1;
  backward(7, 9) += 1;
  forward(8, 8) += 1;
  backward(8, 0) += 1;
  backward(8, 0) += 1;
  backward(8, 1) += 1;
  forward(9, 8) += 1;
  backward(9, 10) += 1;
  forward(10, 8) += 1;
  backward(10, 11) += 1;
  forward(11, 8) += 1;
  backward(11, 13) += 1;
  forward(12, 8) += 1;
  backward(12, 15) += 1;
  backward(12, 16) += 1;
  forward(13, 8) += 1;
  backward(13, 17) += 1;
  forward(14, 8) += 1;
  backward(14, 18) += 1;
  forward(15, 8) += 1;
  backward(15, 19) += 1;
  forward(16, 8) += 1;
  backward(16, 20) += 1;
  forward(17, 8) += 1;
  backward(17, 21) += 1;
  forward(18, 8) += 1;
  backward(18, 23) += 1;
  forward(19, 5) += 1;
  backward(19, 12) += 1;
  forward(20, 5) += 1;
  backward(20, 0) += 1;
  backward(20, 1) += 1;
  forward(21, 7) += 1;
  backward(21, 24) += 1;
  forward(22, 7) += 1;
  backward(22, 0) += 1;
  backward(22, 1) += 1;
  forward(23, 4) += 1;
  backward(23, 24) += 1;

  std::vector<std::string> labels;
  for (unsigned int i = 0; i < nCompounds; ++i) {
    labels.push_back(std::to_string(i));
  }
  Eigen::VectorXd masses = Eigen::VectorXd::Constant(nCompounds, 1);
  Network network(masses, {lhsRateConstants.sparseView(), rhsRateConstants.sparseView()},
                  {forward.sparseView(), backward.sparseView()}, labels);
  CashKarp5 solver(network);
  Eigen::VectorXd rFlux, rForwardFlux, rBackwardFlux;
  const Eigen::MatrixXd results =
      solver.runIntegration(startConcentrations, 0.0, 1e-12, rFlux, rForwardFlux, rBackwardFlux, 100, 100, 1e-10);
  double additivity = ((rForwardFlux - rBackwardFlux).cwiseAbs() - rFlux).array().abs().sum();
  ASSERT_NEAR(additivity, 0.0, 5e-1);
}

} /* namespace Tests */
} /* namespace Kinetx */
} /* namespace Scine */
