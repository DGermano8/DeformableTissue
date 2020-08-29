/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#define _CXXTEST_HAVE_EH
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ErrorPrinter.h>

#include "CommandLineArguments.hpp"
int main( int argc, char *argv[] ) {
 CommandLineArguments::Instance()->p_argc = &argc;
 CommandLineArguments::Instance()->p_argv = &argv;
 return CxxTest::ErrorPrinter().run();
}
#include "projects/DeformableTissue/test/Test3DModel/TestMyBendingForceGhosts_Check.hpp"

static Test3dBoxModel suite_Test3dBoxModel;

static CxxTest::List Tests_Test3dBoxModel = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_Test3dBoxModel( "projects/DeformableTissue/test/Test3DModel/TestMyBendingForceGhosts_Check.hpp", 51, "Test3dBoxModel", suite_Test3dBoxModel, Tests_Test3dBoxModel );

static class TestDescription_Test3dBoxModel_TestPeriodicCubeWithGhosts_SmallRegion : public CxxTest::RealTestDescription {
public:
 TestDescription_Test3dBoxModel_TestPeriodicCubeWithGhosts_SmallRegion() : CxxTest::RealTestDescription( Tests_Test3dBoxModel, suiteDescription_Test3dBoxModel, 87, "TestPeriodicCubeWithGhosts_SmallRegion" ) {}
 void runTest() { suite_Test3dBoxModel.TestPeriodicCubeWithGhosts_SmallRegion(); }
} testDescription_Test3dBoxModel_TestPeriodicCubeWithGhosts_SmallRegion;

static class TestDescription_Test3dBoxModel_TestPeriodicCubeWithGhosts_WholeTissue : public CxxTest::RealTestDescription {
public:
 TestDescription_Test3dBoxModel_TestPeriodicCubeWithGhosts_WholeTissue() : CxxTest::RealTestDescription( Tests_Test3dBoxModel, suiteDescription_Test3dBoxModel, 365, "TestPeriodicCubeWithGhosts_WholeTissue" ) {}
 void runTest() { suite_Test3dBoxModel.TestPeriodicCubeWithGhosts_WholeTissue(); }
} testDescription_Test3dBoxModel_TestPeriodicCubeWithGhosts_WholeTissue;

static class TestDescription_Test3dBoxModel_TestPeriodicCubeWithGhosts_PoppedUpCell : public CxxTest::RealTestDescription {
public:
 TestDescription_Test3dBoxModel_TestPeriodicCubeWithGhosts_PoppedUpCell() : CxxTest::RealTestDescription( Tests_Test3dBoxModel, suiteDescription_Test3dBoxModel, 642, "TestPeriodicCubeWithGhosts_PoppedUpCell" ) {}
 void runTest() { suite_Test3dBoxModel.TestPeriodicCubeWithGhosts_PoppedUpCell(); }
} testDescription_Test3dBoxModel_TestPeriodicCubeWithGhosts_PoppedUpCell;

static class TestDescription_Test3dBoxModel_TestPeriodicCubeWithGhosts_PoppedDownCell : public CxxTest::RealTestDescription {
public:
 TestDescription_Test3dBoxModel_TestPeriodicCubeWithGhosts_PoppedDownCell() : CxxTest::RealTestDescription( Tests_Test3dBoxModel, suiteDescription_Test3dBoxModel, 918, "TestPeriodicCubeWithGhosts_PoppedDownCell" ) {}
 void runTest() { suite_Test3dBoxModel.TestPeriodicCubeWithGhosts_PoppedDownCell(); }
} testDescription_Test3dBoxModel_TestPeriodicCubeWithGhosts_PoppedDownCell;

#include <cxxtest/Root.cpp>
