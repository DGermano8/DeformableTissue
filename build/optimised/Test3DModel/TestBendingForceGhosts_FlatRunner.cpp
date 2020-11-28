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
#include "projects/DeformableTissue/test/Test3DModel/TestBendingForceGhosts_Flat.hpp"

static Test3dBoxModel suite_Test3dBoxModel;

static CxxTest::List Tests_Test3dBoxModel = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_Test3dBoxModel( "projects/DeformableTissue/test/Test3DModel/TestBendingForceGhosts_Flat.hpp", 63, "Test3dBoxModel", suite_Test3dBoxModel, Tests_Test3dBoxModel );

static class TestDescription_Test3dBoxModel_TestCubeWithGhosts : public CxxTest::RealTestDescription {
public:
 TestDescription_Test3dBoxModel_TestCubeWithGhosts() : CxxTest::RealTestDescription( Tests_Test3dBoxModel, suiteDescription_Test3dBoxModel, 72, "TestCubeWithGhosts" ) {}
 void runTest() { suite_Test3dBoxModel.TestCubeWithGhosts(); }
} testDescription_Test3dBoxModel_TestCubeWithGhosts;

#include <cxxtest/Root.cpp>
