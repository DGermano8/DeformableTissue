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
#include "projects/DeformableTissue/test/Test3DModel/TestMy3dBoxModel_MyBendingForce.hpp"

static Test3dBoxModel suite_Test3dBoxModel;

static CxxTest::List Tests_Test3dBoxModel = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_Test3dBoxModel( "projects/DeformableTissue/test/Test3DModel/TestMy3dBoxModel_MyBendingForce.hpp", 49, "Test3dBoxModel", suite_Test3dBoxModel, Tests_Test3dBoxModel );

static class TestDescription_Test3dBoxModel_TestSimpleCubeBending : public CxxTest::RealTestDescription {
public:
 TestDescription_Test3dBoxModel_TestSimpleCubeBending() : CxxTest::RealTestDescription( Tests_Test3dBoxModel, suiteDescription_Test3dBoxModel, 99, "TestSimpleCubeBending" ) {}
 void runTest() { suite_Test3dBoxModel.TestSimpleCubeBending(); }
} testDescription_Test3dBoxModel_TestSimpleCubeBending;

#include <cxxtest/Root.cpp>
