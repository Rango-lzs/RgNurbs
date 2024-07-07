#include "gtest/gtest.h"
#include "bezier_curve.h"

class BezierCurveTest : public testing::Test {

protected:

	virtual void SetUp()
	{
		std::vector<Point3d> ctrlPts{
		Point3d(0,0,0),
		Point3d(5,10,0),
		Point3d(15,10,0),
		Point3d(20,0,0)
		};
		m_curve = BeizierCurve<Point3d>(3, ctrlPts);
	};

	virtual void TearDown() {
	};

	void test_EvalPoint()
	{
		Point3d eval = m_curve.EvalPoint(0.3);
	}

private:
	BeizierCurve<Point3d> m_curve;
};


TEST_F(BezierCurveTest, EvalPoint) {
	test_EvalPoint();
}


//// 宏定义，用于注册测试用例
//#define REGISTER_TEST_CASE(test_case_name, test_name, test_func) \
//    TEST_F(test_case_name, test_name) { \
//        test_func(); \
//    }
//
//// 注册所有测试用例
//REGISTER_TEST_CASE(BezierCurveTest, EvalPoint, test_EvalPoint)
