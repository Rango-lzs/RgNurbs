#include "gtest/gtest.h"
#include "bspline_curve.h"


//// MyClassTest.cpp
//#include <gtest/gtest.h>
//#include "MyClass.h"
//
//// 测试类，继承自 ::testing::Test
//class MyClassTest : public ::testing::Test {
//protected:
//	void SetUp() override {
//		// 在每个测试之前设置初始状态
//		my_class_ = new MyClass(0);
//	}
//
//	void TearDown() override {
//		// 在每个测试之后清理资源
//		delete my_class_;
//		my_class_ = nullptr;
//	}
//
//	// 测试 GetValue 接口
//	void TestGetValue() {
//		EXPECT_EQ(my_class_->GetValue(), 0);
//		my_class_->SetValue(5);
//		EXPECT_EQ(my_class_->GetValue(), 5);
//	}
//
//	// 测试 SetValue 接口
//	void TestSetValue() {
//		my_class_->SetValue(10);
//		EXPECT_EQ(my_class_->GetValue(), 10);
//	}
//
//	// 测试 Increment 接口
//	void TestIncrement() {
//		EXPECT_EQ(my_class_->Increment(), 1);
//		EXPECT_EQ(my_class_->Increment(), 2);
//		my_class_->SetValue(10);
//		EXPECT_EQ(my_class_->Increment(), 11);
//	}
//
//	MyClass* my_class_;
//};
//
//// 定义测试用例，并调用测试类中的成员函数
//TEST_F(MyClassTest, GetValue) {
//	TestGetValue();
//}
//
//TEST_F(MyClassTest, SetValue) {
//	TestSetValue();
//}
//
//TEST_F(MyClassTest, Increment) {
//	TestIncrement();
//}