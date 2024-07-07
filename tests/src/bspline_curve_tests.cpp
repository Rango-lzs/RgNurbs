#include "gtest/gtest.h"
#include "bspline_curve.h"


//// MyClassTest.cpp
//#include <gtest/gtest.h>
//#include "MyClass.h"
//
//// �����࣬�̳��� ::testing::Test
//class MyClassTest : public ::testing::Test {
//protected:
//	void SetUp() override {
//		// ��ÿ������֮ǰ���ó�ʼ״̬
//		my_class_ = new MyClass(0);
//	}
//
//	void TearDown() override {
//		// ��ÿ������֮��������Դ
//		delete my_class_;
//		my_class_ = nullptr;
//	}
//
//	// ���� GetValue �ӿ�
//	void TestGetValue() {
//		EXPECT_EQ(my_class_->GetValue(), 0);
//		my_class_->SetValue(5);
//		EXPECT_EQ(my_class_->GetValue(), 5);
//	}
//
//	// ���� SetValue �ӿ�
//	void TestSetValue() {
//		my_class_->SetValue(10);
//		EXPECT_EQ(my_class_->GetValue(), 10);
//	}
//
//	// ���� Increment �ӿ�
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
//// ������������������ò������еĳ�Ա����
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