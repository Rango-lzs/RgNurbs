#ifndef BSPLINE_BASIS_HH
#define BSPLINE_BASIS_HH

#include "vector"

//the class that provide some algorithm of b-spline basis
//���еļ��㶼�ǻ��ڽڵ���������
class BSplineBasis
{
	//���ֲ��Ҽ���
	//ע��������[u(p),u(m-p)]�������Ҫ��p+1���������� �����i<p ����i>m-p ����Ҫ����Ļ���������p+1�� 
	//�������ʱ��ֱ�Ӽ��㿪ͷ��p+1�������� �� ��β�� p+1�������� һ���ǶԵ�
	static int FindSpan(const std::vector<double>& knots, double u)
	{

	}
};

#endif