#ifndef BSPLINE_BASIS_HH
#define BSPLINE_BASIS_HH

#include "vector"

//the class that provide some algorithm of b-spline basis
//所有的计算都是基于节点向量计算
class BSplineBasis
{
	//二分查找计算
	//注意在区间[u(p),u(m-p)]都是最多要算p+1个基函数， 如果在i<p 或者i>m-p 则需要计算的基函数少于p+1个 
	//所以这个时候直接计算开头的p+1个基函数 和 结尾的 p+1个基函数 一定是对的
	static int FindSpan(const std::vector<double>& knots, double u)
	{

	}
};

#endif