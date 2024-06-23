#ifndef BSPLINE_CURVE_HH
#define BSPLINE_CURVE_HH

#include "bspline_basis.h"
#include <vector>
#include <memory>

//BSpline curve was define by the degree, knots and control points
/*
* 1. 参数点求值(直接法 和 DeBoor法)
* 2. 导数计算，支持一阶和二阶
* 3. 节点插入
* 4. 升阶，降阶
* 5. Subdivide
*/

template<class Point>
class BSplineCurve
{
	using PointArray = std::vector<Point>;

public:
	BSplineCurve();
	BSplineCurve(int degree, const std::vector<double>& knots, const std::vector<Point>& ctrlPts);
	~BSplineCurve();
	BSplineCurve(const BSplineCurve& other);

	//通过基函数计算点, 对于t:[ui,ui+1),仅N(k,p) k = i, i-1, i-2, ... i-p 
	//共p+1个基函数为非零的
	Point EvalPointByBasis(double t);

	//ref: https://en.wikipedia.org/wiki/De_Boor%27s_algorithm
	/*Evaluates S(x).
	Arguments
		-------- -
		span : Index of knot interval that contains x.
		u : Position.
		knots : Array of knot positions, needs to be padded as described above.
		c : Array of control points.
		p : Degree of B - spline.
	*/
	//类似Beizier的deCas算法，不停的进行cut cornor操作，deCas方法可以看做是DeBoor方法的特例
	//Deboor方法中切割比例在每个区间是变化的，deCas不变，这是由于Beizie的节点向量是重复的，计算出来的
	//切割系数是一样的
	Point EvalPointByDeBoor(double u);

	//求一阶，二阶导数
	PointArray CalDerivative(int order, double u);

	// 节点插入算法,插入一个节点，保持曲线形状不变
	PointArray KnotInsertion(double u, int num = 1);

	// 移除节点，保持曲线形状不变，注意，不一定能成功移除
	BSplineCurve<Point> KnotRemoval(double u, int times);

	BSplineCurve<Point> KnotRefine(std::vector<double>& us);

	void Decompose(std::vector<BSplineCurve<Point>>& bezierSegs) const;

	BSplineCurve<Point> BSplineCurve<Point>::degreeElevate(int t);
	
private:
	struct DataRep;
	std::unique_ptr<DataRep> m_dataRep;
};

#endif
