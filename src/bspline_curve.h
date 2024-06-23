#ifndef BSPLINE_CURVE_HH
#define BSPLINE_CURVE_HH

#include "bspline_basis.h"
#include <vector>
#include <memory>

//BSpline curve was define by the degree, knots and control points
/*
* 1. ��������ֵ(ֱ�ӷ� �� DeBoor��)
* 2. �������㣬֧��һ�׺Ͷ���
* 3. �ڵ����
* 4. ���ף�����
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

	//ͨ�������������, ����t:[ui,ui+1),��N(k,p) k = i, i-1, i-2, ... i-p 
	//��p+1��������Ϊ�����
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
	//����Beizier��deCas�㷨����ͣ�Ľ���cut cornor������deCas�������Կ�����DeBoor����������
	//Deboor�������и������ÿ�������Ǳ仯�ģ�deCas���䣬��������Beizie�Ľڵ��������ظ��ģ����������
	//�и�ϵ����һ����
	Point EvalPointByDeBoor(double u);

	//��һ�ף����׵���
	PointArray CalDerivative(int order, double u);

	// �ڵ�����㷨,����һ���ڵ㣬����������״����
	PointArray KnotInsertion(double u, int num = 1);

	// �Ƴ��ڵ㣬����������״���䣬ע�⣬��һ���ܳɹ��Ƴ�
	BSplineCurve<Point> KnotRemoval(double u, int times);

	BSplineCurve<Point> KnotRefine(std::vector<double>& us);

	void Decompose(std::vector<BSplineCurve<Point>>& bezierSegs) const;

	BSplineCurve<Point> BSplineCurve<Point>::degreeElevate(int t);
	
private:
	struct DataRep;
	std::unique_ptr<DataRep> m_dataRep;
};

#endif
