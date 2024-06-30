#ifndef NURBS_CURVE_H
#define NURBS_CURVE_H

#include <vector>

template<class Point>
class NurbsCurve
{
	using PointArray = std::vector<Point>;

public:
	NurbsCurve();
	NurbsCurve(int degree, const std::vector<double>& knots, 
		const std::vector<Point>& ctrlPts, std::vector<double>& weights);
	~NurbsCurve();
	NurbsCurve(const NurbsCurve& other);
	NurbsCurve& operator = (const NurbsCurve& other);

	Point EvalPointByBasis(double t);

	Point EvalPointByDeBoor(double u);

	//求一阶，二阶导数
	PointArray CalDerivative(int order, double u);

	// 节点插入算法,插入一个节点，保持曲线形状不变
	PointArray KnotInsertion(double u, int num = 1);

	// 移除节点，保持曲线形状不变，注意，不一定能成功移除
	NurbsCurve<Point> KnotRemoval(double u, int times);

	NurbsCurve<Point> KnotRefine(std::vector<double>& us);

	void Decompose(std::vector<NurbsCurve<Point>>& bezierSegs) const;

	NurbsCurve<Point> NurbsCurve<Point>::degreeElevate(int t);

	bool BezDegreeReduce(NurbsCurve<Point>& result);

	bool BezDegreeReduce(const PointArray& ctrlPts, PointArray& ctrlPts_reduce);

	NurbsCurve<Point> DegreeReduce();

	void ReparamLinear(double low, double high, NurbsCurve<Point>& result);

	void TessellateEqualKnot(std::vector<Point>& tessPts, std::vector<double>& tessUs, int sampNums);

	double ParamOfPoint(const Point& pt);

private:
	struct DataRep;
	std::unique_ptr<DataRep> m_dataRep;
};

#endif