#ifndef NURBS_CURVE_H
#define NURBS_CURVE_H

#include <vector>
#include <memory>

template<class HPoint>
class NurbsCurve
{
	using Point = typename HPoint::Point;
	using HPointArray = std::vector<HPoint>;
	using PointArray = std::vector<Point>;

public:
	NurbsCurve();
	NurbsCurve(int degree, const std::vector<double>& knots,
		const std::vector<HPoint>& ctrlPts, std::vector<double>& weights);
	~NurbsCurve();
	NurbsCurve(const NurbsCurve& other);
	NurbsCurve& operator = (const NurbsCurve& other);

	Point EvalPoint(double t);

	Point EvalHPointByDeBoor(double u);

	//求一阶，二阶导数
	PointArray CalDerivative(int order, double u);

	// 节点插入算法,插入一个节点，保持曲线形状不变
	PointArray KnotInsertion(double u, int num = 1);

	// 移除节点，保持曲线形状不变，注意，不一定能成功移除
	NurbsCurve<HPoint> KnotRemoval(double u, int times);

	NurbsCurve<HPoint> KnotRefine(std::vector<double>& us);

	void Decompose(std::vector<NurbsCurve<HPoint>>& bezierSegs) const;

	NurbsCurve<HPoint> NurbsCurve<HPoint>::degreeElevate(int t);

	bool BezDegreeReduce(NurbsCurve<HPoint>& result);

	bool BezDegreeReduce(const HPointArray& ctrlPts, HPointArray& ctrlPts_reduce);

	NurbsCurve<HPoint> DegreeReduce();

	void ReparamLinear(double low, double high, NurbsCurve<HPoint>& result);

	void TessellateEqualKnot(std::vector<Point>& tessPts, std::vector<double>& tessUs, int sampNums);

	double ParamOfPoint(const Point& pt);

private:
	struct DataRep;
	std::unique_ptr<DataRep> m_dataRep;
};

#endif