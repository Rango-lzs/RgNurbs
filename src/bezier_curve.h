#ifndef BEIZIER_CURVE_HH
#define BEIZIER_CURVE_HH

#include "point.h"
#include <vector>
#include <memory>

template<class Point>
class BeizierCurve
{
	static_assert(std::is_same_v<Point, Point3d>||std::is_same_v<Point, Point2d>, "only support point3d and point2d");
public:
	//beizier curve was total defined by control points
	BeizierCurve(int degree, std::vector<Point> ctrlPts);

	BeizierCurve();

	BeizierCurve(const BeizierCurve& other);

	BeizierCurve& operator = (const BeizierCurve& other);

	~BeizierCurve();
    //
	Point EvalPoint(double param);

	Point EvalPointDirect(double param);

	//µÝÍÆÈý½ÇÐÎ
	Point EvalPointByDeCasteljau(double param);

	// De Casteljau's algorithm to compute a point on a Bezier curve
	Point deCasteljau(const std::vector<Point>& controlPoints, double t);

	// Compute the derivative up to second order using De Casteljau's algorithm
	Point EvalDerivation(int order, double t);

	bool SubDivision(double param, std::vector<Point>& left, std::vector<Point>& right);

	bool SubDivision_1(double param, std::vector<Point>& cv_left, std::vector<Point>& cv_right);

	bool DegreeElevation(std::vector<Point>& nweCtrlPts);

private:
	struct DataRep;
	std::unique_ptr<DataRep> m_dataRep;
};

#endif
