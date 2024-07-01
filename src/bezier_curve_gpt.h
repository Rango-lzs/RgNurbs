#pragma once

template <class T>
class BezierCurve_GPT {
public:
	using Point = T;
	BezierCurve_GPT(const std::vector<Point>& controlPoints) : controlPoints(controlPoints) {}

	// Calculate the point on the Bezier curve at parameter t
	Point calculate(double t) const {
		int n = controlPoints.size() - 1;
		Point result = { 0.0, 0.0, 0.0 };

		for (int i = 0; i <= n; ++i) {
			double coefficient = binomialCoefficient(n, i) * pow(1 - t, n - i) * pow(t, i);
			result = result + coefficient * controlPoints[i];
		}

		return result;
	}

	// Calculate the first derivative of the Bezier curve at parameter t
	Point firstDerivative(double t) const {
		int n = controlPoints.size() - 1;
		Point result = { 0.0, 0.0 , 0.0 };

		for (int i = 0; i < n; ++i) {
			double coefficient = n * binomialCoefficient(n - 1, i) * pow(1 - t, (n - 1) - i) * pow(t, i);
			result = result + coefficient * (controlPoints[i + 1] + (-1.0) * controlPoints[i]);
		}

		return result;
	}

	// Calculate the second derivative of the Bezier curve at parameter t
	Point secondDerivative(double t) const {
		int n = controlPoints.size() - 1;
		Point result = { 0.0, 0.0 , 0.0 };

		for (int i = 0; i < n - 1; ++i) {
			double coefficient = n * (n - 1) * binomialCoefficient(n - 2, i) * pow(1 - t, (n - 2) - i) * pow(t, i);
			result = result + coefficient * (controlPoints[i + 2] + (-2.0) * controlPoints[i + 1] + controlPoints[i]);
		}

		return result;
	}

private:
	std::vector<Point> controlPoints;

	// Calculate the binomial coefficient "n choose k"
	int binomialCoefficient(int n, int k) const {
		if (k > n) return 0;
		if (k == 0 || k == n) return 1;
		k = std::min(k, n - k); // take advantage of symmetry
		int c = 1;
		for (int i = 0; i < k; ++i) {
			c = c * (n - i) / (i + 1);
		}
		return c;
	}
};