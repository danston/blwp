// Copyright Dmitry Anisimov danston@ymail.com (c) 2017.

// README:
/*

    A convex combination of Wachspress and any other 2D closed-form barycentric coordinates.

    This class depends on:
    1. BarycentricCoordinatesR2.hpp
    2. SegmentCoordinatesR2.hpp
    3. VertexExpressionsR2.hpp
    4. VertexR2.hpp

*/
    
#ifndef GBC_BLENDEDWACHSPRESSR2_HPP
#define GBC_BLENDEDWACHSPRESSR2_HPP

// STL includes.
#include <vector>
#include <cassert>
#include <cmath>

// Local includes.
#include "../extra/VertexR2.hpp"
#include "../extra/BarycentricCoordinatesR2.hpp"

namespace gbc {

    // A convex combination of Wachspress and any other 2D closed-form barycentric coordinates
    // given by the type Coordinate and complied with the requirements of this type (see for example MeanValueR2.hpp)
    template<class Coordinate>
    class BlendedWachspressR2 : public BarycentricCoordinatesR2 {

    public:
        // Constructor.
        BlendedWachspressR2(const std::vector<VertexR2> &v, const double tol = 1.0e-10) : super(v, tol), _bc(v, tol) { }

        // Return name of the coordinate function.
        inline std::string name() const {
            return "BlendedWachspressR2";
        }

        // Function that computes coordinates b at a point p. This implementation is based on the following dissertation:
        // D. Anisimov. Analysis and new constructions of generalized barycentric coordinates in 2D. To appear, 2017 (see Chapter 4).
        // For the pseudocode see Apeendix A of the same dissertation.
        void compute(const VertexR2 &p, std::vector<double> &b) {

            b.clear();

            const size_t n = _v.size();
            b.resize(n, 0.0);

            // Boundary.
            if (computeBoundaryCoordinates(p, b)) return;

            // Interior.
            std::vector<double> wp_b;
            std::vector<double> bc_b;

            computeWP(p, wp_b);
            const double mu = computeBlendingFunction(p);

            _bc.compute(p, bc_b);

            assert(b.size() == wp_b.size() && wp_b.size() == bc_b.size());
            for (size_t i = 0; i < n; ++i) b[i] = mu * wp_b[i] + (1.0 - mu) * bc_b[i];
        }

        // Compute the coordinates at p using the internal storage from the VertexR2 class.
        inline void compute(VertexR2 &p) {
            compute(p, p.b());
        }

        // Compute coordinates bb at all points p in the vector.
        void compute(const std::vector<VertexR2> &p, std::vector<std::vector<double> > &bb) {

            const size_t numP = p.size();
            
            bb.resize(numP);
            for (size_t i = 0; i < numP; ++i) compute(p[i], bb[i]);
        }

        // Compute coordinates at all points p in the vector using the internal storage from the VertexR2 class.
        void compute(std::vector<VertexR2> &p) {

            const size_t numP = p.size();
            for (size_t i = 0; i < numP; ++i) compute(p[i], p[i].b());
        }

        // Implementation of the virtual function to compute all coordinates.
        inline void bc(std::vector<VertexR2> &p) {
            compute(p);
        }

    private:
        // Some typedefs.
        typedef BarycentricCoordinatesR2 super;

        // Templated coordinates.
        Coordinate _bc;

        // Compute blending function.
        double computeBlendingFunction(const VertexR2 &p) const {

            const double W = getW(p);

            const double dWx = getdWx(p);
            const double dWy = getdWy(p);

            double sigma_1 = dWx * dWx + dWy * dWy;
            double sigma_2 = 1.0;

            for (size_t i = 0; i < _v.size(); ++i) 
                sigma_2 *= (_v[i] - p).squaredLength();

            sigma_1 = (W * W * W * W) / (sigma_1 * sigma_1);
            const double sigma = sigma_1 + sigma_2;

            double mu = 0.0;
            if (fabs(sigma) > 0.0) mu = sigma_1 / sigma;

            return mu;
        }

        // Get Wachspress denominator W.
        inline double getW(const VertexR2 &p) const {

            std::vector<double> w;
            return computeWPWeights(p, w);
        }

        // Get x partial derivative dWx of the Wachspress denominator W.
        double getdWx(const VertexR2 &p) const {

            const int n = 8;
            const double h = 1.0 / std::pow(4.0, n);

            VertexR2 tmp(p.x() + h, p.y());
            const double f_xph = getW(tmp);

            tmp = VertexR2(p.x() - h, p.y());
            const double f_xmh = getW(tmp);

            tmp = VertexR2(p.x() + 2.0 * h, p.y());
            const double f_xp2h = getW(tmp);

            tmp = VertexR2(p.x() - 2.0 * h, p.y());
            const double f_xm2h = getW(tmp);

            return (f_xph - f_xmh) / (2.0 * h) - (f_xp2h - 2.0 * f_xph + 2.0 * f_xmh - f_xm2h) / (12.0 * h);
        }

        // Get y partial derivative dWy of the Wachspress denominator W.
        double getdWy(const VertexR2 &p) const {

            const int n = 8;
            const double h = 1.0 / std::pow(4.0, n);

            VertexR2 tmp(p.x(), p.y() + h);
            const double f_yph = getW(tmp);

            tmp = VertexR2(p.x(), p.y() - h);
            const double f_ymh = getW(tmp);

            tmp = VertexR2(p.x(), p.y() + 2.0 * h);
            const double f_yp2h = getW(tmp);

            tmp = VertexR2(p.x(), p.y() - 2.0 * h);
            const double f_ym2h = getW(tmp);

            return (f_yph - f_ymh) / (2.0 * h) - (f_yp2h - 2.0 * f_yph + 2.0 * f_ymh - f_ym2h) / (12.0 * h);
        }

        // Compute Wachspress coordinates.
        void computeWP(const VertexR2 &p, std::vector<double> &b) const {

            b.clear();

            const size_t n = _v.size();
            b.resize(n, 0.0);

            // Boundary.
            if (computeBoundaryCoordinates(p, b)) return;

            // Interior.
            std::vector<double> w;
            const double W = computeWPWeights(p, w);

            assert(w.size() == n);
            assert(fabs(W) > 0.0);

            const double invW = 1.0 / W;

            for (size_t i = 0; i < n; ++i) b[i] = w[i] * invW;
        }

        // Compute Wachspress weights.
        double computeWPWeights(const VertexR2 &p, std::vector<double> &w) const {

            w.clear();

            const size_t n = _v.size();
            w.resize(n, 0.0);

            std::vector<VertexR2> s(n);
            std::vector<VertexR2> e(n);

            for (size_t i = 0; i < n; ++i) {
                const size_t ip = (i + 1) % n;

                s[i] = _v[i] - p;
                e[i] = _v[ip] - _v[i];
            }

            std::vector<double> A(n, 1.0);
            std::vector<double> C(n);

            for (size_t i = 0; i < n; ++i) {
                const size_t im = (i + n - 1) % n;

                for (size_t j = 0; j < n; ++j) {
                    if (j != im && j != i) {

                        const size_t jp = (j + 1) % n;
                        A[i] *= 0.5 * s[j].crossProduct(s[jp]);
                    }
                }

                C[i] = 0.5 * (e[i].crossProduct(-e[im]));
            }

            double W = 0.0;
            for (size_t i = 0; i < n; ++i) {

                w[i] = C[i] * A[i];
                W += w[i];
            }

            return W;
        }

        // Compute cotangent.
        inline double cotangent(const VertexR2 &a, const VertexR2 &b) const {
            return a.scalarProduct(b) / fabs(a.crossProduct(b));
        }
    };

} // namespace gbc

#endif // GBC_BLENDEDWACHSPRESSR2_HPP
