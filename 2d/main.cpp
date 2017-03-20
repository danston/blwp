// Copyright Dmitry Anisimov danston@ymail.com (c) 2017.

// An example on how to compute a convex combination of Wachspress and mean value coordinates.

// Local includes.
#include "./coords/MeanValueR2.hpp"
#include "./coords/BlendedWachspressR2.hpp"

// Main function.
int main() {

    using namespace gbc;

    // Polygon.
    std::vector<VertexR2> poly(8);

    poly[0] = VertexR2(0.0 , 0.0 );
    poly[1] = VertexR2(0.34, 0.32);
    poly[2] = VertexR2(0.68, 0.08);
    poly[3] = VertexR2(1.12, 0.31);
    poly[4] = VertexR2(1.12, 1.11);
    poly[5] = VertexR2(0.49, 1.07);
    poly[6] = VertexR2(0.72, 0.49);
    poly[7] = VertexR2(0.14, 0.58);

    // Evaluation point.
    VertexR2 query(0.56, 0.555);

    // Storage for the computed coordinates.
    std::vector<double> b;

    // Compute coordinates.
    BlendedWachspressR2<MeanValueR2> wmc(poly);
    wmc.compute(query, b);

    // Output the resulting coordinates.
    std::cout << "\nResult: ";
    for (size_t i = 0; i < b.size(); ++i) std::cout << b[i] << " ";
    std::cout << "\n\n";
}
