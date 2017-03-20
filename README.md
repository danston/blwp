# A convex combination of Wachspress and any other 2D closed-form barycentric coordinates Ver. 1.0.0

### Description

This set of C++ classes provides you with an implementation of a convex combination of Wachspress and any other 2D closed-form barycentric coordinates as discussed in Chapter 4 of the dissertation: D. Anisimov. Analysis and new constructions of generalized barycentric coordinates in 2D. To appear, 2017. The implementation is based on the corresponding pseudocode from Appendix A of this dissertation.

##### NOTE: This code has been tested only on Mac OS!

### Run the code

In order to run the code, do the following:

1. Install [macports](https://www.macports.org/install.php)
2. Open terminal and type the following:

```bash
  sudo port install cmake
```
```bash
  cd path_to_the_folder/blwp/2d/
```
```bash
  mkdir bin
```
```bash
  cd bin
```
```bash
  cmake -DCMAKE_BUILD_TYPE=Debug ..
```
```bash
  make
```
```bash
  ./blwp
```

For the release version use instead: 

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
```

### Example

```C++
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
```

##### NOTE: For the complete example see main.cpp!

### Bugs

If you find any bugs, please report them to me, and I will try to fix them as soon as possible!
