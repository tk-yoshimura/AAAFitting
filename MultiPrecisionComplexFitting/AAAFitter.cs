using MultiPrecision;
using MultiPrecisionAlgebra;
using MultiPrecisionComplex;
using MultiPrecisionComplexAlgebra;

namespace MultiPrecisionComplexFitting {
    // [Reference]
    // The AAA Algorithm for Rational Approximation
    // Yuji Nakatsukasa, Olivier Sete, and Lloyd N. Trefethen
    // SIAM Journal on Scientific Computing 2018 40:3, A1494-A1522
    // https://doi.org/10.1137/16M1106122
    public static class AAAFitter<N> where N : struct, IConstant {
        public static BarycentricRational<N> Approx(
            ComplexVector<N> z, ComplexVector<N> f,
            MultiPrecision<N> reltol, MultiPrecision<N> abstol, int max_points = 128) {

            if (z.Dim != f.Dim) {
                throw new ArgumentException("mismatch count", $"{nameof(z)},{nameof(f)}");
            }
            if (!(reltol >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(reltol));
            }
            if (!(abstol >= 0)) {
                throw new ArgumentOutOfRangeException(nameof(abstol));
            }
            ArgumentOutOfRangeException.ThrowIfLessThan(max_points, 2, nameof(max_points));

            Vector<N> eps = f.Select(v => v.val.Magnitude * reltol + abstol).ToArray();
            List<Complex<N>> nodes = [], values = [];

            ComplexVector<N> r = ComplexVector<N>.Fill(f.Dim, f.Mean);

            List<int> indexes = Enumerable.Range(0, f.Dim).ToList();

            while (nodes.Count < max_points && nodes.Count <= indexes.Count) {
                List<MultiPrecision<N>> error_norms = (f - r).Select(v => v.val.Norm).ToList();
                int index_maxerr = error_norms.IndexOf(error_norms.Max());

                nodes.Add(z[index_maxerr]);
                values.Add(f[index_maxerr]);

                indexes.Remove(index_maxerr);

                ComplexMatrix<N> c = ComplexMatrix<N>.Zero(indexes.Count, nodes.Count);
                ComplexMatrix<N> a = ComplexMatrix<N>.Zero(indexes.Count, nodes.Count);

                for (int i = 0; i < indexes.Count; i++) {
                    for (int j = 0; j < nodes.Count; j++) {
                        c[i, j] = Complex<N>.Inverse(z[indexes[i]] - nodes[j]);
                        a[i, j] = (f[indexes[i]] - values[j]) * c[i, j];
                    }
                }

                (ComplexVector<N> weights, _) = MatrixUtil<N>.SmallestSingularValueVector(a);

                BarycentricRational<N> approx = new(nodes, values, (Complex<N>[])weights);
            }

            throw new NotImplementedException();
        }
    }
}
