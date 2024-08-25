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
        public static BarycentricRational<N> ExecuteFitting(
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

            BarycentricRational<N> approx = new(Enumerable.Empty<(Complex<N>, Complex<N>, Complex<N>)>());

            while (nodes.Count < max_points && nodes.Count <= indexes.Count) {
                List<MultiPrecision<N>> errors = (f - r).Select(v => v.val.Magnitude).ToList();
                int index_maxerror = -1;
                MultiPrecision<N> maxerror = 0;

                bool is_allok = true;
                
                for (int i = 0; i < errors.Count; i++) {
                    if (maxerror < errors[i]) {
                        maxerror = errors[i];
                        index_maxerror = i;
                    }
                    if (errors[i] > eps[i]) {
                        is_allok = false;
                    }
                }

                if (is_allok && nodes.Count >= 1) {
                    break;
                }

                nodes.Add(z[index_maxerror]);
                values.Add(f[index_maxerror]);

                indexes.Remove(index_maxerror);

                ComplexMatrix<N> a = ComplexMatrix<N>.Zero(indexes.Count, nodes.Count);

                for (int i = 0; i < indexes.Count; i++) {
                    for (int j = 0; j < nodes.Count; j++) {
                        a[i, j] = (f[indexes[i]] - values[j]) / (z[indexes[i]] - nodes[j]);
                    }
                }

                (ComplexVector<N> weights, _) = MatrixUtil<N>.SmallestSingularValueVector(a);

                approx = new(nodes, values, (Complex<N>[])weights);

                r = approx.FittingValue(z);
            }

            return approx;
        }
    }
}
