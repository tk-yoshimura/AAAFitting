using MultiPrecision;
using MultiPrecisionComplex;
using MultiPrecisionComplexAlgebra;
using System.Collections.ObjectModel;

namespace AAAFitting {
    public class BarycentricRational<N> where N : struct, IConstant {
        public ReadOnlyCollection<(Complex<N> node, Complex<N> value, Complex<N> weight)> Parameters { get; }

        public int Points => Parameters.Count;

        public BarycentricRational(IEnumerable<(Complex<N> node, Complex<N> value, Complex<N> weight)> parameters) {
            Parameters = new ReadOnlyCollection<(Complex<N> node, Complex<N> value, Complex<N> weight)>(parameters.ToArray());
        }

        public BarycentricRational(IEnumerable<Complex<N>> nodes, IEnumerable<Complex<N>> values, IEnumerable<Complex<N>> weights) {
            if (nodes.Count() != values.Count() || nodes.Count() != weights.Count()) {
                throw new ArgumentException("mismatch length", $"{nameof(nodes)}, {nameof(values)}, {nameof(weights)}");
            }

            Parameters = new ReadOnlyCollection<(Complex<N> node, Complex<N> value, Complex<N> weight)>(
                nodes.Zip(values, weights).ToArray()
            );
        }

        public Complex<N> FittingValue(Complex<N> z) {
            Complex<N> n = Complex<N>.Zero, d = Complex<N>.Zero;

            foreach ((Complex<N> node, Complex<N> value, Complex<N> weight) in Parameters) {
                Complex<N> v = weight / (z - node);

                if (!Complex<N>.IsFinite(v)) {
                    return value;
                }

                n += value * v;
                d += v;
            }

            Complex<N> r = n / d;

            return r;
        }

        public ComplexVector<N> FittingValue(ComplexVector<N> z) {
            ComplexVector<N> r = (FittingValue, z);

            return r;
        }
    }
}
