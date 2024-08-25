using MultiPrecision;
using MultiPrecisionComplex;
using MultiPrecisionComplexAlgebra;
using System.Diagnostics;

namespace MultiPrecisionComplexFitting {
    internal static class MatrixUtil<N> where N : struct, IConstant {
        public static (ComplexVector<N>, Complex<N>) SmallestSingularValueVector(ComplexMatrix<N> m) {
            Debug.Assert(m.Rows >= m.Columns);

            ComplexMatrix<Plus4<N>> m_ex = m.Convert<Plus4<N>>();

            ComplexMatrix<Plus4<N>> r = (m_ex.T.Conj * m_ex).Inverse;
            int n = r.Size;

            MultiPrecision<Plus4<N>> prev_diff = MultiPrecision<Plus4<N>>.PositiveInfinity;

            (ComplexVector<Plus4<N>> vec, ComplexVector<Plus4<N>> prev_vec, Complex<Plus4<N>> lambda)[] vs
                = new (ComplexVector<Plus4<N>>, ComplexVector<Plus4<N>>, Complex<Plus4<N>>)[n * 2];

            for (int i = 0; i < n; i++) {
                ComplexVector<Plus4<N>> v0 = ComplexVector<Plus4<N>>.Zero(r.Size), v1 = ComplexVector<Plus4<N>>.Zero(r.Size);

                v0[i] = 1;
                v1[i] = (0, 1);

                vs[i].vec = vs[i].prev_vec = v0;
                vs[i + n].vec = vs[i + n].prev_vec = v1;

                vs[i].lambda = vs[i + n].lambda = 0;
            }

            for ((int iter, bool convergenced) = (1, false); iter <= 65536 && !convergenced; iter++) {
                for (int i = 0; i < vs.Length; i++) {
                    vs[i].vec = r * vs[i].vec;
                }

                if ((iter % 4) == 0) {
                    for (int i = 0; i < vs.Length; i++) {
                        vs[i].lambda =
                            ComplexVector<Plus4<N>>.Dot(vs[i].vec, vs[i].prev_vec) /
                            ComplexVector<Plus4<N>>.Dot(vs[i].vec, vs[i].vec);

                        vs[i].vec = vs[i].vec.Normal;
                    }

                    vs = [.. vs.OrderBy(item => item.lambda.Norm)];

                    MultiPrecision<Plus4<N>> diff = (vs[0].prev_vec - vs[0].vec).Norm;

                    if (diff < 1e-4) {
                        vs = vs[..1];
                    }

                    if (diff < prev_diff) {
                        prev_diff = diff;
                    }
                    else {
                        convergenced = true;
                        break;
                    }
                }

                for (int i = 0; i < vs.Length; i++) {
                    vs[i].prev_vec = vs[i].vec;
                }
            }

            ComplexVector<Plus4<N>> eigen_vec = vs[0].vec;
            Complex<Plus4<N>> eigen_val = vs[0].lambda;

            if (eigen_vec[0].R.Sign == Sign.Minus) {
                eigen_vec = -eigen_vec;
            }

            long truncate_bits = -MultiPrecision<N>.Bits;

            eigen_vec = eigen_vec.Select(
                v => new Complex<Plus4<N>>(
                    v.val.R.Exponent >= truncate_bits ? v.val.R : 0,
                    v.val.I.Exponent >= truncate_bits ? v.val.I : 0
                )
            ).ToArray();

            eigen_val = (
                eigen_val.R.Exponent >= truncate_bits ? eigen_val.R : 0,
                eigen_val.I.Exponent >= truncate_bits ? eigen_val.I : 0
            );

            return (eigen_vec.Convert<N>(), eigen_val.Convert<N>());
        }
    }
}
