using MultiPrecision;
using MultiPrecisionComplex;
using MultiPrecisionComplexAlgebra;
using System.Diagnostics;

namespace AAAFitting {
    internal static class MatrixUtil<N> where N : struct, IConstant {
        public static (ComplexVector<N>, Complex<N>) SmallestSingularValueVector(ComplexMatrix<N> m) {
            Debug.Assert(m.Rows >= m.Columns);

            ComplexMatrix<Plus4<N>> m_ex = m.Convert<Plus4<N>>();

            ComplexMatrix<Plus4<N>> r = m_ex.H * m_ex;

            (Complex<Plus4<N>>[] eigen_vals, ComplexVector<Plus4<N>>[] eigen_vecs)  = ComplexMatrix<Plus4<N>>.EigenValueVectors(r);

            return (eigen_vecs[^1].Convert<N>(), eigen_vals[^1].Convert<N>());
        }
    }
}
