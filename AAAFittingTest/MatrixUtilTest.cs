using MultiPrecision;
using MultiPrecisionComplex;
using MultiPrecisionComplexAlgebra;

using AAAFitting;

namespace AAAFittingTest {
    [TestClass]
    public class MatrixUtilTest {
        [TestMethod]
        public void SmallestSingularValueVectorNxMTest() {
            ComplexMatrix<Pow2.N4> m = new(new Complex<Pow2.N4>[,]
                {{ ( 1, + 2), ( 3, - 1), ( 5, + 0) },
                 { ( 4, + 0), ( 6, + 2), ( 8, - 3) },
                 { ( 7, - 1), ( 9, + 1), (11, - 2) },
                 { (10, + 0), (12, - 1), (14, + 2) } }
            );

            (ComplexVector<Pow2.N4> vec, Complex<Pow2.N4> val) = MatrixUtil<Pow2.N4>.SmallestSingularValueVector(m);

            ComplexMatrix<Pow2.N4> r = m.T.Conj * m;

            ComplexVector<Pow2.N4> a = r * vec;
            ComplexVector<Pow2.N4> b = val * vec;

            Console.WriteLine(vec);
            Console.WriteLine(val);
            Console.WriteLine(a);
            Console.WriteLine(b);

            Assert.IsTrue((a - b).Norm < 1e-32);
        }

        [TestMethod]
        public void SmallestSingularValueVectorNxNTest() {
            ComplexMatrix<Pow2.N4> m = new(new Complex<Pow2.N4>[,]
                {{ ( 1, + 2), ( 3, - 1), ( 5, + 0), ( 4, + 2) },
                 { ( 4, + 0), ( 6, + 2), ( 8, - 3), ( 3, + 1) },
                 { ( 7, - 1), ( 9, + 1), (11, - 2), (12, - 3) },
                 { (10, + 0), (12, - 1), (14, + 2), ( 7, - 4) } }
            );

            (ComplexVector<Pow2.N4> vec, Complex<Pow2.N4> val) = MatrixUtil<Pow2.N4>.SmallestSingularValueVector(m);

            ComplexMatrix<Pow2.N4> r = m.T.Conj * m;

            ComplexVector<Pow2.N4> a = r * vec;
            ComplexVector<Pow2.N4> b = val * vec;

            Console.WriteLine(vec);
            Console.WriteLine(val);
            Console.WriteLine(a);
            Console.WriteLine(b);

            Assert.IsTrue((a - b).Norm < 1e-32);
        }
    }
}