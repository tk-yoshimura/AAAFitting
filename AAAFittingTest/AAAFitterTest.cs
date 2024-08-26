using MultiPrecision;
using MultiPrecisionAlgebra;
using MultiPrecisionComplex;
using MultiPrecisionComplexAlgebra;

using AAAFitting;

namespace AAAFittingTest {
    [TestClass]
    public class AAAFitterTest {
        [TestMethod]
        public void FittingRealTest() {
            Vector<Pow2.N4> z = Enumerable.Range(0, 101).Select(i => MultiPrecision<Pow2.N4>.Div(i, 100)).ToArray();
            Vector<Pow2.N4> f = (v => MultiPrecision<Pow2.N4>.Exp(v) * MultiPrecision<Pow2.N4>.Sin(2 * MultiPrecision<Pow2.N4>.PI * v), z);

            BarycentricRational<Pow2.N4> approx = AAAFitter<Pow2.N4>.ExecuteFitting(z, f, 1e-16, 1e-16);

            ComplexVector<Pow2.N4> r = approx.FittingValue(z);

            for (int i = 0; i < z.Dim; i++) {
                Console.WriteLine($"{f[i]:e25}");
                Console.WriteLine($"{r[i]:e20}");
            }

            for (int i = 0; i < z.Dim; i++) {
                Assert.IsTrue(MultiPrecision<Pow2.N4>.Abs(f[i] - r[i].R) < 1e-16);
                Assert.IsTrue(r[i].I == 0);
            }
        }
    }
}