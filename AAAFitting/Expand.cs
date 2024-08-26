using MultiPrecision;

namespace AAAFitting {
    internal struct Plus4<N> : IConstant where N : struct, IConstant {
        public readonly int Value => checked(default(N).Value + 4);
    }
}
