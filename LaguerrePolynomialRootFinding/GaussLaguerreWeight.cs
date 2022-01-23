using MultiPrecision;

namespace LaguerrePolynomialRootFinding {
    internal static class GaussLaguerreWeight<N> where N : struct, IConstant {
        public static MultiPrecision<N> Eval(MultiPrecision<N> x_root, int n, Polynomial<N> poly_np1) {
            MultiPrecision<N> w = x_root / MultiPrecision<N>.Square((n + 1) * poly_np1.Value(x_root));

            return w;
        }
    }
}
