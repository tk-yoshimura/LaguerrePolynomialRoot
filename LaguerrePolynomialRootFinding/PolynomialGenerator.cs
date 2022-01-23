using MultiPrecision;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace LaguerrePolynomialRootFinding {
    internal static class PolynomialGenerator {
        static readonly List<BigInteger[]> coef_table = new();
        static readonly List<BigInteger> frac_table = new();

        static PolynomialGenerator() {
            coef_table.Add(new BigInteger[] { 1 });
            coef_table.Add(new BigInteger[] { 1, -1 });
            frac_table.Add(1);
            frac_table.Add(1);
        }

        public static Polynomial<N> Table<N>(int n) where N : struct, IConstant {
            if (n < coef_table.Count) {
                return new Polynomial<N>(coef_table[n].Select(c => (MultiPrecision<N>)c / frac_table[n]).ToArray());
            }

            for (int i = coef_table.Count; i <= n; i++) {
                BigInteger[] coef = new BigInteger[i + 1];

                coef[0] = coef_table[i - 1][0] * (2 * i - 1) - coef_table[i - 2][0] * (i - 1) * (i - 1);

                for (int k = 1; k <= i - 2; k++) {
                    coef[k] = coef_table[i - 1][k] * (2 * i - 1) - coef_table[i - 1][k - 1] - coef_table[i - 2][k] * (i - 1) * (i - 1);
                }

                coef[i - 1] = coef_table[i - 1][i - 1] * (2 * i - 1) - coef_table[i - 1][i - 2];
                coef[i] = -coef_table[i - 1][i - 1];

                coef_table.Add(coef);
                frac_table.Add(frac_table.Last() * i);
            }

            return new Polynomial<N>(coef_table[n].Select(c => (MultiPrecision<N>)c / frac_table[n]).ToArray());
        }
    }
}
