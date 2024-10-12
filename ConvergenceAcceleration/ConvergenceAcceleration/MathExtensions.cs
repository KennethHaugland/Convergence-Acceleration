using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace MathExtensions
{
    public static class Convergence
    {

        /// <summary>
        /// Calculates the Aitken summation
        /// </summary>
        /// <param name="S_n">The partial sums of a series</param>
        /// <returns>A extrapolated value for the partial sum</returns>
        /// <exception cref="Exception">If you have less then 3 items this will return with an exception</exception>
        public static double Aitken(double[] S_n)
        {
            if (S_n.Length > 3)
            {
                double S_2 = S_n[S_n.Length - 1];
                double S_1 = S_n[S_n.Length - 2];
                double S_0 = S_n[S_n.Length - 3];

                return (S_2 * S_0 - S_1 * S_1) / (S_2 - 2 * S_1 + S_0);
            }
            else
                throw new Exception("Must have at least 3 partial sums to give the correct value");
        }

        /// <summary>
        /// Returns all the partial sums with n iterations of Pi
        /// </summary>
        /// <param name="n">Number of iterations</param>
        /// <returns></returns>
        public static double[] PartialSumPi_Leibnitz(int n = 6)
        {
            double sum = 0;
            double sign = 1;
            double[] result = new double[n];
            int ind = 0;
            for (int i = 1; i < n * 2; i = i + 2)
            {
                sum = sum + sign * 4 / i;
                sign *= -1;
                result[ind] = sum;
                ind += 1;
            }
            return result;
        }

        /// <summary>
        /// Partial serie sums of pi^2/6 = 1.644934066848
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double[] BaselProblem(int n = 6)
        {
            List<double> result = new List<double>();
            double sum = 0;
            for (int i = 1; i < n + 1; i++)
            {
                sum += 1 / System.Math.Pow((double)i, 2);
                result.Add(sum);
            }

            return result.ToArray();
        }

        /// <summary>
        /// Levins algorithm L_{k}^{n}(S_n) but k + n has to be less then 
        /// the length of S_n minus one
        /// </summary>
        /// <param name="S_n">The partial truncated sums</param>
        /// <param name="k">The 'order' of Levins algorithm</param>
        /// <param name="n">The n'th point approximation</param>
        /// <returns>The accelerated value</returns>
        public static double LevinsAlgorithm(double[] S_n, double k = 4, double n = 5, double b = 1)
        {

            double numerator = 0;
            double denominator = 0;

            if (k + n > S_n.Length - 1)
                throw new ArgumentException("Not ennough elements to sum over with the n'th element with order k");

            for (int j = 0; j < n; j++)
            {

                double rest = System.Math.Pow(-1d, j) * GetBinCoeff((double)k, (double)j);

                double C_jkn_U = System.Math.Pow((double)(n + j + b), k - 1);
                double C_jkn_L = System.Math.Pow((double)(n + k + b), k - 1);

                double C_njk = C_jkn_U / C_jkn_L;

                double S_nj = S_n[(int)n + (int)j];

                // t transform that calculates a_n
                double g_n = S_n[(int)n + j] - S_n[(int)n + j - 1];
                // u transform that calcualtes (n+k) * a_n
                // g_n = (n + k)* (S_n[(int)n + j] - S_n[(int)n + j - 1]);

                numerator += rest * C_njk * S_nj / g_n;
                denominator += rest * C_njk / g_n;
            }

            return numerator / denominator;
        }

        /// <summary>
        /// Peter Wynns epsilon algorithm for calculating accelerated convergence
        /// </summary>
        /// <param name="S_n">The partial sums</param>
        /// <returns>The best accelerated sum it finds</returns>
        public static double EpsilonAlgorithm(double[] S_n, bool Logaritmic = false)
        {

            int m = S_n.Length;

            double[,] r = new double[m + 1, m + 1];

            // Fill inn the partial sums in the 1 column
            for (int i = 0; i < m; i++)
                r[i, 1] = S_n[i];


            // Epsilon algorithm 
            for (int column = 2; column <= m; column++)
            {
                for (int row = 0; row < m - column + 1; row++)
                {
                    //Check for divisions of zero (other checks should be done here)
                    double divisor = (r[row + 1, column - 1] - r[row, column - 1]);

                    // Epsilon
                    double numerator = 1;

                    if (Logaritmic)
                        numerator = column + 1;

                    if (divisor != 0)
                        r[row, column] = r[row + 1, column - 2] + numerator / divisor;
                    else
                        r[row, column] = 0;
                }
            }



            // Clean up, only interested in the odd number columns
            int items = (int)System.Math.Floor((double)((m + 1) / 2));
            double[,] e = new double[m, items];

            for (int row = 0; row < m; row++)
            {
                int index = 0;
                for (int column = 1; column < m + 1; column = column + 2)
                {
                    if (row + index == m)
                        break;

                    //e[row + index, index] = r[row, column];
                    e[row, index] = r[row, column];
                    index += 1;
                }
            }
            return e[0, e.GetLength(1) - 1];
        }


        /// <summary>
        /// Generalized Richardson extrapolation
        /// </summary>
        /// <param name="S_n">The partial sums</param>
        /// <returns>Best approximation of S</returns>
        public static double RichardsonExtrapolation(double[] S_n)
        {
            double S = 0d;

            double m = S_n.Length;
            for (int j = 1; j <= (int)m; j++)
            {
                S += System.Math.Pow(-1, m + (double)j) *
                    S_n[j - 1] *
                    System.Math.Pow((double)j, m - 1)
                    / (Factorial(j - 1) * Factorial(m - j));
            }

            return S;
        }

        /// <summary>
        /// Slow convergence Log(2) series
        /// </summary>
        /// <param name="n">Number of iterations</param>
        /// <returns></returns>
        public static double[] Log2(int n = 6)
        {
            // Code for log2
            List<double> LOG_2 = new List<double>();
            for (int i = 1; i <= n; i++)
            {
                LOG_2.Add(System.Math.Pow(-1d, (double)(i - 1)) * 1 / (double)i);
            }
            return LOG_2.ToArray();
        }

        /// <summary>
        /// Euler transform that transforms the alternating series a_0 into a faster convergence with no negative coefficients
        /// </summary>
        /// <param name="a_0">The alternating power series</param>
        /// <returns></returns>
        public static double[] EulerTransformation(double[] a_0)
        {
            // Each series item
            List<double> a_k = new List<double>();

            // finite difference of each item
            double delta_a_0 = 0;


            for (int k = 0; k < a_0.Length; k++)
            {
                delta_a_0 = 0;
                for (int m = 0; m <= k; m++)
                {
                    double choose_k_over_m = (double)GetBinCoeff(k, m);
                    delta_a_0 += System.Math.Pow(-1d, (double)m) * choose_k_over_m * System.Math.Abs(a_0[(int)(m)]);
                }

                a_k.Add(System.Math.Pow(1d / 2d, (double)(k + 1)) * delta_a_0);
            }
            return a_k.ToArray();
        }

        /// <summary>
        /// Binominal coefficient
        /// </summary>
        /// <param name="N"></param>
        /// <param name="K"></param>
        /// <returns>N choose K</returns>
        public static double GetBinCoeff(double N, double K)
        {
            double r = 1;

            for (double d = 1; d <= K; d++)
            {
                r *= ((N--) / d);
            }
            return r;
        }

        /// <summary>
        /// Returns factorial, integers only
        /// </summary>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double Factorial(double n)
        {
            double result = 1;
            for (int i = 1; i <= n; i++)
            {
                result *= i;
            }
            return result;
        }

        /// <summary>
        /// Returns a double tabulated table with all values of the matix. 
        /// Zero entries are not shown
        /// </summary>
        /// <param name="A">The input matrix</param>
        /// <returns></returns>
        static public string ToString(this double[,] A)
        {
            int row = A.GetLength(0);
            int col = A.GetLength(1);
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < col; j++)
            {
                sb.Append("S^" + j.ToString());
                sb.Append("\t\t\t");
            }
            sb.AppendLine();

            for (int i = 0; i < row; i++)
            {

                for (int j = 0; j < col; j++)
                {
                    if (A[i, j] != 0)
                        sb.Append(A[i, j].ToString("N10"));

                    sb.Append("\t\t");
                }
                sb.AppendLine();
            }

            return sb.ToString();
        }
    }
}

