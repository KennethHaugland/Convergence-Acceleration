using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace ConvergenceAcceleration
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();

            LeibnitzFomula.Formula = @"\pi = 4 \cdot \sum_{n=0}^\infty (-1)^n \frac{1}{2n+1} = 4 - \frac{4}{3} + \frac{4}{5} - \cdots = " + Math.PI.ToString("N10");
            BaselFormula.Formula= @"\frac{\pi^2}{6} = \sum_{n=1}^\infty \frac{1}{n^2} = \frac{1}{1^2}+\frac{1}{2^2}+\frac{1}{3^2}+\frac{1}{4^2}+\cdots = " + (Math.PI*Math.PI/6d).ToString("N10");
            Log2Formula.Formula = @"\log(2) = \sum_{n=1}^\infty \frac{(-1)^{k-1}}{k} = " + Math.Log(2).ToString("N10");//           
        }

        private void btnCalculatePi(object sender, RoutedEventArgs e)
        {

            try
            {
                int terms = int.Parse(txtPiTerms.Text);

                double[] pi_leibnitz = MathExtensions.Convergence.PartialSumPi_Leibnitz(terms);
                txtPiStandardResult.Text = pi_leibnitz.Last().ToString("N10");
                txtPiAitkin.Text = MathExtensions.Convergence.Aitken(pi_leibnitz).ToString("N10");
                txtPiEpsilon.Text = MathExtensions.Convergence.EpsilonAlgorithm(pi_leibnitz).ToString("N10");
                txtPiLevin.Text = MathExtensions.Convergence.LevinsAlgorithm(pi_leibnitz, 4, 5).ToString("N10");
            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.Message);
                
            }
            


        }

        private void btnCalculateLog2(object sender, RoutedEventArgs e)
        {
            try
            {
                int terms = int.Parse(txtLog2Terms.Text);
                double[] LOG_2 = MathExtensions.Convergence.Log2(terms);

                double LOG2_SUM = LOG_2.Sum();
                double LOG2 = System.Math.Log(2);
                double LOG2_EULER = MathExtensions.Convergence.EulerTransformation(LOG_2.ToArray()).Sum();

                txtStandardLog2Summation.Text = LOG_2.Sum().ToString("N10");
                txtEulerLog2Summation.Text = LOG2_EULER.ToString("N10");

                double Log2_error = LOG2_SUM - LOG2;
                double Log2_euler_error = LOG2_EULER - LOG2;
            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.Message);
            }

        }

        private void btnCalculateBasel(object sender, RoutedEventArgs e)
        {

            try
            {
                int terms = int.Parse(txtBaselTerms.Text);
                double[] BaselNumbers = MathExtensions.Convergence.BaselProblem(terms);
                double RichardsonBaselNumber = MathExtensions.Convergence.RichardsonExtrapolation(BaselNumbers);
                txtAcceleratedBaselSereis.Text = RichardsonBaselNumber.ToString("N10");
                txtBaselSum.Text = BaselNumbers.Last().ToString("N10");
            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.Message);
            }

        }
    }
}
