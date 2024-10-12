<ul class="download">
	<li><a href="ConvergenceAcceleration.zip">Download source code - 444.6 KB</a></li>
</ul>

<p><img src="ConvergenceScreenShot.png" style="height: 646px; width: 640px" /></p>

<h2>Introduction</h2>

<p>Many applications implements some calculation of a series or sequence&nbsp;that can have a very slow convergence rate or might require more precision than you could get from a brute force implementation. There are several candidates to choose from when implementing these acceleration techniques and I will present the most common ones in this article.</p>

<p>However, there is the the issue that you will not find&nbsp;one&nbsp;universal acceleration technique that will work for all sequences, which is proven mathematically. This proof is based on a property of the sequence called remanence if you are interested in looking it up for yourself. This gives you the hint of either using a mathematical test/prove before releasing it to your costumers or at least your own judgement to test the acceleration technique on the series first. And you should not expect the same magic to work universally.</p>

<h2>Background</h2>

<p>Assume that we have sequence or series which have the actual (theoretical) sum is <span class="math">\(S\)</span> but we only have the ability to calculate the approximation <span class="math">\(S_n\)</span> :</p>

<div class="math">$S_n = \sum_{n=0}^k{a_n}$</div>

<p>This analysis can also be used on an infinite product:</p>

<div class="math">$S_n = \Pi_{n=0}^k{a_n}$</div>

<p>So, describing the problem at hand in the most general way, we simply want to minimize the difference <span class="math">\(S - S_n = e\)</span> or if possible to eliminate the error <span class="math">\(e\)&nbsp;</span>. Estimating the form that the error function has and use this information to create a correction is called the kernel, and it&nbsp;is actually the only way we know of making these faster converging or more precise algorithm. If we know how the error grows after each iterative process, we can in fact use this knowledge to predict a better real value then a brute force calculation. And also, since we must guess values that lie beyond what we currently know, we call all these implementations with kernel correcting errors extrapolation methods in mathematics.</p>

<h2>Aitken&#39;s Summation Formula</h2>

<p>Many series that we perform summations on grows&nbsp;geometrical by nature, so at first it seems natural to assume that the error and the correcting kernel can be described by a formula&nbsp;&nbsp;<span class="math">\(b \cdot q^n\). The kernel we want to use&nbsp;</span>in order to try to correct and adjust the complete sum is:</p>

<div class="math">$S_n \approx S + b \cdot q^n$</div>

<p>We assume here that <span class="math">\(q &lt; 1\)</span> so that the error actually decrease for each iterative term we sum over in the original series. With the general kernel above we have 3 unknowns <span class="math">\(S, b, q\)</span> in the equation to solve for, which means that we need 3 equations in order to get a unique solution. The result is what is known as Aitken&#39;s formula:</p>

<div class="math">$S \approx \frac{(S_{n+1} S_{n-1} - S_n^2)}{S_{n+1} - 2 S_n + S_{n-1}}$</div>

<p>The first discoverer of this formula seems to be the Japanese mathematician Takakazu Seki in around 1680. (He is sometimes referred to by mistake as Seki Kowa due to some spelling similarities in Japanese). It has also been discovered independently by others before the name&nbsp;of this way of&nbsp;accelerate the&nbsp;sequence was given to Aitken.&nbsp;</p>

<p>The code&nbsp;is pretty straight forward to implement:</p>

<pre lang="cs">
        /// &lt;summary&gt;
        /// Calculates the Aitken summation
        /// &lt;/summary&gt;
        /// &lt;param name=&quot;S_n&quot;&gt;The partial sums of a series&lt;/param&gt;
        /// &lt;returns&gt;A extrapolated value for the partial sum&lt;/returns&gt;
        /// &lt;exception cref=&quot;Exception&quot;&gt;If you have less than 3 items, 
        /// this will return with an exception&lt;/exception&gt;
        public double Aitken(double[] S_n)
        {
            if (S_n.Length &gt; 3)
            {
                double S_2 = S_n[S_n.Length - 1];
                double S_1 = S_n[S_n.Length - 2];
                double S_0 = S_n[S_n.Length - 3];

                return (S_2 * S_0 - S_1 * S_1) / (S_2 -2*S_1 + S_0);
            }
            else
                throw new Exception(&quot;Must have at least 3 partial sums 
                                     to give the correct value&quot;);
        }</pre>

<p>Aitkens formula can indeed accelerate all linearly converging series, but in addition, it can also be used to accelerate the much more difficult logarithmic sequences.&nbsp;This algorithm can also be used in multiple stages, meaning that we apply&nbsp;it over and over again on the same summation series.</p>

<h2>Shanks Transform and Wynns Epsilon Algorithm</h2>

<p>The general acceleration kernel of a geometric series was later formulated into the Shanks transform, were Shanks formulated the correction term with more variables:</p>

<div class="math">$S_n \approx S + \sum_{i=1}^m{b_i \cdot q_i^n}$</div>

<p>This will lead&nbsp;to a 2m +1 unknowns with the general determinant ration found by Shanks in 1947:</p>

<div class="math">$S \approx \frac{\begin{vmatrix}S_{n-m} &amp; S_{n-m +1} &amp; \cdots &amp; S_{n} \\ S_{n-m+1} - S_{n-m} &amp; S_{n-m+2} - S_{n-m+1} &amp; \cdots &amp; S_{n+1} - S_{n} \\ S_{n-m+2} - S_{n-m+1} &amp; S_{n-m+3} - S_{n-m+2} &amp; \cdots &amp; S_{n+2} - S_{n+1} \\ \vdots &amp; \vdots &amp; \ddots &amp; \vdots \\ S_{n} - S_{n-1} &amp; S_{n+1} - S_{n} &amp; \cdots &amp; S_{n+m} - S_{n+m+1} \\\end{vmatrix} }{\begin{vmatrix}1 &amp; 1 &amp; \cdots &amp; 1 \\ S_{n-m+1} - S_{n-m} &amp; S_{n-m+2} - S_{n-m+1} &amp; \cdots &amp; S_{n+1} - S_{n} \\ S_{n-m+2} - S_{n-m+1} &amp; S_{n-m+3} - S_{n-m+2} &amp; \cdots &amp; S_{n+2} - S_{n+1} \\ \vdots &amp; \vdots &amp; \ddots &amp; \vdots \\ S_{n} - S_{n-1} &amp; S_{n+1} - S_{n} &amp; \cdots &amp; S_{n+m} - S_{n+m+1} \\\end{vmatrix} }$</div>

<p>Most of the elements in the&nbsp;matrix consists of many finite difference items of consecutive terms, so these elements are usually abbreviated by <span class="math">\( \Delta S_p = S_{p+1} - S_p \). The absolute mark in front and behind the matrix means that you should calculate the determinant of it, and al</span>though the determinants can be found by creating a lower (or upper) triangular matrix and simply multiplying&nbsp;the resulting diagonal elements, it does not in general a stable method.&nbsp;The issue here is that the Henkel determinants of matrices on this form uses multiple arithmetical operations&nbsp; that can lead&nbsp;to significant truncation errors when implemented directly on computers. So using a different method is desirable.</p>

<p>There is also a more hidden pattern behind the full Shanks transform of a sequence since it is equivalent to the Pade approximation. One can actually use the Shanks transform to calculate the full Pade approximation, which is also the fastest way of finding a convergence series due to the links to&nbsp;continued fractions. But the same calculation&nbsp;problem remains, we need a better way of calculating this transformation.&nbsp;</p>

<p>Peter Wynn&#39;s <span class="math">\(\epsilon\)</span>-algorithm produces the same results by using an iterative process and in fact it is very similar to calculating the continued fraction of a series. We start off by assuming we have a table designated with X,Y were we define 0,0 at the top left corner.&nbsp;In the first column, we fill it with zeros <span class="math">\(\epsilon(0,y) = 0\)</span> and in the second column, we place the partial sum of each iteration <span class="math">\(\epsilon(1,y) = S_n\)</span>.&nbsp;We now use the iteration schema below:</p>

<div class="math">$\epsilon(x,y) = \epsilon(x-2,y+1) + \frac{1}{\epsilon(x-1,y+1) - \epsilon(x-1,y)}$</div>

<p>The convergence speedup is used&nbsp;on the slow converging Leibnitz formula for <span class="math">\(\pi\) and&nbsp;</span>can be illustrated with the table in the picture below:</p>

<p><img src="Wynns_epsilonIterator2.png" style="width: 640px; height: 445px" /></p>

<p>A better, or at least a more common way to illustrate the relation, is to use a lozenge diagram (image taken from <a href="https://www.sciencedirect.com/science/article/pii/S0377042700003551#FIG1">https://www.sciencedirect.com/science/article/pii/S0377042700003551#FIG1</a>). The reason is that all values in the lozenge diagram that are used in each of the iterative step is right next to each other in the diagram.</p>

<p><img height="198" src="lozange_diagram_of_epsilon.gif" width="493" /></p>

<p>From the calculations, you can see that every even column just contains some intermediary results that is used to calculate the higher degree transformation. These intermediary results are skippable by Wynn cross rule that allows you to only store the odd numbered columns. The issue here is that you need 3 starting columns with the parameters&nbsp;&nbsp;<span class="math">\( \epsilon_{-2}^{(n)} = \infty \),</span><span class="math">\( \epsilon_{-1}^{(n)} = 0 \)</span> and finally <span class="math">\( \epsilon_{0}^{(n)} = S_n \) to implement the iterative process below:&nbsp;</span></p>

<div class="math">$\left(\epsilon_{k+2}^{(n)}-\epsilon_{k}^{(n+1)}\right)^{-1}+\left(\epsilon_{k-2}^{(n+2)}-\epsilon_{k}^{(n+1)}\right)^{-1}=\left(\epsilon_{k}^{(n+2)}-\epsilon_{k}^{(n+1)}\right)^{-1}-\left(\epsilon_{k}^{(n)}-\epsilon_{k}^{(n+1)}\right)^{-1}$</div>

<p>There is&nbsp;no real advantage to implement this algorithm on a computer,&nbsp;except that it eliminates columns that have the intermediate calculation results, but it has the drawback or disadvantage of needing one more starting condition. The more interesting thing that could be mentioned is that&nbsp;Wynn also found an implementation that can calculate the accelerated sums with just one vector and two intermediary variables. But on the modern computer this is not an&nbsp;issue and I will skip this derivation.&nbsp; You might have guessed it already, but Peter Wynn&#39;s <span class="math">\(\epsilon\)</span>-algorithm could also be used to calculate the Hankel determinant, the Pade approximation&nbsp;and surprisingly also be used to solve the differential wave type equation called Korteweg&ndash;De Vries (KdV) equation.</p>

<p>Now we switch to the actual implementation of the code, and first off, we need the test series that&nbsp;the acceleration algorithm shall accelerate. &nbsp;We have already previously illustrated the algorithms usefulness&nbsp;with the slowly converging series that Leibnitz found for pi given by the series:</p>

<div class="math">$\frac{\pi}{4} = \sum_{k=0}^{k=n}\frac{(-1)^k}{2k +1}$</div>

<p>The sum following each new term of the truncated series is saved in a list and returned.</p>

<pre lang="cs">
        public static double[] PartialSumPi_Leibnitz(int n = 6)
        {
            double sum = 0;
            double sign = 1;
            double[] result = new double[n];
            int ind = 0;
            for (int i = 1; i &lt; n*2; i= i+2)
            {
                sum = sum + sign * 4 / i;
                sign *= -1;
                result[ind] = sum;
                ind += 1;                
            }
            return result;
        }</pre>

<p>The actual code for Wynn&#39;s algorithm is given below:</p>

<pre lang="cs">
        /// &lt;summary&gt;
        /// Peter Wynns epsilon algorithm for calculating accelerated convergence
        /// &lt;/summary&gt;
        /// &lt;param name=&quot;S_n&quot;&gt;The partial sums&lt;/param&gt;
        /// &lt;returns&gt;The best accelerated sum it finds&lt;/returns&gt;
        public static double EpsilonAlgorithm(double[] S_n, bool Logaritmic = false)
        {

            int m = S_n.Length;

            double[,] r = new double[m + 1, m + 1];

            // Fill in the partial sums in the 1 column
            for (int i = 0; i &lt; m; i++)
                r[i, 1] = S_n[i];

            // Epsilon algorithm 
            for (int column = 2; column &lt;= m; column++)
            {
                for (int row = 0; row &lt; m - column + 1; row++)
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

            for (int row = 0; row &lt; m; row++)
            {
                int index = 0;
                for (int column = 1; column &lt; m + 1; column = column + 2)
                {
                    if (row + index == m)
                        break;

                    //e[row + index, index] = r[row, column];
                    e[row, index] = r[row, column];
                    index += 1;
                }
            }
            return e[0, e.GetLength(1) - 1];
        }</pre>

<p>Other algorithms based on the <span class="math">\(\epsilon\)</span>-algorithm have been proposed like the <span class="math">\(\rho\)</span>-algorithm also by Wynn:</p>

<div class="math">$\rho(x,y) = \rho(x-2,y+1) + \frac{x_{n+k+1} - x_n}{\rho(x-1,y+1) - \rho(x-1,y)}$</div>

<p>Here the <span class="math">\( x_n \)</span> is the interpolation points used in the sequence. In the usual case of integer sums the simplification <span class="math">\( x_n = n + 1\)&nbsp;</span>is used.</p>

<div class="math">$\rho(x,y) = \rho(x-2,y+1) + \frac{k + 1}{\rho(x-1,y+1) - \rho(x-1,y)}$</div>

<p>These algorithms generally works better for accelerating series that have logarithmically increasing errors.&nbsp;</p>

<h2>Richardson Extrapolation and Euler&#39;s Alternating Series Transformation</h2>

<p>Two similar methods for faster convergence that are often not grouped together, even though Richardson extrapolation is calculated with Taylor expansion and the Euler transform by finite difference approximation by a devilishly clever trick from its inventor Euler. We start off with the well known Taylor expansion of the unknown <span class="math">\( f(x) \)</span>&nbsp;were the function is adjusted by a small amount <span class="math">\( h \)</span>. And by the way.&nbsp;Taylor expansion is every mathematicians favorited method in applied mathematics. Not a day goes by without using it. Well,&nbsp;almost none anyway, so here it is incase you had forgotten it:</p>

<div class="math">$f(x \pm h) = f(x) \pm h \cdot f&#39;(x) + \frac{h^2}{2!} \cdot f&#39;&#39;(x) \pm \frac{h^3}{3!} \cdot f^{(3)}(x) + \frac{h^4}{4!} \cdot f^{(4)}(x) $</div>

<p>In the Richardson extrapolation we utilize the structure of the Taylor series with different multiple values of <span class="math">\( h \)</span> as <span class="math">\( 2h, 3h, \cdots \). This information is u</span>sed by a system of equations that allows us to eliminate&nbsp;higher and higher order derivatives in the Taylor expansion. In the generalized Richardson extrapolation, we take the idea of the Taylor expansion and use the integer sum the value <span class="math">\( \frac{1}{n} \)</span> is used instead of <span class="math">\( h \) and&nbsp;</span>reformulate the sum with the error estimate:</p>

<div class="math">$S_n = S + \frac{Q_1}{n^{1}} + \frac{Q_2}{n^{2}} + \cdots + \frac{Q_m}{n^{m}}$</div>

<p>The resulting closed form solution equations, that will give us a unique solution, are&nbsp;set up in the equations below. (Derivation and result is taken from &quot;Advanced mathematical methods for scientists and engineers&quot; by Bender ang Orszag)</p>

<div class="math">$\begin{split}S_n &amp;= S + \frac{Q_1}{n^{1}} + \frac{Q_2}{n^{2}} + \cdots + \frac{Q_m}{n^{m}} \\S_{n+1} &amp;= S + \frac{Q_1}{(n+1)^{1}} + \frac{Q_2}{(n+1)^{2}} + \cdots + \frac{Q_m}{(n+1)^{m}} \\ \vdots &amp;= \vdots \\ S_{n+m} &amp;= S + \frac{Q_1}{(n+m)^{1}} + \frac{Q_2}{(n+m)^{2}} + \cdots + \frac{Q_m}{(n+m)^{m}} \end{split}$</div>

<p>Solving for <span class="math">\(S\)</span> gives the solution to the m equations:</p>

<div class="math">$ S = \sum_{k=0}^{k = m}{\frac{S_{n+k} \cdot (n+k)^m \cdot (-1)^{k+m}}{k! (m-k)!}} $</div>

<p>The code for the Richardson extrapolation have the weakness in the dependence of the factorial, effectively limiting higher order implementations on the computer:</p>

<pre lang="cs">
        /// &lt;summary&gt;
        /// Generalized Richardson extrapolation
        /// &lt;/summary&gt;
        /// &lt;param name=&quot;S_n&quot;&gt;The partial sums&lt;/param&gt;
        /// &lt;returns&gt;Best approximation of S&lt;/returns&gt;
        public static double RichardsonExtrapolation(double[] S_n)
        {
            double S = 0d;

            double m = S_n.Length;
            for (int j = 1; j &lt;= (int)m; j++)
            {
                S += System.Math.Pow(-1, m + (double)j) * 
                    S_n[j - 1] * 
                    System.Math.Pow((double)j, m - 1) 
                    / (Factorial(j - 1) * Factorial(m - j));
            }

            return S;
        }</pre>

<p>The next item to showcase is the Euler series transformation that are usually applied to an alternating series on the form:</p>

<div class="math">$\sum_{n=0}^\infty (-1)^{n}a_n$</div>

<p>The Euler transformation of the series above is given below:</p>

<div class="math">$\sum_{n=0}^\infty \frac{\Delta^n a_0}{2^{n+1}}$</div>

<p>Here the general difference formula is defined as:</p>

<div class="math">$\Delta^n a_0 = \sum_{k=0}^n(-1)^k \left(\begin{array}{c}n \\ k \end{array}\right ) a_k$</div>

<p>What is actually happening here, and Euler&#39;s brilliant insight, is that the alternating series could be written out as:</p>

<div class="math">$S = \frac{a_0}{2} + \left(\frac{a_0}{2} - \frac{a_1}{2}\right) - \left(\frac{a_1}{2} - \frac{a_2}{2}\right) + \left(\frac{a_2}{2} - \frac{a_3}{2}\right) - \left(\frac{a_3}{2} - \frac{a_4}{2}\right) + \cdots $</div>

<p>By the way, you can always count on Euler beings astoundingly simple and brilliant. It there is a simple solution a problem&nbsp;he will find it in most cases. Doing the transformation on the alternating series multiple times effectively implies that we use higher and higher order of finite difference and getting a better and better approximation each time.</p>

<p>Actually the difference between Taylors formula and the&nbsp;finite difference equivalent can be illustrated by the exponential number <span class="math">\( e \)</span> from calculus and the difference equivalent 2:</p>

<div class="math">$2^x = 1 + x + \frac{x^{\underline{2}}}{2!} + \frac{x^{\underline{3}}}{3!}+ \frac{x^{\underline{4}}}{4!} + \frac{x^{\underline{5}}}{5!} + \cdots $</div>

<p>The underline in the power <span class="math">\(X^{\underline{2}}\)</span> is called falling powers:</p>

<div class="math">$x^{\underline{3}} = x (x-1)(x-2)$</div>

<p>The formula is quite general, you can even plug in rational numbers like <span class="math">\( 2^{\frac{1}{2}} = \sqrt{2}\)</span> and the formula still works. This was Newtons big insight into the binomial coefficients allowing him to calculate the square root of numbers.</p>

<p>Now for the code that performs the Euler transformation is given below:</p>

<pre lang="cs">
        /// &lt;summary&gt;
        /// Euler transform that transforms the alternating series 
        /// a_0 into a faster convergence with no negative coefficients
        /// &lt;/summary&gt;
        /// &lt;param name=&quot;a_0&quot;&gt;The alternating power series&lt;/param&gt;
        /// &lt;returns&gt;&lt;/returns&gt;
        public static double[] EulerTransformation(double[] a_0)
        {
            // Each series item
            List&lt;double&gt; a_k = new List&lt;double&gt;();

            // finite difference of each item
            double delta_a_0 = 0;

            for (int k = 0; k &lt; a_0.Length; k++)
            {
                delta_a_0 = 0;
                for (int m = 0; m &lt;= k; m++)
                {
                    double choose_k_over_m = (double)GetBinCoeff(k, m);
                    delta_a_0 += System.Math.Pow(-1d, (double)m) * 
                                 choose_k_over_m * System.Math.Abs(a_0[(int)(m)]);
                }

                a_k.Add(System.Math.Pow(1d / 2d, (double)(k + 1)) * delta_a_0);
            }
            return a_k.ToArray();
        }</pre>

<p>It is also possible to use Euler&#39;s series transformation on series that are not alternating, with the help of some conversion algorithm:</p>

<div class="math">$\sum_{n=1}^\infty v_n = \sum_{n=1}^\infty (-1)^{n-1} w_n$</div>

<p>where the altered series is the following transformation:</p>

<div class="math">$w_n = v_n + 2 \cdot v_{2n} + 4 \cdot v_{4n} + 8 \cdot v_{8n} \cdots$</div>

<h2>Levin&#39;s Transformation</h2>

<p>Levin&#39;s transformation is based on the assumption that one chooses b larger than zero and the remainder estimate <span class="math">\( \omega_m \)</span> is suitably chosen for the series that one wants to accelerate:</p>

<div class="math">$S_n = S - \omega_m \sum_{j=0}^\infty \frac{d_j}{(m+b)^j}$</div>

<p>The resulting formula is given below:</p>

<div class="math">$ \cal{L}_{k}^{(n)}(s)=\frac{\sum_{j=0}^{k}(-1)^{j} {k \choose j} c_{j,k,n}\frac{s_{n+j}}{a_{n+j+1}}}{\sum_{j=0}^{k}(-1)^{j} {k \choose j}c_{j,k,n}/a_{n+j+1}}$</div>

<p>The code implementation of the Levin algorithm:</p>

<pre lang="cs">
        /// &lt;summary&gt;
        /// Levins algorithm L_{k}^{n}(S_n) but k + n has to be less than 
        /// the length of S_n minus one
        /// &lt;/summary&gt;
        /// &lt;param name=&quot;S_n&quot;&gt;The partial truncated sums&lt;/param&gt;
        /// &lt;param name=&quot;k&quot;&gt;The &#39;order&#39; of Levins algorithm&lt;/param&gt;
        /// &lt;param name=&quot;n&quot;&gt;The n&#39;th point approximation&lt;/param&gt;
        /// &lt;returns&gt;The accelerated value&lt;/returns&gt;
        public static double LevinsAlgorithm(double[] S_n, double k = 4, double n = 5)
        {
            double numerator = 0;
            double denominator = 0;

            if (k + n &gt; S_n.Length - 1)
                throw new ArgumentException(&quot;Not enough elements to 
                      sum over with the n&#39;th element with order k&quot;);

            for (int j = 0; j &lt; n ; j++)
            {
                double rest = System.Math.Pow(-1d, j) * 
                              GetBinCoeff((double)k, (double)j);

                double C_jkn_U = System.Math.Pow((double)(n + j + 1), k - 1);
                double C_jkn_L = System.Math.Pow((double)(n + k + 1), k - 1);

                double C_njk = C_jkn_U / C_jkn_L;

                double S_nj = S_n[(int)n + (int)j];

                // t transform that calculates a_n
                double g_n = S_n[(int)n + j] - S_n[(int)n + j - 1];
                // u transform that calculates (n+k) * a_n
                // g_n = (n + k)* (S_n[(int)n + j] - S_n[(int)n + j - 1]);

                numerator += rest * C_njk * S_nj / g_n;
                denominator += rest * C_njk / g_n;
            }

            return numerator / denominator;
        }</pre>

<h2>References</h2>

<h3>Books</h3>

<ul>
	<li><a href="https://www.amazon.com/Advanced-Mathematical-Methods-Scientists-Engineers/dp/0387989315/ref=sr_1_1?keywords=advanced+mathematical+methods+for+scientists+and+engineers&amp;qid=1675954384&amp;sprefix=advanced+mathe%2Caps%2C163&amp;sr=8-1">Advanced Mathematical Methods Scientists Engineers</a></li>
	<li><a href="https://www.amazon.com/Asymptotic-Analysis-Perturbation-William-Paulsen/dp/1466515112/ref=sr_1_2?crid=4U8UOHF6EFQU&amp;keywords=asymtotic+analysis+and+perturbation+theory&amp;qid=1675954460&amp;sprefix=asymtotic+analysis+and+pertubation+theory%2Caps%2C138&amp;sr=8-2">Asymptotic Analysis Perturbation - William Paulsen</a></li>
	<li><a href="https://www.amazon.com/Numerical-Analysis-Historical-Developments-Century/dp/0444506179/ref=sr_1_1?crid=30IJJGQHZUTEO&amp;keywords=history+of+numerical+20th+century&amp;qid=1675954510&amp;sprefix=history+of+numerical+20th+century%2Caps%2C149&amp;sr=8-1">Numerical Analysis - Historical Developments in the 20<sup>th</sup> Century</a></li>
</ul>

<h3>Other Online Sources</h3>

<ul>
	<li><a href="https://download.uni-mainz.de/mathematik/Algebraische%20Geometrie/Euler-Kreis%20Mainz/E212%20Kapitel%201RevisedVersion2.pdf">Eulers acceleration transform</a></li>
	<li><a href="https://kconrad.math.uconn.edu/blurbs/analysis/series_acceleration.pdf">Accelerating Convergence of Series</a></li>
	<li><a href="https://dlmf.nist.gov/3.9">NIST Digital Library of Mathematical Functions</a></li>
	<li><a href="https://encyclopediaofmath.org/wiki/Euler_transformation">Euler transformation from Encyclopedia of Math - Springer</a></li>
</ul>

<h2>History</h2>

<ul>
	<li>9<sup>th</sup> February, 2023: Initial version</li>
</ul>

<p>&lt;quillbot-extension-portal&gt;</p>
