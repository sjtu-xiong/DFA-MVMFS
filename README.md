# DFA-MVMFS
Detrended fluctuation analysis-based multivariate multifractal spectrum (DFA-MVMFS)

(1) bmc.m
Generates two BMC signals with probabilities  p1 and p2

(2) DFA_MVMFS.m
Performs DFA-BVMFS analysis on the two input signals, signal1 and signal2. The parameter m represents the polynomial fitting order, scale refers to the scale 𝑠 used in the DFA analysis, and q denotes the range of moments used, with the same moment range applied to both input signals.

(3) mylegendreM.m
Performs a bivariate Legendre transform using the scaling function and moment step size obtained from DFA-BVMFS and WL-BVMFS analyses to compute the bivariate multifractal spectrum.

(4) dwtleaders.m
MATLAB code for wavelet leaders computation in fractal analysis. References:
[1] Wendt, Herwig, and Patrice Abry. “Multifractality Tests Using Bootstrapped Wavelet Leaders.” IEEE Transactions on Signal Processing 55, no. 10 (October 2007): 4811–20. https://doi.org/10.1109/TSP.2007.896269.
[2] Jaffard, Stéphane, Bruno Lashermes, and Patrice Abry. “Wavelet Leaders in Multifractal Analysis.” In Wavelet Analysis and Applications, edited by Tao Qian, Mang I Vai, and Yuesheng Xu, 201–46. Basel: Birkhäuser Basel, 2007. https://doi.org/10.1007/978-3-7643-7778-6_17.

(5) MultiWL.m
Performs WL-BVMFS analysis on the two input signals, signal1 and signal2. The parameters Q and P represent the maximum values of the moments for the two inputs, and detaq specifies the step size.
