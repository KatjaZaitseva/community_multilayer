# Framework of community detection methods comparison for multilayer network data structures

The current version of the framework uses the algorithms available in R from the libraries [multinet](https://cran.r-project.org/web/packages/multinet/multinet.pdf) [1], [blockmodeling](https://cran.r-project.org/web/packages/blockmodeling/blockmodeling.pdf) [2] and SC-ML algorithm self-written based on the article [Clustering on Multi-Layer Graphs via Subspace Analysis on
Grassmann Manifolds](https://web.media.mit.edu/~xdong/paper/globalsip13.pdf) [3]. The framework provides code with functions. To get a comparison of algorithms for the specific network it is only necessary to load the network in a certain format, the details of which are specified in the file. The function will then output a modularity-based algorithms comparison table along with graphs, which can be used to determine which algorithm performs better on the loaded multilayer network. It is also possible to examine the performance time of the algorithms or their accuracy, if you already know the communities of the network.

The majority of the code is made with the most exhaustive library available for multilayer structures in R - multinet, details of which can be found in the following [repository](https://bitbucket.org/uuinfolab/r_multinet/src/master/). 

Below is an aggregation of available libraries with algorithms for __R__ and __Python__:
* R
* Python
  * [Community detection with node attributes in multilayer networks](https://github.com/mcontisc/MTCOV/tree/master/code) [4]
  * [A multilevel clustering technique for community detection](https://github.com/ijdutse/mct) [5]

- - - -

References:

[1] Dickison, Magnani, and Rossi, 2016. Multilayer Social Networks. Cambridge University Press. ISBN: 978-1107438750

[2] Žiberna, A. (2014). Blockmodeling of multilevel networks. Social Networks, 39(1), 46-61. doi: 10.1016/j.socnet.2014.04.002

[3] Dong, X., Frossard, P., Vandergheynst, P., & Nefedov, N. (2013). Clustering on multi-layer graphs via subspace analysis on Grassmann manifolds. IEEE Transactions on signal processing, 62(4), 905-918.

[4] Contisciani M., Power E. & De Bacco C. (2020). Community detection with node attributes in multilayer networks, Scientific Reports 10, 15736 (2020)

[5] Isa Inuwa-Dutse, Mark Liptrott, Ioannis Korkontzelos, A multilevel clustering technique for community detection, Neurocomputing, Volume 441, 2021, Pages 64-78, ISSN 0925-2312, doi: 10.1016/j.neucom.2021.01.059
