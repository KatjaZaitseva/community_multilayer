# Framework of community detection methods comparison for multilayer network data structures

The current version of the framework uses the algorithms available in R from the libraries [multinet](https://cran.r-project.org/web/packages/multinet/multinet.pdf) [1], [blockmodeling](https://cran.r-project.org/web/packages/blockmodeling/blockmodeling.pdf) [2] and [multilayer_extraction](https://github.com/jdwilson4/MultilayerExtraction) [3]. The framework provides code with functions in which it is only necessary to load the network in a certain format, the details of which are specified in the file. The function will then output a modularity-based algorithm comparison table along with graphs, which can be used to determine which algorithm performs better on the loaded multilayer network.

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

[3] Wilson, J.D., Palowitch, J., Bhamidi, S., and Nobel, A.B. (2017) Significance based extraction in multilayer networks with heterogeneous community structure, Journal of Machine Learning Research (18) 1-49

[4] Contisciani M., Power E. & De Bacco C. (2020). Community detection with node attributes in multilayer networks, Scientific Reports 10, 15736 (2020)

[5] Isa Inuwa-Dutse, Mark Liptrott, Ioannis Korkontzelos, A multilevel clustering technique for community detection, Neurocomputing, Volume 441, 2021, Pages 64-78, ISSN 0925-2312, doi: 10.1016/j.neucom.2021.01.059
