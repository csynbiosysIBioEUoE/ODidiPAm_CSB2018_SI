## Optimally designed vs intuition-driven inputs: The study case of promoter activity modelling ##


The scripts used to generate, analyse and visualise the data presented in the paper can be found in this repository.To run these scripts, you will need AMIGO2 Toolbox, available at the link [https://sites.google.com/site/amigo2toolbox/download](https://sites.google.com/site/amigo2toolbox/download)

The scripts are organised in multiple subfolders:

- --MPLac contains:
  - Run\_MPLac\_in\_silico\_experiments.m, which calls MPLac\_in\_silico.m. The latter runs a simulation of the MPLac model presented by Gnugge et al. [1] in response to a step in IPTG. The same script is used to generate a dose-response curve;
  - MPLac\_load\_model.m and MPLac\_Compute\_SteadyState.m are used to specify the model (structure and parameter values) and compute the analytical steady state for a specified parameter vector and IPTG concentration;
  - MPLac\_SSE\_Computation.m runs a simulation of the model under the experimental conditions considered in [1] and computes the SSE between numerical and experimental expression levels.
- --MPLacr contains:
  - Run\_MPLacr\_in\_silico\_experiment.m, which calls MPLacr\_in\_silico.m. The latter runs a simulation of the MPLacr model (obtained by re-fitting the model structure presented by Gnugge et al. [1]  to the experimental data acquired by the same authors) in response to a step in IPTG. The same script is used to generate a dose-response curve;
  - PLac\_load\_model.m and PLac\_Compute\_SteadyState.m are used to specify the model (structure and parameter values) and compute the analytical steady state for a specified parameter vector and IPTG concentration;
  - MPLacr\_SSE\_Computation.m runs a simulation of the model under the experimental conditions considered in [1] and computes the SSE between numerical and experimental expression levels.
  - Run\_fit\_to\_MPLacr.m calls fit\_to\_MPLacr.m, which is used to perform parameter estimation, using cross validation and starting from 100 different initial conditions, on the experimental datasets available in [1].
  - ExtractionOfBestEstimate.m extracts, from the 100 iterations, the vector of parameter estimates yielding the minimum SSE over the test set.
- --MIP contains:
  - Run\_MIP\_in\_silico\_experiment.m, which calls MIP\_in\_silico.m. The latter runs a simulation of the MIP model (obtained by fitting the 3D model structure to the experimental data acquired in [1]) in response to a step in IPTG. The same script is used to generate a dose-response curve;
  - M3D\_load\_model\_Experimental.m and M3D\_steady\_state\_Experimental.m are used to specify the model (structure and parameter values) and compute the analytical steady state for a specified parameter vector and IPTG concentration;
  - MIP\_SSE\_Computation.m runs a simulation of the MIP model under the experimental conditions considered in [1] and computes the SSE between numerical and experimental expression levels.
  - Run\_fit\_to\_MIP.m calls fit\_to\_MIP.m, which is used to perform parameter estimation, using cross validation and starting from 100 different initial conditions, on the experimental datasets available in [1].
  - ExtractionOfBestEstimate.m extracts, from the 100 iterations, the vector of parameter estimates yielding the minimum SSE over the test set.
  - Ref7\_Data\_StructureCreation.m creates, starting from data provided by the authors in [1], a structure to be used in the fitting. Please refer to reference [1] for access to the data.
- --MIPr contains:
  - Run\_fit\_to\_MIPr.m calls fit\_to\_MIPr.m, which is used to perform parameter estimation, using cross validation and starting from 100 different initial conditions, on the pseudo experimental dataset.
  - M3D\_load\_model.m and M3D\_steady\_state.m are used to specify the model (structure and parameter values) and compute the analytical steady state for a specified parameter vector and IPTG concentration;
  - ExtractionOfBestEstimate.m extracts, from the 100 iterations, the vector of parameter estimates yielding the minimum SSE over the test set.
  - ExtractionOfDataForImage3.m simulates the response of the fitted MIPr to inputs prototype from the pseudo-data datasets and extracts the csv files required to generate image 3 of the paper.
  - bounds_for_parameter_estimation.docx specifies the bounds on parameters during parameter estimation.
- --InputComparison contains:
  - Run\_in\_silico\_experiment\_parfor\_IntuitionDriven.m, which calls (fit\_to\_InduciblePromoter\_Step/Pulse/Random.m). These scripts simulate the response of MIPr to the specified classes of inputs and runs parameter estimation using 100 different initial conditions.
  - Run\_in\_silico\_experiment\_parfor\_Optimised.m, which calls (fit\_to\_InduciblePromoter\_Optimised.m). This script runs OID, simulates the response of MIPr to the designed input and runs parameter estimation. The procedure is repeated starting from 100 different initial conditions.  Note that the same script is used to generate the data on which the comparison between on-line and off-line OID is based: you just need to specify the NumLoops you wish to execute. 
  - InduciblePromoter\_load\_model.m and InduciblePromoter\_load\_model\_optimised.m load the model used for the intuition driven and 	optimised inputs. The models differ because of the presence of a constraint on the input design. While in here the constraint ensures non-negative expression levels, the script was organised in this way to ensure flexibility of use for in-vivo experiments (in that case the constraint would support the selection of input values compatible to a detectable expression level, for example).
  - InduciblePromoter\_steady\_state.m computes the analytical steady state of MIPr for the specified parameter vector and IPTG concentration.
- --Images contains:
  - MPLac\_MPLacr\_MIP\_comparison.ipynb required to generate Fig. 2 of the paper.
  - MIPr\_Fitting\_Fig3.ipynb used to generate Fig. 3 of the paper
  - MIPr\_InputComparison\_ParameterUncertainty.ipynb used to generate Fig.4 of the paper.
  - Analysis_OIDLoops.ipynb is used to generate Fig. 6 of the paper. 

To run the scripts in AMIGO2, please copy the folder of interest in the folder Examples.

The data associated with these scripts can be found at:

[https://datasync.ed.ac.uk/index.php/s/xGL4oSbJMyu91Bt](https://datasync.ed.ac.uk/index.php/s/xGL4oSbJMyu91Bt)
(Pwd: ODidiPAm\_CSB2018\_Data)

References:

[1] Gnügge, R., Dharmarajan, L., Lang, M. and Stelling, J., 2016. An Orthogonal Permease–Inducer–Repressor Feedback Loop Shows Bistability. _ACS synthetic biology_, _5_(10), pp.1098-1107.
