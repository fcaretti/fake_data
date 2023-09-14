# fake_data
Repo to create data to simulate a deconvolution pipeline
A brief description of what each of the notebook does:
1. aml_data_deconvolution: plays with a few datasets found online
2. Clinical_Kernel_NB: given the statistics of the clinical data we will work with, it is a super simple way to create a patients-wise kernel of similarity 
3. test_notebook is an attempt to generate fae data for different patients. Notice that the choice of a Poissonian is probably wrong, as for RNA it seems better to use a Zero Inflated Negative Binomial (ZINB)
4. Run links fake RNA data to fake clinical data using simulated annealing

