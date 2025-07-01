# Microbiome associations with riboflavin intermediates


## associations_mbiome.qmd

Raw ion counts are processed (relevant intermediates selected, blank subtraction, normalization, scaling).
Latent Dirichlet allocation (LDA) is performed for 9 topics, 
associations between the weights of the samples for these topics and the intermediates is then tested with 
linear models.
A differential abundance is performed on module (as defined in the Isala study) level (also see next paragraph). 


## mda.qmd

Script in which differential abundance at species level is performed using multiple methods.
Same procedure was followed on the CLI for the modules.
