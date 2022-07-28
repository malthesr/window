"""
Rules for getting the log-likelihood of the estimating SFS (and tidying these)
"""


include: "log_likelihood/iterations.smk"
include: "log_likelihood/blocks.smk"


rule all_log_likelihood:
    """ Phony rule to expand log-likelihoods based on config. """
    input:
        expand(rules.gather_log_likelihoods.output.csv, csv="loglik", **config["run"]),


rule all_block_log_likelihood:
    """ Phony rule to expand log-likelihoods based on config. """
    input:
        expand(rules.gather_block_log_likelihoods.output.csv, **config["run"]),
