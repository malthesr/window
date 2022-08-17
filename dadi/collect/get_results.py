import dadi
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import sys
sys.path.append("/home/genis/other/sfsMalthe/dadi_2022/scripts")
import dadi_models

indir="/home/genis/other/sfsMalthe/dadi_2022/results/fit"
outdir="/home/genis/other/sfsMalthe/dadi_2022/collect_results/output"


models = {"NOMIG":"no_mig",
         "IM": "sym_mig",
         "IM_ASYM": "asym_mig",
         "IM_ASYM_SIZE": "asym_mig_size"}


parslists = {"NOMIG":["N1", "N2", "T"],
         "IM": ["N1", "N2", "M", "T"],
         "IM_ASYM": ["N1", "N2", "M12", "M21", "T"],
         "IM_ASYM_SIZE": ["N1A", "N2A", "N1B","N2B", "M12", "M21", "T1", "T2"]}


sfs_versions={"shangani_masaimara_realsfs":"/home/genis/other/sfsMalthe/dadi_2022/results/2dsfs_dadi/shangani_masaimara_realsfs_unfolded_dadi.sfs",
              "shangani_masaimara_winsfs_2it":"/home/genis/other/sfsMalthe/dadi_2022/results/2dsfs_dadi/shangani_masaimara_winsfs_2it_unfolded_dadi.sfs"}


for modelname in models.keys():
    
    func = getattr(dadi_models, models[modelname]) # extract function defining model
    func_ex = dadi.Numerics.make_extrap_func(func)


    for sfs_version in sfs_versions.keys():

        print("doing model " + modelname + " for " + sfs_version)

        f = "{dd}/{p}/{m}/collected_pars_allruns.txt".format(dd=indir, p=sfs_version, m=modelname)
        sfs_path = sfs_versions[sfs_version]

        d = pd.read_table(f)
        d = d.sort_values(['likelihood'], ascending=False)

        fs = dadi.Spectrum.from_file(sfs_path)
        fs = fs.fold()
        fs.pop_ids = ["Shangani", "MasaiMara"]
        ns = fs.sample_sizes
        pts = [30,40,50]
        pars = parslists[modelname]
        npars = len(pars)


        ml_pars = list(d.iloc[0][1:-1])
        model = func_ex(ml_pars, ns, pts)
        model = dadi.Inference.optimally_scaled_sfs(model,fs)

        model = model.fold()
        masked_model, masked_data = dadi.Numerics.intersect_masks(model, fs)
        vmin = min(masked_model.min(), masked_data.min())

        resid = dadi.Inference.Anscombe_Poisson_residual(masked_model, masked_data,
                                                         mask=vmin)


        fs.to_file("{dd}/folded_observed_{p}.sfs".format(dd=outdir, p=sfs_version))
        model.to_file("{dd}/fitted_{m}_{p}.sfs".format(dd=outdir, m=modelname, p=sfs_version))
        resid.to_file("{dd}/residuals_{m}_{p}.sfs".format(dd=outdir, m=modelname, p=sfs_version))

        dadi.Plotting.plot_2d_comp_multinom(model, fs, show=False)
        plt.savefig("{dd}/fit_maxlike_{m}_{p}.png".format(dd=outdir, m=modelname, p=sfs_version))
        plt.clf()


        # reescale to get true params
        # quantities we need for the rescaling
        mu = 1.41e-08
        g = 5.7
        model = func_ex(ml_pars, ns, pts)
        theta = dadi.Inference.optimal_sfs_scaling(model, fs)
        L = fs.data.sum()


        rescaled_pars = ["Na"]
        [rescaled_pars.append(x) for x in pars]

        Na = theta / (4 * mu * L)
        rescaled_ml_pars = [Na]
        for i in range(npars):
            par = pars[i]
            ml_par = ml_pars[i]
            if par[0] == "N":
                p = ml_par * Na
            elif par[0] == "M":
                p = ml_par / (2*Na)
            elif par[0] == "T":
                p = ml_par * Na * 2 * g
            else:
                p = None
            rescaled_ml_pars.append(p)
            
        pars = ["Na"] + pars + ["likelihood"]
        rescaled_ml_pars.append(d.iloc[0]["likelihood"])
            
        with open("{dd}/rescaled_maxlike_estimates_{m}_{p}.txt".format(dd=outdir, m=modelname, p=sfs_version), "w+") as fh:
            fh.write(" ".join(pars) + "\n")
            fh.write(" ".join(map(str, rescaled_ml_pars)) + "\n")

