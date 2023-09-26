#! /usr/bin/env python3

import numpy as np
import pandas as pd
from Simulate.simulation_helpers.Sim_generator import pheno_simulator
from pandas_plink import read_plink1_bin, write_plink1_bin
import numpy as np

if __name__ == "__main__" :
    rng = np.random.default_rng(123)
    sim = pheno_simulator(rng = rng, nsubjects= 100, nSNPs = 100)
    sim.sim_sites(nsites =1)
    sim.sim_pops(nclusts= 2)
    sim.sim_genos()
    
    N = sim.df.shape[0]
    G = sim.genotypes.shape[1]
    
    temp = sim.genotypes.astype(str).T
    temp[temp=="0"] = "0/0"
    temp[temp=="1"] = "0/1"
    temp[temp=="2"] = "1/1"

    df = pd.DataFrame(temp, columns = ["sub" + str(i) for i in range(N)])
    df2 = pd.DataFrame({"CHROM" : np.repeat(1, G), "POS" : np.arange(G), "ID" : np.arange(G), "REF" : np.repeat("A", G), "ALT" : np.repeat("T", G),
                        "QUAL" : np.repeat(100, G), "FILTER" : np.repeat("PASS", G), "INFO" : np.repeat(".", G), "FORMAT" : np.repeat("GT", G)})
    
    pd.merge(df2, df, left_index = True, right_index = True).to_csv("toySim/toy.vcf", index = False, sep = '\t', header = True)
    sim.df.to_csv("toySim/ancestries.csv", index= False, header= True, sep = ",")

