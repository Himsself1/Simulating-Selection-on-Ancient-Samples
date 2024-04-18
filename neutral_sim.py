# * Libraries

import sys
import math
from pathlib import Path
import msprime
import numpy
import tskit
import argparse
import random
import time

# import matplotlib.pyplot as plt
from more_itertools import collapse
import gc
import threading
from IPython.display import SVG, display

# * Starting Variables

gen_time = 28

Starting_N_Anatolians = 14516
N_WHG = 2000
N_Steppe = 6500
N_Parent_Anat_WHG = 2340
N_Parent_All = 2340
N_Parent_pre_WHG = 12000

# * Figure out Samples

primitive_samples = numpy.genfromtxt(
    "/home/himsself/Documents/Simulating-Selection-on-Ancient-Samples/Distribution_of_retained_samples_by_gen_to_plot_no_single_gens",
    delimiter=" ",
)[1:, :]

list_of_samples = [
    msprime.SampleSet(
        num_samples=int(i[1]),
        population="Anatolians",
        time=i[0],
    )
    if i[0] * 28 < 4500
    else msprime.SampleSet(
        num_samples=int(i[1]),
        population="Anat_Pre_steppe_admix",
        time=i[0],
    )
    for i in primitive_samples
]

# * Argument for running multiple times

cli = argparse.ArgumentParser(
    description="msprime wrapper for simulations", prog="neutral_sim.py"
)
cli.add_argument(
    "-times", type=int, required=True, help="How many replicated will be simulated"
)
cli.add_argument(
    "-out_folder", type=str, required=True, help="Full path to output folder"
)
argue = cli.parse_args()

# ** Make Names
Path(argue.out_folder).mkdir(parents=True, exist_ok=True)  # Makes Output directories
out_names = [
    argue.out_folder + "/neutral_sim_output_" + str(i) + ".vcf"
    for i in range(argue.times)
]

# * Adding Populations and Events

for i in range(argue.times):
    demography = msprime.Demography()
    demography.add_population(name="Anatolians", initial_size=Starting_N_Anatolians)
    demography.add_population(name="WHG", initial_size=N_WHG)
    demography.add_population(name="Steppe", initial_size=N_Steppe)
    demography.add_population(
        name="Anat_Pre_steppe_admix", initial_size=N_Parent_pre_WHG
    )
    demography.add_population(name="Anat_Pre_WHG_admix", initial_size=N_Parent_pre_WHG)
    demography.add_population(name="Parent_Anat_WHG", initial_size=N_Parent_Anat_WHG)
    demography.add_population(name="Parent_All", initial_size=N_Parent_All)
    demography.add_population_split(
        time=36999 / gen_time,
        derived=["WHG", "Anat_Pre_WHG_admix"],
        ancestral="Parent_Anat_WHG",
    )
    demography.add_population_split(
        time=37000 / gen_time,
        derived=["Parent_Anat_WHG", "Steppe"],
        ancestral="Parent_All",
    )
    demography.add_admixture(
        time=4500 / gen_time,
        derived="Anatolians",
        ancestral=["Anat_Pre_steppe_admix", "Steppe"],
        proportions=[0.66, 0.34],
    )
    demography.add_admixture(
        time=8500 / gen_time,
        derived="Anat_Pre_steppe_admix",
        ancestral=["Anat_Pre_WHG_admix", "WHG"],
        proportions=[0.5, 0.5],
    )
    demography.add_population_parameters_change(
        time=95800 / gen_time, initial_size=29100, population="Parent_All"
    )
    demography.sort_events()
    ts = msprime.sim_ancestry(
        demography=demography,
        samples=list_of_samples,
        recombination_rate=1e-8,
        sequence_length=1e6,
        ploidy=2,
    )
    mutated_ts = msprime.sim_mutations(
        ts, rate=1e-8, model="binary", discrete_genome=True, keep=False
    )
    f = open(out_names[i], "w")
    mutated_ts.write_vcf(f)
    f.close()
    print("Simulation " + str(i) + " finished")
