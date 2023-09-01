# *i*MFP2023 model validation and analysis

The notebooks in this repository recreate the analyses from sections 2.4.1, 2.4.2, and 2.4.3.
We used a combination of bibliomic and experimental data to validate model *i*MFP2023 in three conditions:
- Methanotrophic
- Heterotrophic (Propape, Isopropanol, Acetone)
- Autotrophic (CO<sub>2</sub>+H<sub>2</sub>)

Although notebooks can run independetly, we recommend following the proposed workflow

---
**Important note:** We suggest using jupyter-notebook instead of jupyter-lab to execute the notebooks.
Jupyterlab 3 and the most recent versions of IPwidgets are incompatible with Escher maps. Nonetheless, for users of jupyter-lab, we provide a yml file with a list of the compatible dependencies used in the creation of the notebooks.

To install the dependencies in a new environment run these commands in a terminal:

`conda env create -f environment.yml`

`conda activate saldivar_etal_2023`

`jupyter labextension install @jupyter-widgets/base && jupyter labextension install escher && jupyter labextension install @jupyter-widgets/jupyterlab-manager`

---

## 1. Sensitivity analysis to Growth and Non-Growth Associated Maintenance

Run the following notebook to generate Figure S3 and calculate the sensitivity of the growth rate to changes in GAM and NGAM

[1-sensitivity_analysis_to_gam_ngam.ipynb](1-sensitivity_analysis_to_gam_ngam.ipynb)

## 2. Effect of reaction FALDHpp on methanotrophic growth

The next notebook generates Figure S4 and Figure 3A, B, as well as table 2. In this section we study the changes in growth parameters as a function of flux through reaction FALDHpp

[2-simulations_methane.ipynb](2-simulations_methane.ipynb)

## 3. Effect of reaction HYD4pp on autotrophic growth

The third notebook generates Figure 3C and Figure 4. In this section we study the trade-off between using hydroplasmic hydrogenases vs. cytoplasmic hydrogenases for hydrogen oxidation.

[3-simulations_h2.ipynb](3-simulations_h2.ipynb)

## 4. Metabolic changes on heteretrophic growth

This section is divided in three notebooks:
- Notebook 4A. Generates Figure 3B, and calculates growth parameters in heterotrophic conditions

[4A-simulations_c3.ipynb](4A-simulations_c3.ipynb)

- Notebooks 4B and 4C are the only two that need to be run sequentially
- Notebook 4B. Runs flux sampling and analysis of the results. It produces files summary_sampling.xlsx, pagerank_Acetone.csv, pagerank_Methane.csv, pagerank_Propane.csv, pagerank_Isopropanol.csv

[4B-flux_sampling_simulations.ipynb](4B-flux_sampling_simulations.ipynb)

- Notebook 4C. Generates Figures 5A,B,C and Figure 5D.

[4C-flux_sampling_results.ipynb](4C-flux_sampling_results.ipynb)
