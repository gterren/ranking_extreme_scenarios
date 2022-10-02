# Identification of Extreme Scenarios in Day-ahead Operational Planning under Stocasticity

This repository contains the jupyter notebooks for the exploratory analysis to find a statistical method to indentify extreme scenarios in day-ahead forecasting. In particular, we are interest in finding extreme scenarios in solar, wind generation and energy demard for the operation planning of power networks. If it is possible to identify which scenarios are the most extreme, the computation cost of a operation contrained optimization problem under stochasticity can be simplify by only solving the problem for the most extreme scenarios. The method study in this reserach are included in the following notebooks:

## Kernel Principal Componenet Analsyis

**netload_extreme_scenarios.ipynb** analyzes the performance of Kernel Principal Componenent Analysis (KPCA) when applied to this problem.


## Functional Principal Componenet Analsyis

**depth_n_fPCA.R**

**sample_fPCA.R**


## Functional Depth

**extremality_metrics.ipynb** implements in Python multiple depth metrics. Some of the functions are from a library priviouly implemented in R. This notebook includes the state-of-the-art of depth notions from a recetly publish review. The depth notions aim to statistically define *outlierness* or *extremality*. See notebook write-ups for information and references.

**extremality_detection.ipynb** analyzes the performance of different depth notions rooted in functional data analysis methods. In addition, different heuristic are implemented and the performance are compared to depth notions. A method that exploits both the statistical and energy definition of extremality is proposed identify the most extrame scenarios when are defined as the most costly as a consecuence of the the error in the day-ahead operational planning.
