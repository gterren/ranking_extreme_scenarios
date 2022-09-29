# Identification Extreme Scenarios in Day-ahead Operational Planning under Stocasticity

This repository contains the jupyter notebooks for the exploratory analysis to find a statistical method to indentify extreme scenarios in day-ahead forecasting. In particular, we are interest in finding extreme scenarios in solar, wind generation and energy demard for the operation planning of power networks. If it is possible to identify which scenarios are the most extreme, the computation cost of a operation contrained optimization problem under stochasticity can be simplify by only solving the problem for the most extreme scenarios. The method study in this reserach are included in the following notebooks:

* $netload_extreme_scenarios.ipynb$ analyzes the performance of Kernel Principal Componenent Analysis (KPCA) when applied to this problem.
* implements in Python multiple depth metrics. Some of the functions are from a library priviouly implemented in R. This implemented the state-of-the-art of depth notions from a recetly publish review. See notebook write-ups for information and references.
* analyze the performance of different depth notions rooted in functional data analysis methods.
