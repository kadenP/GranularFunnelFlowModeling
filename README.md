# GranularFunnelFlowModeling

A collection of files used to model the funnel flow process in particle-based high temperature thermal energy storage.

This repository contains three versions of the grannular funnel flow TES model:
1) The FF.m class definition, which contains all base definitions and functions used to simulate thermal storage cycles. This can be downloaded and used directly along with a run script (FF_test_1MW_baseline.m).
2) The TES.m class definition, which contains all base definitions with syntax and function definitions that are designed to integrate with code generation in Simulink. The base code in this class was adapted from FF.m with minor modifications.
3) GTS.mlapp is a packaged application that can be loaded directly into matlab and simulated with a user-friendly GUI.

It also contains a reduced class definition (Insulation.m) that focuses on the bin insulation (wall, roof, and base). This can be downloaded and simulated using a run script (InsulationTest_1MW.m).
