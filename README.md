# BDMsolver

This is a BDM finite-element implementation of an ice-sheet model. The ice-sheet model is a full Stokes model with Glen's flow law and sliding boundary conditions.
The impenetrability condition is implemented as dirichlet condition in a BDM_1 finite element space. For the solving algorithm the interior penalty method is applied.
This model is validated by performing test E of the ISMIP-HOM benchmark experiments and comparing the results against the ones obtained by a P2P1 finite-element solver of the full Stokes model. The here uploaded files are the ingredients to perform the plots of velocity and basal shear stress of the BDM solver and the P2P1 solver.

The full Stokes ice-sheet model is solved in BDMsolver_Experiment1.py and BDMsolver_Experiment2.py for two different basal sliding scenarios. In in the first file no basal sliding occurs, hence, the slipperiness factor beta is set very high. In the second file partial basal sliding occures.

The plots for velocity and basal shear stress can be performed with the files Plot_velocity.py and Plot_TauMeshStudy.py

The mesh files are available in .geo and also already converted to .xml.
The files oga1e000.txt and oga1e001.txt obtain the solution of the test E1/E2 of the ISMIP-HOM benchmark experiments obtained with the P2P1 solver.
