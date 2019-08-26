The MeltingSimulation_NewSetup contains the C++ code for the new setup and MeltingSimulation_OldSetup contains the C++ code for the old setup.

The Houdini folder contains the .hipnc file which shows the two different particle setup assets and the exported files for the finalRender, solidVsFluid and stiffSolid. The alembic files for these can also be found in the folder, as well as the two asset files. 

It should be possible to run the system by setting the simulation parameters in Houdini, right clicking on the asset and choosing save geometry. Then set the position of the file in the SimulationController and running the code. The code will export to the file specified in the SimulationController. 


