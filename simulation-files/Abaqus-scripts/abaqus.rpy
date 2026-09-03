# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2024 replay file
# Internal Version: 2023_09_21-08.55.25 RELr426 190762
# Run by cwei2 on Fri Aug 21 10:33:28 2026
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=225.668411254883, 
    height=95.4980316162109)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
import os
os.chdir(r"/groups/mhollan5/commandlab/Chongyi/Repositories/Test_astrocytes")
execfile(
    '/groups/mhollan5/commandlab/Chongyi/Repositories/Test_astrocytes/QuarterEllipse_PUSH_script.py', 
    __main__.__dict__)
#: The model "QuarterEllipse_pull_gammahat5em2" has been created.
#: Warning: findAt could not find a geometric entity at (17.25, 24.6817240078565, 0.5)
#: The interaction property "IntProp-1" has been created.
#: The interaction property "IntProp-2" has been created.
#: The interaction "Int-1" has been created.
#: The interaction "Int-2" has been created.
#: The interaction "Int-3" has been created.
#: Warning: findAt could not find a geometric entity at (12.097431, 26.376844, 0.0)
#: Warning: findAt could not find a geometric entity at (31.092919, 11.744441, 1.0)
#: Warning: findAt could not find a geometric entity at (34.2, 0.0, 0.25)
#: Warning: findAt could not find a geometric entity at (0.0, 28.2, 0.75)
#: Warning: parallelizationMethodExplicit member has been removed from Job objects.
#: Job QuarterEllipse_pull_gammahat5em2: Analysis Input File Processor completed successfully.
#: Job QuarterEllipse_pull_gammahat5em2: Abaqus/Explicit Packager completed successfully.
#: Job QuarterEllipse_pull_gammahat5em2 crashed
#: The model "QuarterEllipse_pull_gammahat10em2" has been created.
#: Warning: findAt could not find a geometric entity at (17.25, 24.6817240078565, 0.5)
#: The interaction property "IntProp-1" has been created.
#: The interaction property "IntProp-2" has been created.
#: The interaction "Int-1" has been created.
#: The interaction "Int-2" has been created.
#: The interaction "Int-3" has been created.
#: Warning: findAt could not find a geometric entity at (12.097431, 26.376844, 0.0)
#: Warning: findAt could not find a geometric entity at (31.092919, 11.744441, 1.0)
#: Warning: findAt could not find a geometric entity at (34.2, 0.0, 0.25)
#: Warning: findAt could not find a geometric entity at (0.0, 28.2, 0.75)
#: Warning: parallelizationMethodExplicit member has been removed from Job objects.
#: Job QuarterEllipse_pull_gammahat10em2: Analysis Input File Processor completed successfully.
#: Job QuarterEllipse_pull_gammahat10em2: Abaqus/Explicit Packager completed successfully.
