DSC Visual Studio 2013

Information
------------
The DSC solution contains the following projects:
CGLA
DEMO
ISMesh
SOIL
src
SCGenerator
TetGen

Directions:
------------
1. Install Glew and Glut
2. Open the DSC.visualstudio.sln file with Visual Studio 2013
3. Go to DSC.visualstudio Solution Properties > Common Properties
4. Set 'Single startup project' as DEMO or SCGenerator depending on what you want to run 
5. Go to Demo Properties > Configuration Properties > Debugging
6. Set 'Configurations' to All Configurations
7. Set 'Working Directory' to the DSC folder which should be $(ProjectDir)\..\..\..\
