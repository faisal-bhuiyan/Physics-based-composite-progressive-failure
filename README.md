# CFRP-KTF-Based-Progressive-Damage-Subroutines
This repository contains an implementation of KTF-based progressive damage modeling of continuous fiber-reinforced polymers under fatigue loding conditions. The user-defined material (UMAT) subroutines found here are compatible with the commercial FE code Abaqus. The work archived here was developed in the University of Wyoming at the [Fertig Research Group](http://www.uwyo.edu/mechanical/faculty-staff/Ray-Fertig/fertig.html) circa 2016 - 2021. Find more about the project [here](https://www.researchgate.net/project/A-Physics-Based-Combined-Creep-and-Fatigue-Prediction-Methodology-for-Fiber-Reinforced-Polymer-Composites) and consider citing them in your work if you find it useful.
  
## Instructions

1. Clone the project repository: `git clone https://github.com/faisal-bhuiyan/Physics-based-composite-progressive-failure.git`
2. Be familiar with the physics-based, multiscale progressive damage model based off KTF from the papers: [Bhuiyan 18](https://www.sciencedirect.com/science/article/pii/S0142112318303530) and [Bhuiyan 22](https://www.sciencedirect.com/science/article/pii/S0263822321011764). These two papers provide a comprehensive description of the material properties required by this routine, among other things.
3. Please read through the usage instructions on individual SUBROUTINES in this file as well as the included functions. This is super important to get an overview of the logic and how to modify the UMAT for your needs.
4. Define materials properties in `MODULE KTF_CONSTANTS` and `MODULE MATERIAL_PROPS`, in addition to the constituent and homogenized properties supplied via the text file called `MAT_FILE`. These properties can be computed using the provided Matlab routines, which require lamina level elastic properties and fatigue test data as inputs.
5. Run subroutine on the provided input files - especially the one element and the 140 element models to ensure things are working as intended. This code was not developed with good SW development practices in mind, so it can be fragile and painful to work with. Be patient while you're trying to run it for the first time.
6. Debugging in Fortran can require a lot of patience. Uncomment the existing `WRITE` statements (and add new ones) as required to print the debug statements throughout the UMAT as the first line of defence. Use debugging capabilities of the IDE of choice.
7. Because Github restricts large files in the repos, some of the larger input files to test the UMAT can be found [elsewhere](https://www.dropbox.com/scl/fo/i652cclt1f8b5o3gxb7j4/h?dl=0&rlkey=eqkpknacibwhk5h6xh4rr7faf).

## Directories

`src` - contains the source files required to run the UMAT. `subroutines` contain several versions of the same UMAT, name of the file indicates what is different in each one.

`input-files` - contains several input files to test the UMATs on.

`utils` - contains a text file containing a complete set of calibrated material properties as required by the material model.
