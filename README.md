# 9N WUSTL Sort Code

WUSTL sort code for use with Gobbi 28 array for measurement of 9N->5p+α and 12O->4p+2α @ FRIB. This particular Gobbi variant is Si-CsI, and has had the innermost CsI quadrants split into four smaller CsI crystals for 7 crystals per telescope and 28 total. I adapted this code from a previous code I wrote for a 6Li(2+2)->α+p+n experiment @ TAMU, before which I had adapted from Nicolas Dronchi, who adapted his code from Robert Charity and Kyle Brown. So continues the passing down of the Wash U Radiochemistry sort code. All is as it should be.

Writen by Henry Webb (h.s.webb@wustl.edu)

# Instructions

To set up project:

1. Run `git clone <URL>`, where `<URL>` should be replaced with the desired URL of the project's GitHub repository
2. Run `cd n9sort` to enter the project directory
3. Proceed to the build steps below

To build:

1. In root project directory, run `mkdir build`
2. Then `cd build`
3. Then `cmake cmake -DROOT_DIR=$ROOTSYS/cmake -DCMAKE_C_COMPILER=/bin/gcc -DCMAKE_CXX_COMPILER=/bin/g++ ..`
	- Make sure to replace the paths to gcc and g++ with the desired compiler, or ommit entirely if your system does not have multiple potentially conflicting compilers
4. Finally, `make -j`
5. Run the code with `./sort` inside the build directory

Before running the code, make sure to change all of the configuration settings in `sort.config` to the desired values, along with all the other various required input files for the different parts of the code (CsI channel map, Si strip ranges, calibrations, etc.).

To run supplimentary analysis macros:

1. I usually put macros inside the `macros/` directory of the project. Examples include `O12_4p2a_process.C` and `N9_5pa_process.C`.
2. These kinds of supplimentary analysis macros usually require the use of classes and dictionaries generated for/from the base sort code. For this to work, one must run the macros from inside the build directory AFTER building the base sort code, next to where the `.so` shared library and ROOT dictionary files are placed. This looks something like `root -l -q "../macros/N9_5pa_process.C+"`. The extra `+` enables ACLiC, which I believe means that ROOT uses g++ or whatever compiler to parse the macro instead of using the cling interpreter like it otherwise would.

# CsI input files

CsIChannelMap.dat column ordering:
* Each CsI crystal across the whole Gobbi 28 array (all four telescopes) gets its own line here
* Column 1 (`chan`) is the ADC channel, starting from 0
* Column 2 (`tel`) is the telescope ID, starting from 0
* Column 3 (`id`) is the crystal ID within the telescope (ranging from 0 to 6 for Gobbi 28)

```
chan | tel | id
```

CsIStripExtents.dat column ordering:
* Each CsI crystal within a single telescope gets its own line here, and all telescopes are then assumed to be the same and rotationally invariant
* The lines in the file are read in order of crystal ID, from 0 to 6, but are not directly labeled as such
* Column 1 (`Fmin`) is the minimum front silicon strip corresponding to the silicon detector area in front of the crystal
* Column 2 (`Fmax`) is the maximum front silicon strip corresponding to the silicon detector area in front of the crystal
* Column 3 (`Bmin`) is the minimum back silicon strip corresponding to the silicon detector area in front of the crystal
* Column 4 (`Bmax`) is the maximum back silicon strip corresponding to the silicon detector area in front of the crystal

```
Fmin | Fmax | Bmin | Bmax
```

# Old instructions from Nic

zlines:
To draw zlines, open sort.root, navigate to the DeltaE-E plots (DEEplots/DEE0).
Run the following lines in the root command line:

.L zline/banana.C
banana()

//here you need to click on the TBrowser interface to draw the gates

//Copy paste the data from banana.dat into a .zline file for the particle identified

.L zline/readPIDzline.C
readPIDzline(1)
//where 1 is the quadrant number of interest

**Note from Henry:** There is also a macro somewhere in this repository which we used for CsI calibrations during the experiment, I think it's like `Cal/CsICal.C` or something. This macro allows you to calibrate easily by just clicking on the histogram to set the centroid, which it uses to automatically perform a peak fit. This is an older macro which I think was written by Kyle or something.
