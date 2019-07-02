# The TECA, Toolkit for Extreme Climate Analaysis
![Storm Tracks Generated by TECA](doc/images/tracks_crop_2.gif)

TECA(Toolkit for Extreme Climate Analysis) is a collection of climate analysis algorithms geared toward extreme event detection and tracking implemented in a scalable parallel framework. The core is written in modern c++ and uses MPI+thread for parallelism. The framework supports a number of parallel design patterns including distributed data parallelism and map-reduce. Python bindings make the high performance c++ code easy to use. TECA has been used up to 750k cores.

[![Build Status](https://travis-ci.com/LBL-EESA/TECA.svg?token=zV3LhFtYvjcvo67W2uji&branch=master)](https://travis-ci.com/LBL-EESA/TECA)

# Documentation
For more information please see the [TECA User's Guide](https://teca.readthedocs.io/en/latest/).

#Copyright Notice#
TECA, Copyright (c) 2015, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.

NOTICE.  This software is owned by the U.S. Department of Energy.  As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, prepare derivative works, and perform publicly and display publicly.  Beginning five (5) years after the date permission to assert copyright is obtained from the U.S. Department of Energy, and subject to any subsequent five (5) year renewals, the U.S. Government is granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
