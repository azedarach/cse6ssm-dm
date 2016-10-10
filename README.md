# CSE6SSM Spectrum and Dark Matter

Analysis code for studying the mass spectrum and dark matter in the CSE6SSM,
used for the numerical analysis reported in [arXiv:1610.XXXXX][2-preprint].

## Requirements

To obtain analytic results using [SARAH][] requires

  * Mathematica (version 7.0 or higher)
  * SARAH (version 4.5.6 or higher)

The C++ spectrum generators were generated using [FlexibleSUSY][].  To compile
them requires

  * C++ compiler (g++ >= 4.7.2 or clang++ >= 3.1 or icpc >= 12.1)
  * Fortran compiler (gfortran, ifort)
  * [Boost][] (version 1.37.0 or higher)
  * [Eigen 3][] (version 3.1 or higher)
  * [GNU scientific library][] (GSL)
  * [Lapack / Blas][]
  
Calculation of the lightest CP-even Higgs mass using the EFT framework
requires

  * [SUSYHD][] version 1.0.2

To calculate dark matter predictions in the models requires

  * [micrOMEGAs][] version 4.1.8
  
The results reported in arXiv:1610.XXXXX were obtained using

  * SARAH-4.5.6 running on Mathematica 10.0
  * FlexibleSUSY-1.1.0
  * g++ and gfortran, version 4.8.2
  * Boost version 1.60.0
  * Eigen version 3.2.8
  * GSL version 1.15
  * SUSYHD-1.0.2
  * micrOMEGAs-4.1.8

[SARAH]: https://sarah.hepforge.org/   "SARAH"
[FlexibleSUSY]: https://flexiblesusy.hepforge.org/   "FlexibleSUSY"
[Boost]: https://www.boost.org   "Boost"
[Eigen 3]: https://eigen.tuxfamily.org   "Eigen 3"
[GNU scientific library]: https://www.gnu.org/software/gsl/   "GSL"
[Lapack / Blas]: http://www.netlib.org/lapack/   "Lapack"
[SUSYHD]: https://users.ictp.it/~susyhd/   "SUSYHD"
[micrOMEGAs]: https://lapth.cnrs.fr/micromegas/   "micrOMEGAs"

## Analytic Expressions

Analytic expressions for the mass matrices, loop corrections and
renormalization group equations can be generated using the
SARAH model files found under `./sarah`.  To run the SARAH model, a new model
must first be created in your SARAH installation, for example,

```shell
    $ mkdir /path/to/SARAH/Models/SE6SSM
    $ cp ./sarah/* /path/to/SARAH/Models/SE6SSM/
```

The model can then be run within a SARAH session using `Start["SE6SSM"];`.

## Mass Spectrum Calculation

Spectrum generators for the CSE6SSM and MSSM are located in the `./models`
directory.  The compilation procedure is the same as documented for an ordinary
FlexibleSUSY model.  The Makefile is generated and the code built by running

```shell
    $ ./configure --with-models=CSE6SSM,MSSM
    $ make
```

The generated executables are located in the directory for the corresponding
model.  See `./configure --help` for more options.

By default, only the two-scale spectrum generators are compiled.  To also
compile the semi-analytic solvers, the package should be configured using

```shell
    $ ./configure --with-models=CSE6SSM,MSSM --with-algorithms=semianalytic
```

Optionally, the lightest CP-even Higgs mass can be calculated using SUSYHD.
To do so, install the FlexibleSUSY addon [susyhd_call][] to the addons
directory (see the documentation provided with the addon for how to do this).
Executables calling SUSYHD can then be built using

```shell
    $ ./configure --with-addons=susyhd_call
```

After compiling with `make`, the executable `./models/CSE6SSM/run_susyhd_CSE6SSM.x`
can be used to compute the CP-even Higgs mass.  Note that the approximation
of using the MSSM EFT calculation is expected to be reasonable only when
the additional exotic contributions to the mass are small.  In general it
is advisable to use a model specific EFT calculation, see for example the
calculation provided by [FlexibleEFTHiggs][].

[susyhd_call]: https://github.com/dylan-harries/susyhd-call   "susyhd_call"
[FlexibleEFTHiggs]: https://flexiblesusy.hepforge.org/models.html#FlexibleEFTHiggs   "FlexibleEFTHiggs"

## Dark Matter Calculations

micrOMEGAs models for the CSE6SSM and CMSSM are contained in the directories 
`./micro-models/CSE6SSM` and `./micro-models/CMSSM`, respectively.  To run the 
dark matter calculation, the desired model should be added to your micrOMEGAs
installation and compiled.  The Makefiles for the models should first be generated
by configuring with microMEGAs enabled, using

```shell
   $ ./configure --enable-micromegas --with-micromegas-dir=/path/to/micromegas
```

The models can then be copied to your micrOMEGAs installation and used.  For
example, assuming the main micrOMEGAs package has already been compiled,

```shell
    $ cp -r ./micro-models/CSE6SSM /path/to/micromegas
    $ cd /path/to/micromegas/CSE6SSM
    $ make
```

After building the model, the executable `run_slha_file.x` is produced in
the `./src` sub-directory.  It reads an SLHA file and, by default, computes
the relic density for the point.  For example,

```shell
    $ ./src/run_slha_file.x LesHouches.out
```

See `./src/run_slha_file.x --help` for more options.

## References

If you use this code in your work please cite

  * P. Athron, D. Harries, R. Nevzorov, and A. G. Williams,
    [Phys. Lett. B760, 19 (2016)][1],
    [arXiv:1512.07040 \[hep-ph\]][1-preprint]
  * P. Athron, D. Harries, R. Nevzorov, and A. G. Williams, (2016),
    [arXiv:1610.XXXXX \[hep-ph\]][2-preprint]

The analytic results and micrOMEGAs model files were produced
using SARAH.  Please also cite the appropriate references for
SARAH,

  * F. Staub, [Comput. Phys. Commun. 181, 1077 (2010)][3],
    [arXiv:0909.2863 \[hep-ph\]][3-preprint]
  * F. Staub, [Comput. Phys. Commun. 182, 808 (2011)][4],
    [arXiv:1002.0840 \[hep-ph\]][4-preprint]
  * F. Staub, [Comput. Phys. Commun. 184, 1792 (2013)][5],
    [arXiv:1207.0906 \[hep-ph\]][5-preprint]
  * F. Staub, [Comput. Phys. Commun. 185, 1773 (2014)][6],
    [arXiv:1309.7223 \[hep-ph\]][6-preprint]

The C++ spectrum generators were produced using FlexibleSUSY, which
itself also depends on SARAH and components of SOFTSUSY.  If you use
this aspect of the code, please cite the above as well as

  * P. Athron, J.-h. Park, D. St&ouml;ckinger, and A. Voigt,
    [Comput. Phys. Commun. 190, 139 (2015)][7],
    [arXiv:1406.2319 \[hep-ph\]][7-preprint]
  * B. C. Allanach, [Comput. Phys. Commun. 143, 305 (2002)][8],
    [hep-ph/0104145][8-preprint]
  * B. C. Allanach, P. Athron, L. C. Tunstall, and A. G. Williams,
    [Comput. Phys. Commun. 185, 2322 (2014)][9],
    [arXiv:1311.7659 \[hep-ph\]][9-preprint]

If you use SUSYHD to compute the Higgs mass in the model, please cite

  * J. P. Vega and G. Villadoro, [JHEP 07, 159 (2015)][10],
    [arXiv:1504.05200 \[hep-ph\]][10-preprint]

If you use this code to calculate dark matter predictions for the models,
please cite the above as well as

  * G. B&eacute;langer, F. Boudjema, A. Pukhov, and A. Semenov,
    [Comput. Phys. Commun. 149, 103 (2002)][11],
    [hep-ph/0112278][11-preprint]
  * G. B&eacute;langer, F. Boudjema, A. Pukhov, and A. Semenov,
    [Comput. Phys. Commun. 174, 577 (2006)][12],
    [hep-ph/0405253][12-preprint]
  * G. B&eacute;langer, F. Boudjema, A. Pukhov, and A. Semenov,
    [Comput. Phys. Commun. 176, 367 (2006)][13],
    [hep-ph/0607059][13-preprint]
  * G. B&eacute;langer, F. Boudjema, A. Pukhov, and A. Semenov,
    [Comput. Phys. Commun. 180, 747 (2009)][14],
    [arXiv:0803.2360 \[hep-ph\]][14-preprint]
  * G. B&eacute;langer, F. Boudjema, P. Brun, A. Pukhov, 
    S. Rosier-Lees, P. Salati, and A. Semenov,
    [Comput. Phys. Commun. 182, 842 (2011)][15],
    [arXiv:1004.1094 \[hep-ph\]][15-preprint]
  * G. B&eacute;langer, F. Boudjema, A. Pukhov, and A. Semenov,
    [Comput. Phys. Commun. 185, 960 (2014)][16],
    [arXiv:1305.0237 \[hep-ph\]][16-preprint]
  * G. B&eacute;langer, F. Boudjema, A. Pukhov, and A. Semenov,
    [Comput. Phys. Commun. 192, 322 (2015)][17],
    [arXiv:1407.6129 \[hep-ph\]][17-preprint]

[1]: http://dx.doi.org/10.1016/j.physletb.2016.06.040   "Phys. Lett. B760, 19 (2016)"
[1-preprint]: https://arxiv.org/abs/1512.07040   "arXiv:1512.07040"
[2-preprint]: https://arxiv.org/abs/1610.XXXXX   "arXiv:1610.XXXXX"
[3]: http://dx.doi.org/10.1016/j.cpc.2010.01.011   "Comput. Phys. Commun. 181, 1077 (2010)"
[3-preprint]: https://arxiv.org/abs/0909.2863   "arXiv:0909.2863"
[4]: http://dx.doi.org/10.1016/j.cpc.2010.11.030   "Comput. Phys. Commun. 182, 808 (2011)"
[4-preprint]: https://arxiv.org/abs/1002.0840   "arXiv:1002.0840"
[5]: http://dx.doi.org/10.1016/j.cpc.2013.02.019   "Comput. Phys. Commun. 184, 1792 (2013)"
[5-preprint]: https://arxiv.org/abs/1207.0906   "arXiv:1207.0906"
[6]: http://dx.doi.org/10.1016/j.cpc.2014.02.018   "Comput. Phys. Commun. 185, 1773 (2014)"
[6-preprint]: https://arxiv.org/abs/1309.7223   "arXiv:1309.7223"
[7]: http://dx.doi.org/10.1016/j.cpc.2014.12.020   "Comput. Phys. Commun. 190, 139 (2015)"
[7-preprint]: https://arxiv.org/abs/1406.2319   "arXiv:1406.2319"
[8]: http://dx.doi.org/10.1016/S0010-4655(01)00460-X   "Comput. Phys. Commun. 143, 305 (2002)"
[8-preprint]: https://arxiv.org/abs/hep-ph/0104145   "hep-ph/0104145"
[9]: http://dx.doi.org/10.1016/j.cpc.2014.04.015   "Comput. Phys. Commun. 185, 2322 (2014)"
[9-preprint]: https://arxiv.org/abs/1311.7659   "arXiv:1311.7659"
[10]: http://dx.doi.org/10.1007/JHEP07(2015)159   "JHEP 07, 159 (2015)"
[10-preprint]: https://arxiv.org/abs/1504.05200   "arXiv:1504.05200"
[11]: http://dx.doi.org/10.1016/S0010-4655(02)00596-9   "Comput. Phys. Commun. 149, 103 (2002)"
[11-preprint]: https://arxiv.org/abs/hep-ph/0112278   "hep-ph/0112278"
[12]: http://dx.doi.org/10.1016/j.cpc.2005.12.005   "Comput. Phys. Commun. 174, 577 (2006)"
[12-preprint]: https://arxiv.org/abs/hep-ph/0405253   "hep-ph/0405253"
[13]: http://dx.doi.org/10.1016/j.cpc.2006.11.008   "Comput. Phys. Commun. 176, 367 (2007)"
[13-preprint]: https://arxiv.org/abs/hep-ph/0607059   "hep-ph/0607059"
[14]: http://dx.doi.org/10.1016/j.cpc.2008.11.019   "Comput. Phys. Commun. 180, 747 (2009)"
[14-preprint]: https://arxiv.org/abs/0803.2360   "arXiv:0803.2360"
[15]: http://dx.doi.org/10.1016/j.cpc.2010.11.033   "Comput. Phys. Commun. 182, 842 (2011)"
[15-preprint]: https://arxiv.org/abs/1004.1092   "arXiv:1004.1092"
[16]: http://dx.doi.org/10.1016/j.cpc.2013.10.016   "Comput. Phys. Commun. 185, 960 (2014)"
[16-preprint]: https://arxiv.org/abs/1305.0237   "arXiv:1305.0237"
[17]: http://dx.doi.org/10.1016/j.cpc.2015.03.003   "Comput. Phys. Commun. 192, 322 (2015)"
[17-preprint]: https://arxiv.org/abs/1407.6129   "arXiv:1407.6129"

## License

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
