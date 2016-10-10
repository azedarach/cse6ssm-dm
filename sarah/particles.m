ParticleDefinitions[GaugeES] = {
    {SdL,      { Description -> "Left Down-Squarks"}},
    {SdR,      { Description -> "Right Down-Squarks"}},
    {SuL,      { Description -> "Left Up-Squarks"}},
    {SuR,      { Description -> "Right Up-Squarks" }},
    {SeL,      { Description -> "Left Selectron"}},
    {SeR,      { Description -> "Right Selectron"}},
    {SvL,      { Description -> "Left Sneutrino"}},
    {SHd0,     { Description -> "Neutral Down-Higgs"}},
    {SHdm,     { Description -> "Charged Down-Higgs"}},
    {SHu0,     { Description -> "Neutral Up-Higgs"}},
    {SHup,     { Description -> "Charged Up-Higgs"}},
    {SsR,      { Description -> "Singlet"}},
    {SsbarR,   { Description -> "Sbar SM Singlet",
                 PDG -> {0},
                 PDG.IX -> 0,
                 Width -> 0,
                 FeynArtsNr -> 668,
                 LaTeX -> "\\overline{S}",
                 OutputName -> "SsbR"}},
    {SphiR,    { Description -> "Pure Singlet",
                 PDG -> {0},
                 PDG.IX -> 0,
                 Width -> 0,
                 FeynArtsNr -> 669,
                 LaTeX -> "\\phi",
                 OutputName -> "SphiR"}},
    {SDxL,     { Description -> "Left SExotics",
                 PDG -> {0, 0, 0},
                 PDG.IX -> {0, 0, 0},
                 FeynArtsNr -> 666,
                 LaTeX -> "\\tilde{Dx}_L",
                 OutputName -> "SDxL"}},
    {SDxbarR,  { Description -> "Right SExotics",
                 PDG -> {0, 0, 0},
                 PDG.IX -> {0, 0, 0},
                 FeynArtsNr -> 667,
                 LaTeX -> "\\tilde{Dx}_R",
                 OutputName -> "SDxbR"}},
    {SH1I0,    { Description -> "Neutral Inert-Down-Higgs",
                 PDG -> {0, 0},
                 PDG.IX -> {0, 0},
                 FeynArtsNr -> 101,
                 LaTeX -> "h^{0Inert}_{1}",
                 OutputName -> "SH1I0"}},
    {SH1Im,    { Description -> "Charged Inert-Down-Higgs",
                 PDG -> {0, 0},
                 PDG.IX -> {0, 0},
                 FeynArtsNr -> 102,
                 LaTeX -> "h^{-Inert}_{1}",
                 OutputName -> "SH1IM"}},
    {SH2I0,    { Description -> "Neutral Inert-Up-Higgs",
                 PDG -> {0, 0},
                 PDG.IX -> {0, 0},
                 FeynArtsNr -> 103,
                 LaTeX -> "h^{0Inert}_{2}",
                 OutputName -> "SH2I0"}},
    {SH2Ip,    { Description -> "Charged Inert-Up-Higgs",
                 PDG -> {0, 0},
                 PDG.IX -> {0, 0},
                 FeynArtsNr -> 104,
                 LaTeX -> "h^{+Inert}_{2}",
                 OutputName -> "SH2IP"}},
    {SSIR,     { Description -> "Inert-Singlet",
                 PDG -> {0, 0, 0},
                 PDG.IX -> {0, 0, 0},
                 FeynArtsNr -> 105,
                 LaTeX -> "S^{I}",
                 OutputName -> "SI"}},
    {SHpd0,    { Description -> "Neutral L4",
                 PDG -> {0},
                 PDG.IX -> {0},
                 FeynArtsNr -> 900,
                 LaTeX -> "L_4^0",
                 OutputName -> "SHpd0"}},
    {SHpdm,    { Description -> "Negative Charged L4",
                 PDG -> {0},
                 PDG.IX -> {0},
                 FeynArtsNr -> 901,
                 LaTeX -> "L_{4}^{-}",
                 OutputName -> "SHpdm"}},
    {SHpup,    { Description -> "Positive Charged L4-Bar",
                 PDG -> {0},
                 PDG.IX -> {0},
                 FeynArtsNr -> 902,
                 LaTeX -> "\\bar{L}_{4}^{+}",
                 OutputName -> "SHpup"  }},
    {SHpu0,    { Description -> "Neutral L4-Bar",
                 PDG -> {0},
                 PDG.IX -> {0},
                 FeynArtsNr -> 903,
                 LaTeX -> "\\bar{L}_{4}^{0}",
                 OutputName -> "SHpu0"}},
    {VB,       { Description -> "B-Boson"}},
    {VG,       { Description -> "Gluon"}},
    {VWB,      { Description -> "W-Bosons"}},
    {VBp,      { Description -> "B'-Boson"}},
    {gB,       { Description -> "B-Boson Ghost"}},
    {gG,       { Description -> "Gluon Ghost" }},
    {gWB,      { Description -> "W-Boson Ghost"}},
    {gBp,      { Description -> "B'-Boson Ghost",
                 LaTeX -> "\\eta^{B^\\prime}"}},
    {Glu,      { Description -> "Gluino"}},
    {Wino,     { Description -> "Wino"}},
    {Bino,     { Description -> "Bino"}},
    {FBp,      { Description -> "Bino'"}},
    {H01,      { Description -> "Neutral down Higgsinos",
                 PDG -> {0},
                 PDG.IX -> {0},
                 FeynArtsNr -> 122,
                 LaTeX -> "\\tilde{H}^0_1",
                 OutputName -> "H01"}},
    {HC1,      { Description -> "Charged down Higgsinos",
                 PDG -> {0},
                 PDG.IX -> {0},
                 FeynArtsNr -> 123,
                 LaTeX -> "\\tilde{H}^-_1",
                 OutputName -> "HC1"}},
    {H02,      { Description -> "Neutral up Higgsinos",
                 PDG -> {0},
                 PDG.IX -> {0},
                 FeynArtsNr -> 124,
                 LaTeX -> "\\tilde{H}^0_2",
                 OutputName -> "H02"}},
    {HC2,      { Description -> "Charged up Higgsinos",
                 PDG -> {0},
                 PDG.IX -> {0},
                 FeynArtsNr -> 125,
                 LaTeX -> "\\tilde{H}^+_2",
                 OutputName -> "HC2"}},
    {H0I1,     { Description -> "Neutral Inert-down-Higgsinos",
                 PDG -> {0, 0},
                 PDG.IX -> {0, 0},
                 Width -> 0,
                 FeynArtsNr -> 721,
                 LaTeX -> "\\tilde{h}^{0,Inert}_1",
                 OutputName -> "HNI1"}},
    {HCI1,     { Description -> "Charged Inert-down-Higgsinos",
                 PDG -> {0, 0},
                 PDG.IX -> {0, 0},
                 Width -> 0,
                 FeynArtsNr -> 722,
                 LaTeX -> "\\tilde{h}^{-,Inert}_1",
                 OutputName -> "HCI1"}},
    {H0I2,     { Description -> "Neutral  Inert-up-Higgsinos",
                 PDG -> {0, 0},
                 PDG.IX -> {0, 0},
                 Width -> 0,
                 FeynArtsNr -> 723,
                 LaTeX -> "\\tilde{h}^{0,Inert}_2",
                 OutputName -> "HNI2"}},
    {HCI2,     { Description -> "Charged Inert-up-Higgsinos",
                 PDG -> {0, 0},
                 PDG.IX -> {0, 0},
                 Width -> 0,
                 FeynArtsNr -> 724,
                 LaTeX -> "\\tilde{h}^{+}_2",
                 OutputName -> "HCI2"}},
    {Hp01,     { Description -> "Neutral Prime-Higgsinos",
                 PDG -> {0},
                 PDG.IX -> {0},
                 Width -> 0,
                 FeynArtsNr -> 950,
                 LaTeX -> "\\tilde{h}^{'0}_1",
                 OutputName -> "HNP1"}},
    {HpC1,     { Description -> "Charged Prime-Higgsinos",
                 PDG -> {0},
                 PDG.IX -> {0},
                 Width -> 0,
                 FeynArtsNr -> 951,
                 LaTeX -> "\\tilde{h}^{'-}_1",
                 OutputName -> "HCP1"}},
    {Hp02,     { Description -> "Neutral  Prime-Bar-Higgsinos",
                 PDG -> {0},
                 PDG.IX -> {0},
                 Width -> 0,
                 FeynArtsNr -> 952,
                 LaTeX -> "\\tilde{\\bar{h}}^{'0}_2",
                 OutputName -> "HNP2"}},
    {HpC2,     { Description -> "Charged Prime-Bar-Higgsinos",
                 PDG -> {0},
                 PDG.IX -> {0},
                 Width -> 0,
                 FeynArtsNr -> 953,
                 LaTeX -> "\\tilde{\\bar{h}}^{'+}_2",
                 OutputName -> "HCP2"}},
    {FSI1,     { Description -> "Dirac Left Inert Singlino",
                 PDG -> {0, 0, 0},
                 PDG.IX -> {0, 0, 0},
                 Width -> 0,
                 FeynArtsNr -> 812,
                 LaTeX -> "\\tilde{s}_L",
                 OutputName -> "FSI1"}},
    {FSI2,     { Description -> "Dirac Right Inert Singlino",
                 PDG -> {0, 0, 0},
                 PDG.IX -> {0, 0, 0},
                 Width -> 0,
                 FeynArtsNr -> 813,
                 LaTeX -> "\\tilde{s}_R",
                 OutputName -> "FSI2"}},
    {FS1,      { Description -> "Dirac Left Singlino",
                 PDG -> {0},
                 PDG.IX -> {0},
                 Width -> 0,
                 FeynArtsNr -> 806,
                 LaTeX -> "\\tilde{s}_L",
                 OutputName -> "FSL"}},
    {FS2,      { Description -> "Dirac Right Singlino",
                 PDG -> {0},
                 PDG.IX -> {0},
                 Width -> 0,
                 FeynArtsNr -> 807,
                 LaTeX -> "\\tilde{s}_R",
                 OutputName -> "FSR"}},
    {FSbar1,   { Description -> "Dirac Left SinglinoBar",
                 PDG -> {0},
                 PDG.IX -> {0},
                 Width -> 0,
                 FeynArtsNr -> 808,
                 LaTeX -> "\\tilde{sbar}_L",
                 OutputName -> "FSbarL"}},
    {FSbar2,   { Description -> "Dirac Right SinglinoBar",
                 PDG -> {0},
                 PDG.IX -> {0},
                 Width -> 0,
                 FeynArtsNr -> 809,
                 LaTeX -> "\\tilde{sbar}_R",
                 OutputName -> "FSbarR"}},
    {Fphi1,    { Description -> "Dirac Left Phinglino",
                 PDG -> {0},
                 PDG.IX -> {0},
                 Width -> 0,
                 FeynArtsNr -> 810,
                 LaTeX -> "\\tilde{phi}_L",
                 OutputName -> "FphiL"}},
    {Fphi2,    { Description -> "Dirac Right Phinglino",
                 PDG -> {0},
                 PDG.IX -> {0},
                 Width -> 0,
                 FeynArtsNr -> 811,
                 LaTeX -> "\\tilde{phi}_R",
                 OutputName -> "FphiR"}},
    {Fd1,      { Description -> "Dirac Left Down-Quark"}},
    {Fd2,      { Description -> "Dirac Right Down-Quark"}},
    {Fu1,      { Description -> "Dirac Left Up-Quark"}},
    {Fu2,      { Description -> "Dirac Right Up-Quark"}},
    {Fe1,      { Description -> "Dirac Left Electron"}},
    {Fe2,      { Description -> "Dirac Right Electron"}},
    {Fv,       { Description -> "Neutrinos" }},
    {FDx1,     { Description -> "Dirac Left Exotics",
                 LaTeX -> "Dx_1",
                 FeynArtsNr -> 660,
                 OutputName -> "Dx1"}},
    {FDx2,     { Description -> "Dirac Right Exotics",
                 LaTeX -> "Dx_2",
                 FeynArtsNr -> 661,
                 OutputName -> "Dx2"}}
};

ParticleDefinitions[TEMP] = {
    {VZ,       { Description -> "Z-Boson" }},
    {gZ,       { Description -> "Z-Boson Ghost" }}
};

ParticleDefinitions[EWSB] = {
    {Sd ,      { Description -> "Down-Squarks"}},
    {Su ,      { Description -> "Up-Squarks"}},
    {Se ,      { Description -> "Sleptons"}},
    {Sv ,      { Description -> "Sneutrinos"}},
    {SDX,      { Description -> "SExotics",
                 PDG -> {1000061, 2000061, 1000062, 2000062, 1000063, 2000063},
                 PDG.IX -> {-200890207, -200890208, -200890209, -200890210, -200890211, -200890212},
                 Mass -> Automatic,
                 ElectricCharge -> -1/3,
                 FeynArtsNr -> 666,
                 LaTeX -> "\\tilde{x}",
                 OutputName -> "SDX"}},
    {hh,       { Description -> "Higgs",
                 PDG -> {25, 35, 45, 65, 75},
                 PDG.IX -> {101000001,101000002,101000003, 101000004, 101000005}}},
    {Ah,       { Description -> "Pseudo-Scalar Higgs",
                 PDG -> {0, 0, 36, 46, 56},
                 PDG.IX -> {0, 0, 102000001, 102000002, 102000003}}},
    {Hpm,      { Description -> "Charged Higgs"}},
    {hhI ,     { Description -> "Neutral Inert-Higgs",
                 PDG -> {68, 69, 70, 71, 72, 73, 74},
                 PDG.IX -> {201000001, 201000002, 201000003, 201000004, 201000005, 201000006, 201000007},
                 FeynArtsNr -> 7,
                 ElectricCharge -> 0,
                 LaTeX -> "h_{I}",
                 OutputName -> "hhI"}},
    {AhI,      { Description -> "Inert Pseudo-Scalar Higgs",
                 PDG -> {9900068, 9900069, 9900070, 9900071, 9900072, 9900073, 9900074},
                 PDG.IX -> {201000008, 201000009, 201000010, 201000011, 201000012, 201000013, 201000014},
                 FeynArtsNr -> 6,
                 ElectricCharge -> 0,
                 LaTeX -> "A_{I}",
                 OutputName -> "AhI"}},
    {SHIPM,    { Description -> "Charged Inert-Higgs",
                 PDG -> {47, 57, 67, 77},
                 PDG.IX -> {-200000607, -200000608, -200000609, -200000610},
                 FeynArtsNr -> 8,
                 ElectricCharge -> -1 ,
                 LaTeX -> {"h^{-,Inert}","h^{+,Inert}"},
                 OutputName -> "HCI"}},
    {SHp0,     { Description -> "Neutral Prime-Higgs",
                 PDG -> {78, 79},
                 PDG.IX -> {200000011, 200000012},
                 FeynArtsNr -> 900,
                 ElectricCharge -> 0,
                 LaTeX -> "H^{'0}",
                 OutputName -> "HPN"}},
    {SHpp,     { Description -> "Charged Prime-Higgs",
                 PDG -> {58, 59},
                 PDG.IX -> {-200000611, -200000612},
                 ElectricCharge -> -1,
                 FeynArtsNr -> 901,
                 LaTeX -> {"H^{'-}","H^{'+}"},
                 OutputName -> "HPA"}},
    {VP,       { Description -> "Photon"}},
    {gP,       { Description -> "Photon Ghost"}},
    {VZ,       { Description -> "Z-Boson" }},
    {gZ,       { Description -> "Z-Boson Ghost" }},
    {VZp,      { Description -> "Z'-Boson"}},
    {gZp,      { Description -> "Z'-Ghost" }},
    {VWm,      { Description -> "W-Boson" }},
    {gWm,      { Description -> "Negative W-Boson Ghost"}},
    {gWmC,     { Description -> "Positive W-Boson Ghost" }},
    {VG,       { Description -> "Gluon" }},
    {gG,       { Description -> "Gluon Ghost" }},
    {Fd,       { Description -> "Down-Quarks"}},
    {Fu,       { Description -> "Up-Quarks"}},
    {Fe,       { Description -> "Leptons" }},
    {Fv,       { Description -> "Neutrinos" }},
    {FDX,      { Description->"Exotics",
                 PDG -> {61, 62, 63},
                 PDG.IX -> {-210890201, -210890202, -210890203},
                 Width -> 0,
                 Mass -> Automatic,
                 ElectricCharge -> -1/3,
                 FeynArtsNr -> 666,
                 LaTeX -> "x",
                 OutputName -> "FDX"}},
    {Glu,      { Description -> "Gluino" }},
    {Chi,      { Description -> "Neutralinos",
                 PDG -> {1000022, 1000023, 1000025, 1000035, 1000045, 1000065, 1000075, 1000032},
                 PDG.IX -> {211000001,211000002,211000003,211000004,211000005,211000006, 211000007, 211000008},
                 FeynArtsNr -> 11 }},
    {ChiI,     { Description -> "Inert Neutralinos",
                 PDG -> {1000068, 1000069, 1000070, 1000071, 1000072, 1000073, 1000074},
                 PDG.IX -> {211000009, 211000010, 211000011, 211000012, 211000013, 211000014, 211000015},
                 Mass -> Automatic,
                 FeynArtsNr -> 13,
                 ElectricCharge -> 0,
                 LaTeX -> "\\tilde{\\chi}^{0,Inert}",
                 OutputName -> "NI"}},
    {ChiP,     { Description -> "Prime Neutralinos",
                 PDG -> {1000078, 1000079},
                 PDG.IX -> {211000016, 211000017},
                 ElectricCharge -> 0,
                 FeynArtsNr -> 900,
                 LaTeX -> "\\tilde{\\chi}^{'0}",
                 OutputName -> "NP"}},
    {Cha,      { Description -> "Charginos",
                 FeynArtsNr -> 12}},
    {ChaI,     { Description -> "Inert Charginos",
                 PDG -> {-1000047, -1000057},
                 PDG.IX -> {-210000603, -210000604},
                 Mass -> Automatic,
                 ElectricCharge -> -1,
                 FeynArtsNr -> 14,
                 LaTeX -> {"\\tilde{\\chi}^{-,Inert}",
                           "\\tilde{\\chi}^{+,Inert}"},
                 OutputName -> "AI"}},
    {ChaP,     { Description -> "Prime Chargino",
                 PDG -> {-1000058},
                 PDG.IX -> {-210000605},
                 ElectricCharge -> -1,
                 Mass -> Automatic,
                 FeynArtsNr -> 901,
                 LaTeX -> {"\\tilde{\\chi}^{'-}",
                           "\\tilde{\\chi}^{'+}"},
                 OutputName -> "AP"}}
};

WeylFermionAndIndermediate = {
    {FHd0,     { Description -> "Neutral Down-Higgsino"}},
    {FHu0,     { Description -> "Neutral Up-Higgsino" }},
    {FHdm,     { Description -> "Charged Down-Higgsino"}},
    {FHup,     { Description -> "Charged Up-Higgsino"}},
    {FHpd0,    { Description -> "Neutral L4-ino",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\tilde{L}_{4}^{0}",
                 OutputName -> ""}},
    {FHpu0,    { Description -> "Neutral L4-Bar-ino",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\tilde{\\bar{L}}_{4}^{0}",
                 OutputName -> ""}},
    {FHpdm,    { Description -> "Charged L4-ino",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\tilde{L}_{4}^{-}",
                 OutputName -> ""}},
    {FHpup,    { Description -> "Charged L4-Bar-ino",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\tilde{\\bar{L}}_{4}^{+}",
                 OutputName -> ""}},
    {FH1I0,    { Description -> "Neutral Inert-down-Higgsino",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\tilde{h}^{0,Inert}_{d}",
                 OutputName -> ""}},
    {FH2I0,    { Description -> "Neutral Inert-up-Higgsino",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\tilde{h}^{0,Inert}_{u}",
                 OutputName -> ""}},
    {FH1Im,    { Description -> "Charged Inert-down-Higgsino",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\tilde{h}^{-,Inert}_{d}",
                 OutputName -> ""}},
    {FH2Ip,    { Description -> "Charged Inert-up-Higgsino",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\tilde{h}^{+,Inert}_{u}",
                 OutputName -> ""}},
    {FSIR,     { Description -> "Weyl-Inert-Singlino",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\tilde{S}^{Inert}",
                 OutputName -> ""}},
    {L0,       { Description -> "Neutralino Weyl-Spinor"}},
    {Lm,       { Description -> "Negative Chargino Weyl-Spinor"}},
    {Lp,       { Description -> "Positive Chargino Weyl-Spinor"}},
    {L0I,      { Description -> "Neutralino Inert-Weyl-Spinor",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\lambda_I^0",
                 OutputName -> ""}},
    {LmI,      { Description -> "Negative Chargino Inert-Weyl-Spinor",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\lambda_I^-",
                 OutputName -> ""}},
    {LpI,      { Description -> "Positive Chargino Inert-Weyl-Spinor",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\lambda_I^+",
                 OutputName -> ""}},
    {L0p,      { Description -> "Neutralino Prime-Weyl-Spinor",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\lambda^{\\prime 0}",
                 OutputName -> ""}},
    {Lmp,      { Description -> "Chargino Prime-Weyl-Spinor",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\lambda^{\\prime\\pm}",
                 OutputName -> ""}},
    {fG,       { Description ->"Gluino Weyl-Spinor"}},
    {fWB,      { Description ->"Wino Weyl-Spinor"}},
    {fW0,      { Description ->"Neutral Wino" }},
    {fWm,      { Description ->"Negative Wino"}},
    {fWp,      { Description ->"Positive Wino"}},
    {fB,       { Description ->"Bino Weyl-Spinor"}},
    {phid,     { Description -> "Scalar Down"}},
    {phiu,     { Description -> "Scalar Up"}},
    {phiS,     { Description -> "Scalar Singlet",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\phi_{S}",
                 OutputName -> ""}},
    {phiSbar,  { Description -> "Scalar SingletBar",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\phi_{\\overline{S}}",
                 OutputName -> ""}},
    {phiPhi,   { Description -> "Scalar Pure Singlet",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\phi_{\\phi}",
                 OutputName -> ""}},
    {phiH1I0,  { Description -> "Scalar Inert Down",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\tilde{\\phi}_{d}",
                 OutputName -> ""}},
    {phiH2I0,  { Description -> "Scalar Inert Up",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\tilde{\\phi}_{u}",
                 OutputName -> ""}},
    {phiSIR,   { Description -> "Scalar Inert Singlet",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\tilde{\\phi}_{S}",
                 OutputName -> ""}},
    {sigmad,   { Description -> "Pseudo Scalar Down"}},
    {sigmau,   { Description -> "Pseudo Scalar Up" }},
    {sigmaS,   { Description -> "Pseudo Scalar Singlet",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\sigma_{S}",
                 OutputName -> ""}},
    {sigmaSbar,{ Description -> "Pseudo Scalar SingletBar",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\sigma_{\\overline{S}}",
                 OutputName -> ""}},
    {sigmaPhi, { Description -> "Pseudo Scalar Pure Singlet",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\sigma_{\\phi}",
                 OutputName -> ""}},
    {sigmaH1I0,{ Description -> "Pseudo Scalar Inert Down",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\tilde{\\sigma}_{d}",
                 OutputName -> ""}},
    {sigmaH2I0,{ Description -> "Pseudo Scalar Inert Up",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\tilde{\\sigma}_{u}",
                 OutputName -> ""}},
    {sigmaSIR, { Description -> "Pseudo Scalar Inert Singlet",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\tilde{\\sigma}_{S}",
                 OutputName -> ""}},
    {SHd,      { Description -> "Down-Higgs"}},
    {SHu,      { Description -> "Up-Higgs"}},
    {Sl,       { Description -> "Left Slepton" }},
    {Sq,       { Description -> "Left Squark" }},
    {FeL,      { Description -> "Left Electron" }},
    {FeR,      { Description -> "Right Electron" }} ,
    {FdL,      { Description -> "Left Down-Quark" }},
    {FdR,      { Description -> "Right Down-Quark" }},
    {FuL,      { Description -> "Left Up-Quark" }},
    {FuR,      { Description -> "Right Up-Quark" }},
    {FEL,      { Description -> "Rotated Left Electron" }},
    {FER,      { Description -> "Rotated Right Electron" }} ,
    {FDL,      { Description -> "Rotated Left Up-Quark" }},
    {FDR,      { Description -> "Rotated Right Up-Quark" }},
    {FUL,      { Description -> "Rotated Left Down-Quark"}},
    {FUR,      { Description -> "Rotated Right Down-Quark" }},
    {FHd,      { Description -> "Down-Higgsino" }},
    {FHu,      { Description -> "Up-Higgsino" }},
    {Fl,       { Description -> "Left Leptons"}},
    {Fq,       { Description -> "Left Quarks"}},
    {FvL,      { Description -> "Left Neutrino"}},
    {FsR,      { Description -> "Weyl Spinor of Singlino"}},
    {FsbarR,   { Description -> "Weyl Spinor of Singlino-Bar",
                 LaTeX -> "\\tilde{\\overline{S}}",
                 OutputName -> "sb"}},
    {FphiR,    { Description -> "Weyl Spinor of Pure Singlino",
                 LaTeX -> "\\tilde{\\phi}",
                 OutputName -> "fphi"}},
    {FHp,      { Description -> "Weyl Spinor of L4",
                 LaTeX -> "\\tilde{L}_{4}",
                 OutputName -> "FHp"}},
    {FHpbar,   { Description -> "Weyl Spinor of L4-Bar",
                 LaTeX -> "\\tilde{\\bar{L}}_{4}",
                 OutputName -> "FHpbar"}},
    {fBp,      { Description -> "Bino prime Weyl-Spinor",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\lambda_{\\tilde{B}^\\prime}",
                 OutputName -> "" }},
    {FDxL,     { Description -> "Left Exotics",
                 LaTeX -> "D x_{L}",
                 OutputName -> "FDxL"}},
    {FDxbarR,  { Description -> "Right Exotics",
                 LaTeX -> "D x_{R}",
                 OutputName -> "FDxR"}},
    {FDXL,     { Description -> "Rotated Left Exotics",
                 LaTeX -> "X_L",
                 OutputName -> ""}},
    {FDXR,     { Description->"Rotated Right Exotics",
                 LaTeX -> "X_R",
                 OutputName -> ""}},
    {B,        { Description -> "B Superfield" }},
    {WB,       { Description -> "W Superfield" }},
    {G,        { Description -> "Gluon Superfield" }},
    {Bp,       { Description -> "U(1)' gauge Superfield",
                 LaTeX -> "\\hat{B}'"}},
    {q,        { Description -> "Left Quark Superfield" }},
    {l,        { Description -> "left Lepton Superfield" }},
    {Hd,       { Description -> "Down-Higgs Superfield"}},
    {Hu,       { Description -> "Up-Higgs Superfield" }},
    {d,        { Description -> "Right Down-Quark Superfield" }},
    {u,        { Description -> "Right Up-Quark Superfield" }},
    {e,        { Description -> "Right Electron Superfield" }},
    {s,        { Description -> "Singlet Superfield" }},
    {sbar,     { Description -> "SingletBar Superfield",
                 LaTeX -> "\\hat{\\overline{S}}"}},
    {H2I,      { Description -> "Inert-Up-Higgs Superfield",
                 LaTeX -> "\\hat{H}^{Inert}_{2}"}},
    {H1I,      { Description -> "Inert-Down-Higgs Superfield",
                 LaTeX -> "\\hat{H}^{Inert}_{1}"}},
    {SI,       { Description -> "Inert-Singlet Superfield",
                 LaTeX -> "\\hat{S}^{Inert}"}},
    {Dx,       { Description -> "Left Exotics Superfield",
                 LaTeX -> "\\hat{Dx}"}},
    {Dxbar,    { Description -> "Right Exotics Superfield",
                 LaTeX -> "\\hat{\\bar{Dx}}"}},
    {Hp,       { Description -> "L4 Superfield",
                 LaTeX -> "\\hat{L}_{4}"}},
    {Hpbar,    { Description -> "L4-Bar Superfield",
                 LaTeX -> "\\hat{\\bar{L}}_{4}"}},
    {phi,      { Description -> "Pure singlet Superfield",
                 LaTeX -> "\\hat{\\phi}"}},
    {SH1I,     { Description-> "Scalar Inert Down Higgs",
                 LaTeX -> "H^{Inert}_1"}},
    {SH2I,     { Description-> "Scalar Inert Up Higgs",
                 LaTeX -> "H^{Inert}_2"}},
    {SHp,      { Description-> "Scalar L4",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "L_4",
                 OutputName -> ""}},
    {SHpbar,   { Description-> "Scalar L4-Bar",
                 PDG -> 0,
                 Width -> 0,
                 Mass -> Automatic,
                 LaTeX -> "\\bar{L}_4",
                 OutputName -> ""}}
};
