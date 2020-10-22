    
 ==============================================
 ==========*********************===============
 ==========***    Bhldem     ***===============
 ==========*********************===============
 ==============================================
    


 ===========================================================================
 =                                                                         =
 =               BBB   B  B  B    B   B  B    B  B                         =
 =               B  B  B  B  B    B   B  BB  BB  B                         =
 =               BBB   BBBB  B    B   B  B BB B  B                         =
 =               B  B  B  B  B    B   B  B    B  B                         =
 =               BBB   B  B  BBBB  BBB   B    B  B                         =
 =                                                                         =
 =               *********************************                         =
 =               *   BHLUMI version 4.04         *                         =
 =               *   June           1991  (2.01) *                         =
 =               *   Sept           1992  (2.02) *                         =
 =               *   January        1995  (4.00) *                         =
 =               *   Febuary        1995  (4.01) *                         =
 =               *   May            1995  (4.02) *                         =
 =               *   July           1995  (4.02a)*                         =
 =               *   June           1996  (4.03) *                         =
 =               *   September      1996  (4.04) *                         =
 =               *         AUTHORS               *                         =
 =               * S. Jadach,      W. Placzek,   *                         =
 =               * E. Richter-Was, B.F.L. Ward,  *                         =
 =               *        and Z. Was             *                         =
 =               *********************************                         =
 ===========================================================================



 ===========================================================================
 =                This program is based on papers                          =
 =                --------------------------------                         =
 =                Phys. Lett. B353 (1995) 362                              =
 =                Phys. Lett. B353 (1995) 349                              =
 =                Comp. Phys. Comm. 70 (1992) 305                          =
 =                Phys. Lett. B268 (1991) 253                              =
 =                Phys. Rev.  D40  (1989) 3582.                            =
 =                Phys. Lett. B260 (1991) 438,                             =
 =                Phys. Lett. B253 (1991) 469,                             =
 =                Nucl. Phys. B228 (1983) 537.                             =
 ===========================================================================



 ===========================================================================
 =               *********************************                         =
 =                 BHLUM4: INPUT PARAMETRES                                =
 =               *********************************                         =
 =             3001                  OPTIONS   switch        KeyOpt        =
 =                1                  rand. numb. switch      KeyRnd        =
 =                0                  weighting switch        KeyWgt        =
 =                0                  photon removal  sw      KeyRem        =
 =             1022                  RADIATION switch        KeyRad        =
 =                2                  vac_pol   switch        KeyPia        =
 =                2                  QED mat. elm. type      KeyMod        =
 =                0                  Test switch, def=0      KeyUpd        =
 =                1                  Z contribution          KeyZet        =
 =      92.30000000                  CMS energy  [GeV]       CMSENE        =
 =    8519.2900                      CMSENE^2  [GeV^2]       SVAR          =
 =    .60110696                      trasf_min [GeV^2]       TRMIN         =
 =    28.626770                      trasf_max [GeV^2]       TRMAX         =
 =    .70558340E-04                  xi_min=TRMIN/SVAR       XIMIN         =
 =    .33602295E-02                  xi_max=TRMAX/SVAR       XIMAX         =
 =      16.80000000                  theta_min  [mrad]       THMIN         =
 =     116.00000000                  theta_max  [mrad]       THMAX         =
 =        .96256910                  theta_min   [deg]       THMIN         =
 =       6.64631042                  theta_max   [deg]       THMAX         =
 =    .10000000E-03                  eps_CM infr. cut        EPSCM         =
 =    .10000000E-05                  delta  infr. cut        DEL           =
 =       -.01480534                  RePi(transf_min)        REPI1         =
 =        .00006655                  error                   dREPI1        =
 =       -.02999273                  RePi(transf_max)        REPI2         =
 =        .00039879                  error                   dREPI2        =
 =      91.18700000                  Z-mass GeV              AMAZ          =
 =       2.49000000                  Z-width GeV             GAMMZ         =
 =        .23190000                  weak mixing angle       SINW2         =
 ===========================================================================

 =====================DUMPS====================
  P2    -.9609181940494   -.6160149114072  46.1358664389564  46.1499838621984
  Q2     .9545403944187    .6104560314024 -45.5574498144363  45.5715375841869
 PHO     .0063560979334    .0055449682609   -.5773812076829    .5774428159041
 SUM    -.0000217016973   -.0000139117439    .0010354168372  92.2989642622895
 =====================DUMPS====================
  P2   -2.4418242325506    .2723107119134  44.6086703983944  44.6762815561043
  Q2    1.6178305345983   -.2364783177992  -8.9404046993193   9.0886816470983
 PHO     .9030588479010   -.0447881482344 -37.0950338729553  37.1060515185421
 PHO    -.0790651499487    .0089557541202   1.4267681738787   1.4289852782521
 SUM     .0000000000000    .0000000000000   -.0000000000016  92.2999999999968
 =====================DUMPS====================
  P2    -.0129875340879   -.8317405695647  46.0867405163538  46.0942470676553
  Q2     .0130324175520    .8314305741683 -46.1388420951031  46.1463345974699
 PHO    -.0000566197983    .0002414399873    .0557593456378    .0557598971038
 SUM    -.0000117363342   -.0000685554090    .0036577668886  92.2963415622291
 =====================DUMPS====================
  P2    -.9083292882332   -.9862115935446  46.1305163415870  46.1499968971370
  Q2     .9081306013081    .9859928842136 -46.1056921324685  46.1251745842978
 PHO    -.0000027766703    .0000006876038   -.0146014691137    .0146014693939
 PHO     .0002014635954    .0002180217271   -.0102227400011    .0102270491436
 SUM     .0000000000000    .0000000000000    .0000000000036  92.2999999999722
 =====================DUMPS====================
  P2     .4746473996194    .6810662454498  25.2185012547013  25.2321609692720
  Q2    -.4758230037785   -.6814076872925 -46.1380646712303  46.1455494694243
 PHO     .0011860451397    .0003503077060  20.9208442693057  20.9208443058582
 SUM     .0000104409806    .0000088658634    .0012808527767  92.2985547445544
1
   9001           bhlum4, weight distribution                                                    
           nent            sum           bmin           bmax
         290817     .30039E+06     .00000E+00     .65321E+05
           undf           ovef           avex
     .00000E+00     .10000E+01     .10329E+01
 -1.0000    .000000D+00 0                                                                 I
  -.8000    .000000D+00 0                                                                 I
  -.6000    .000000D+00 0                                                                 I
  -.4000    .200000D+01 0                                                                 I
  -.2000    .119000D+03 0                                                                 I
   .0000    .599550D+05 0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     I
   .2000    .120290D+05 0XXXXXXXXXXXX                                                     I
   .4000    .108490D+05 0XXXXXXXXXX                                                       I
   .6000    .138010D+05 0XXXXXXXXXXXXX                                                    I
   .8000    .170260D+05 0XXXXXXXXXXXXXXXXX                                                I
  1.0000    .187540D+05 0XXXXXXXXXXXXXXXXXX                                               I
  1.2000    .290290D+05 0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                    I
  1.4000    .639040D+05 0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX I
  1.6000    .653210D+05 0XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  1.8000    .140000D+02 0                                                                 I
  2.0000    .900000D+01 0                                                                 I
  2.2000    .100000D+01 0                                                                 I
  2.4000    .000000D+00 0                                                                 I
  2.6000    .000000D+00 0                                                                 I
  2.8000    .100000D+01 0                                                                 I
  3.0000    .100000D+01 0                                                                 I
  3.2000    .000000D+00 0                                                                 I
  3.4000    .100000D+01 0                                                                 I
  3.6000    .000000D+00 0                                                                 I
  3.8000    .000000D+00 0                                                                 I
  4.0000    .000000D+00 0                                                                 I
  4.2000    .000000D+00 0                                                                 I
  4.4000    .000000D+00 0                                                                 I
  4.6000    .000000D+00 0                                                                 I
  4.8000    .000000D+00 0                                                                 I
  5.0000    .000000D+00 0                                                                 I
  5.2000    .000000D+00 0                                                                 I
  5.4000    .000000D+00 0                                                                 I
  5.6000    .000000D+00 0                                                                 I
  5.8000    .000000D+00 0                                                                 I
  6.0000    .000000D+00 0                                                                 I
  6.2000    .000000D+00 0                                                                 I
  6.4000    .000000D+00 0                                                                 I
  6.6000    .000000D+00 0                                                                 I
  6.8000    .000000D+00 0                                                                 I


 ===========================================================================
 =               *********************************                         =
 =                 BHLUM4:        WINDOW A                                 =
 =               *********************************                         =
 =           100000                  Accepted total          NEVGEN     A1 =
 =           290817                  Raw prior reject.       IEVENT     A2 =
 =    438.34090      +-  .49338970   Xsec M.C. [nb]          XSECMC     A3 =
 =        .00112558                  relat. error            ERELMC     A4 =
 =       1.03290582  +-  .00112558   weight  M.C.            AWT        A5 =
 =              121                  WT<0                    NEVNEG     A6 =
 =                3                  WT>WTMAX                NEVOVE     A7 =
 =       3.00000000                  Maximum WT              WWMX       A8 =
 ===========================================================================



 ===========================================================================
 =               *********************************                         =
 =                 BHLUM4:        WINDOW B                                 =
 =               *********************************                         =
 =        .60977300  +-  .00103037  WT1*WT2*T/TP*T/TQ                   B1 =
 =        .99995125  +-  .00001303  WT3 from KINO4                      B2 =
 =       1.73758801  +-  .00021375  YFS formfac              WT         B4 =
 =       1.03290582  +-  .00112558  TOTAL                               B5 =
 =        .00002771  +-  .00002582  xsec/xtot: WT>WTMAX      WT         B6 =
 =       -.00000313  +- -.00000126  xsec/xtot: WT<0          WT         B7 =
 ===========================================================================



 ===========================================================================
 =               *********************************                         =
 =                           WINDOW C                                      =
 =               Built-in average control weights.                         =
 =               Should equal one +- statist. err.                         =
 =               *********************************                         =
 =       1.00153607  +-  .00129753  <WCTA1>                             C1 =
 =        .99979999  +-  .00130102  <WCTA2>                             C2 =
 =       1.00180566  +-  .00204721  <WCTA1*WCTA2>                       C3 =
 =       1.00055123  +-  .00048902  <WCTB1>                             C4 =
 =       1.00043689  +-  .00048958  <WCTB2>                             C5 =
 =       1.00087433  +-  .00070501  <WCTB1*WCTB2>                       C6 =
 ===========================================================================

1
   1000          Theta distribution (radians)                                                    
           nent            sum           bmin           bmax
         200000     .13180E+04     .14750E+03     .18950E+04
           undf           ovef           avex
     .00000E+00     .00000E+00     .65902E-02
   .0240    .189500D+04 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXI
   .0246    .188700D+04 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXI
   .0251    .176100D+04 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     I
   .0257    .163700D+04 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX         I
   .0263    .165650D+04 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX         I
   .0268    .148750D+04 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX               I
   .0274    .139750D+04 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                  I
   .0280    .132600D+04 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                     I
   .0285    .125300D+04 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                        I
   .0291    .116300D+04 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                           I
   .0297    .110350D+04 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                             I
   .0302    .103400D+04 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                I
   .0308    .104100D+04 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                I
   .0314    .100750D+04 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                 I
   .0319    .927500D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                    I
   .0325    .902000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                     I
   .0331    .865500D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXX                                      I
   .0336    .807000D+03 XXXXXXXXXXXXXXXXXXXXXXXXX                                         I
   .0342    .775000D+03 XXXXXXXXXXXXXXXXXXXXXXXX                                          I
   .0348    .781500D+03 XXXXXXXXXXXXXXXXXXXXXXXX                                          I
   .0353    .689000D+03 XXXXXXXXXXXXXXXXXXXXX                                             I
   .0359    .636000D+03 XXXXXXXXXXXXXXXXXXX                                               I
   .0365    .644000D+03 XXXXXXXXXXXXXXXXXXX                                               I
   .0370    .612500D+03 XXXXXXXXXXXXXXXXXX                                                I
   .0376    .605500D+03 XXXXXXXXXXXXXXXXXX                                                I
   .0382    .577500D+03 XXXXXXXXXXXXXXXXX                                                 I
   .0387    .510500D+03 XXXXXXXXXXXXXX                                                    I
   .0393    .512000D+03 XXXXXXXXXXXXXX                                                    I
   .0399    .480500D+03 XXXXXXXXXXXXX                                                     I
   .0404    .487500D+03 XXXXXXXXXXXXX                                                     I
   .0410    .455000D+03 XXXXXXXXXXXX                                                      I
   .0416    .424500D+03 XXXXXXXXXXX                                                       I
   .0421    .433500D+03 XXXXXXXXXXX                                                       I
   .0427    .407000D+03 XXXXXXXXXX                                                        I
   .0433    .377000D+03 XXXXXXXXX                                                         I
   .0438    .361500D+03 XXXXXXXXX                                                         I
   .0444    .333000D+03 XXXXXXXX                                                          I
   .0450    .323000D+03 XXXXXXX                                                           I
   .0455    .333000D+03 XXXXXXXX                                                          I
   .0461    .296500D+03 XXXXXX                                                            I
   .0467    .287000D+03 XXXXXX                                                            I
   .0472    .322000D+03 XXXXXXX                                                           I
   .0478    .278500D+03 XXXXX                                                             I
   .0484    .269000D+03 XXXXX                                                             I
   .0489    .282000D+03 XXXXXX                                                            I
   .0495    .250000D+03 XXXX                                                              I
   .0501    .246000D+03 XXXX                                                              I
   .0506    .237500D+03 XXXX                                                              I
   .0512    .219500D+03 XXX                                                               I
   .0518    .219500D+03 XXX                                                               I
   .0523    .208500D+03 XXX                                                               I
   .0529    .217000D+03 XXX                                                               I
   .0535    .195500D+03 XX                                                                I
   .0540    .199500D+03 XX                                                                I
   .0546    .161000D+03 X                                                                 I
   .0552    .153000D+03 X                                                                 I
   .0557    .185500D+03 XX                                                                I
   .0563    .165000D+03 X                                                                 I
   .0569    .148500D+03 X                                                                 I
   .0574    .147500D+03 X                                                                 I
1
   1100          Energy distr: x=log10(1-s1/s)                                                   
           nent            sum           bmin           bmax
         100000    -.58604E+07     .20300E+03     .96200E+03
           undf           ovef           avex
     .58610E+04     .00000E+00    -.58604E+02
 -6.0000    .216000D+03 XX                                                                I
 -5.9000    .203000D+03 X                                                                 I
 -5.8000    .238000D+03 XXXX                                                              I
 -5.7000    .243000D+03 XXXX                                                              I
 -5.6000    .239000D+03 XXXX                                                              I
 -5.5000    .268000D+03 XXXXXX                                                            I
 -5.4000    .301000D+03 XXXXXXXXX                                                         I
 -5.3000    .307000D+03 XXXXXXXXXX                                                        I
 -5.2000    .313000D+03 XXXXXXXXXX                                                        I
 -5.1000    .301000D+03 XXXXXXXXX                                                         I
 -5.0000    .314000D+03 XXXXXXXXXX                                                        I
 -4.9000    .341000D+03 XXXXXXXXXXXXX                                                     I
 -4.8000    .369000D+03 XXXXXXXXXXXXXXX                                                   I
 -4.7000    .354000D+03 XXXXXXXXXXXXXX                                                    I
 -4.6000    .364000D+03 XXXXXXXXXXXXXXX                                                   I
 -4.5000    .434000D+03 XXXXXXXXXXXXXXXXXXXXX                                             I
 -4.4000    .410000D+03 XXXXXXXXXXXXXXXXXXX                                               I
 -4.3000    .484000D+03 XXXXXXXXXXXXXXXXXXXXXXXXX                                         I
 -4.2000    .470000D+03 XXXXXXXXXXXXXXXXXXXXXXXX                                          I
 -4.1000    .491000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXX                                        I
 -4.0000    .474000D+03 XXXXXXXXXXXXXXXXXXXXXXXX                                          I
 -3.9000    .442000D+03 XXXXXXXXXXXXXXXXXXXXX                                             I
 -3.8000    .482000D+03 XXXXXXXXXXXXXXXXXXXXXXXXX                                         I
 -3.7000    .465000D+03 XXXXXXXXXXXXXXXXXXXXXXX                                           I
 -3.6000    .484000D+03 XXXXXXXXXXXXXXXXXXXXXXXXX                                         I
 -3.5000    .518000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXX                                      I
 -3.4000    .520000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXX                                      I
 -3.3000    .553000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                   I
 -3.2000    .479000D+03 XXXXXXXXXXXXXXXXXXXXXXXXX                                         I
 -3.1000    .564000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                  I
 -3.0000    .599000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                               I
 -2.9000    .565000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                                  I
 -2.8000    .637000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                            I
 -2.7000    .615000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                              I
 -2.6000    .640000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                           I
 -2.5000    .654000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                          I
 -2.4000    .680000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                        I
 -2.3000    .688000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                       I
 -2.2000    .696000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                       I
 -2.1000    .758000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                 I
 -2.0000    .749000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                  I
 -1.9000    .769000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                I
 -1.8000    .844000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX          I
 -1.7000    .839000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX          I
 -1.6000    .912000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    I
 -1.5000    .875000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX       I
 -1.4000    .884000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX      I
 -1.3000    .950000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX I
 -1.2000    .897000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX     I
 -1.1000    .962000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 -1.0000    .906000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    I
  -.9000    .912000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    I
  -.8000    .908000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    I
  -.7000    .880000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX       I
  -.6000    .855000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX         I
  -.5000    .778000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX               I
  -.4000    .651000D+03 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX                           I
  -.3000    .477000D+03 XXXXXXXXXXXXXXXXXXXXXXXX                                          I
  -.2000    .329000D+03 XXXXXXXXXXX                                                       I
  -.1000    .252000D+03 XXXXX                                                             I
 |||||||||||||||||||||||||||||||||||||||||||||||||||
 Xsec_BARE1 =    169.19520371 Nanob.
 error      =       .67481969 Nanob.
 Xsec_CALO2 =    136.21881786 Nanob.
 error      =       .64151939 Nanob.
 |||||||||||||||||||||||||||||||||||||||||||||||||||
