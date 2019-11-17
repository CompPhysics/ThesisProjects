1

 tbibdc: two-boson-exchange - iterative box diagrams
 ---------------------------------------------------
                  with nd and dd intermediate states
                  ----------------------------------
0input-parameter-set:
 --------------------
0
0tbibd parameters                                                    
0nd          1  0  0
0dd          1  0  0
0is,m,wnd    0  0  1
0is,m,wdd    0  0  1
0ide         1  4  0  0  0  0
0idde        1  4  0  0  0  0
0nj         24 24 24 24 24 24 24  0  0  0  0  0  0  0  0  0  0  0  0  0
0c           800.0000
0wn          938.9190
0ws         1232.0000
0wpi         138.0400
0icmplx      0
0beta          0.0000
0itj         4
1

 obnn:  one-boson-exchange nn-nn interaction (numer. integra.)
 -------------------------------------------------------------
0input-parameter-set:
 --------------------
0MODEL I (BN I), R. Machleidt, Adv. Nucl. Phys. 19, 189 (1989), App. B 
0factor typ  1
0num.int.    4 48
0nucl. mass  938.9190
0
 jp  name      g**2      f/g       mass    iso-spiniprop/spe
         cut typ     c u t - o f f   p a r a m e t e r s
 -------------------------------------------------------------
00-  pion     13.6000    0.0000  138.0400    1.0       0.0
 cut       2.0        0.0   2.0000     1720.0000    0.0000
0end mesons
 -------------------------------------------------------------
 -------------------------------------------------------------
1

 obnd:  one-boson-exchange nn-nd interaction (numer. integra.)
 -------------------------------------------------------------
0input-parameter-set:
 --------------------
0     obnd parameters                                                  
0factor typ  3
0num. int.   4 48
0nucl. mass  938.9190
0n-star mas 1232.0000
0
 jp  name      g**2      f/g       mass    iso-spiniprop/spe
         cut typ     c u t - o f f   p a r a m e t e r s
 -------------------------------------------------------------
00-  pion      2.1832    0.0000  138.0388    1.0       0.0
 cut nn    2.0        0.0   1.0000     1720.0000    0.0000
 cut nd    2.0        0.0   1.0000      850.0000    0.0000
01-  rho       4.0988    6.1042  769.9000    1.0       0.0
 cut nn    2.0        0.0   1.0000     1310.0000    0.0000
 cut nd    2.0        0.0   2.0000     1310.0000    0.0000
0end mesons
 -------------------------------------------------------------
 -------------------------------------------------------------
1

 obdd:  one-boson-exchange nn-dd interaction (numer. integra.)
 -------------------------------------------------------------
0input-parameter-set:
 --------------------
0     obdd parameters                                                  
0factor typ  3
0num. int.   4 48
0nucl. mass  938.9190
0n-star mas 1232.0000
0
 jp  name      g**2      f/g       mass    iso-spiniprop/spe
         cut typ     c u t - o f f   p a r a m e t e r s
 -------------------------------------------------------------
00-  pion      0.3500    0.0000  138.0388    1.0       0.0
 cut       2.0        0.0   2.0000      850.0000    0.0000
01-  rho      20.0000    0.0000  769.9000    1.0       0.0
 cut       2.0        0.0   4.0000     1310.0000    0.0000
0end mesons
 -------------------------------------------------------------
 -------------------------------------------------------------
