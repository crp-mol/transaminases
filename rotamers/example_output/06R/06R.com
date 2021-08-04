%Mem=8192MB
%NProcShared=4
%chk=65307.chk
#P AM1 Opt=(MaxCycle=900,InitialHarmonic=50000)

first AM1 optimization

-1  1
N          -0.78997         8.28096        22.20400
H          -1.60997         8.02596        22.75100
C          -0.13597         7.28196        21.62200
C          -0.62897         5.82696        21.79200
H          -1.47297         5.79496        22.32700
H           0.05503         5.26996        22.26300
H          -0.81097         5.40396        20.90400
C           1.04703         7.58396        20.89800
O           1.71103         6.56696        20.32300
C           1.49703         8.94296        20.73600
C           2.76103         9.34696        19.99700
H           3.11476        10.27572        19.88627
C           0.74103         9.97196        21.38900
C          -0.41397         9.59996        22.09900
H          -0.97197        10.30596        22.53500
C           1.15603        11.43696        21.33300
H           0.58003        11.99996        21.92500
H           2.11503        11.55296        21.59300
O           0.97303        11.84096        19.94500
P           1.84103        13.16096        19.49900
O           1.80103        14.19496        20.54000
O           1.22703        13.55496        18.16000
O           3.36603        12.67796        19.35000
C           8.39767        12.79979        20.42438
C           7.30123        12.34587        19.48689
O           7.31659        10.90741        19.44628
C           6.39636        10.28331        18.65040
C           5.13754         9.94285        19.37736
C           4.47197         8.65915        18.85051
C           5.51095         7.53479        18.88778
N           3.24821         8.35166        19.57662
H           2.87290         7.43335        19.70250
O           6.71969        10.05473        17.47756
H           8.39036        13.91623        20.45975
H           9.39561        12.44985        20.06209
H           8.22595        12.39959        21.45470
H           6.29154        12.65297        19.85385
H           7.45690        12.75151        18.45631
H           5.40130         9.80052        20.46045
H           4.43629        10.82248        19.31580
H           3.88043         8.77073        17.90082
H           6.43247         7.82095        18.31524
H           5.09492         6.60311        18.43540
H           5.80585         7.31368        19.94443

--link1--
%chk=65307.chk
#P HF/6-31G* Opt=(MaxCycle=900,InitialHarmonic=50000) SCF=(Tight,MaxCycle=900) Geom=AllCheck Guess=Read Pop=MK IOp(6/33=2,6/41=10,6/42=17)
