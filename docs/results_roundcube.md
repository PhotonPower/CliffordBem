# Ergebnisse: Würfel mit abgerundeten Kanten im Resonanzfenster

Würfel [−1, 1]³, alle Kanten (und damit Ecken) mit Radius ρ gerundet (Gmsh/OpenCASCADE, `fillet`), ωa = 0,5,
ebene Welle, `scatter_mesh` mit Sauter-Schwab/halbanalytischem Nahfeld, H-Matrix-Toleranz 10⁻³, GMRES 10⁻⁵.

Netz (`tools/make_geometries.py roundcube h out.msh --radius ρ --curv c`): auf den Rundungen Elementgröße
2πρ/c, auf den flachen Seiten mit dem Abstand zur Rundung linear bis h wachsend (Gmsh-Felder Distance/Threshold,
Übergang über max(h, 3ρ)). Verfeinerungsfolge (h, c) = (0,4; 8), (0,3; 12), (0,22; 16). Ohne diese Abstufung
(flache Seiten gleich grob neben der Rundung) konvergiert die Folge für kleines ρ nicht erkennbar (ρ = 0,1:
16,03 → 16,44 → 16,91), weil das an der Rundung konzentrierte Feld auf den angrenzenden flachen Stücken nicht
aufgelöst ist.

## Konvergenz bei festem Radius

| | c = 8 | c = 12 | c = 16 | Differenzen |
|---|---:|---:|---:|---:|
| ρ = 0,2, κ = −2,8 + 0,3i | 14,2237 (1 504) | 13,3077 (3 294) | 12,8649 (5 282) | −0,92 / −0,44 |
| ρ = 0,3, κ = −2,8 + 0,3i | 19,6033 (888) | 18,3948 (1 930) | 17,8708 (3 282) | −1,21 / −0,52 |
| ρ = 0,4, κ = −2,8 + 0,3i | 25,1225 (552) | 24,7659 (1 208) | 24,3995 (2 048) | −0,36 / −0,37 |
| ρ = 0,3, Gold −11 + 1,2i | 10,2095 | 10,3878 | 10,4130 | +0,18 / +0,025 |

(σ_ext, in Klammern Dreieckszahl.) Mit abgerundeten Kanten verhält sich das Fenster wie ein gewöhnliches Problem:
Die Folgen sind monoton, die Differenzen halbieren sich (ρ = 0,2 und 0,3; beobachtete Ordnung etwa 1 in der
Elementgröße an der Rundung). Zum Vergleich der scharfe Würfel bei κ = −3 + 0,3i (v0.5): 18,24 → 18,69 → 18,59 →
18,33, nicht monoton. Außerhalb des Fensters (Gold) konvergiert es deutlich schneller. Bei ρ = 0,4 ist der
asymptotische Bereich mit dieser Folge noch nicht erreicht.

## Abhängigkeit vom Radius

Die Extinktion hängt im Fenster stark vom Rundungsradius ab: 24,4 (ρ = 0,4), 17,9 (0,3), 12,9 (0,2) auf dem jeweils
feinsten Netz; grob extrapoliert (Ordnung 1) etwa 23, 16 und 11,5. Eine Halbierung des Radius von 0,4 auf 0,2
halbiert die Extinktion ungefähr. Der Rundungsradius ist damit im Fenster ein entscheidender physikalischer
Parameter; einen Grenzwert für ρ → 0 gibt es bei kleinem Verlust nicht zu erwarten.

## Bewertung

- Die Abrundung macht das Resonanzfenster rechenbar: Die Lösung konvergiert mit dem Netz, was an scharfen Kanten
  nicht gelang. Der Preis ist die Auflösung des Radius; die Kosten wachsen etwa mit 1/ρ (ρ = 0,2: 5 282 Dreiecke,
  110 s; ρ = 0,1 mit Abstufung wäre ein Vielfaches).
- Die Konvergenz ist mit Ordnung etwa 1 langsamer als bei glatten, gut aufgelösten Körpern (Kugel: Ordnung 2);
  die Rundungen sind mit 2 bis 4 Elementen über den Viertelbogen noch grob. Feinere Folgen brauchen mehr Rechenzeit
  als hier verfügbar.
- Für AP 4 bedeutet das: Vorhersagen für Nanowürfel im Fenster sind nur zusammen mit dem Rundungsradius sinnvoll,
  der aus Elektronenmikroskopie bekannt sein muss. Das ist physikalisch realistisch und experimentell überprüfbar.
