# 5 Durchführung

## Vorbereitung

Laden Sie das bereitgestellte Repository
herunter (<a href="https://github.com/BoernerLab/Praktikum_MD_Biophotonik2" target="_blank">hier</a>) und entpacken Sie
es.
Öffnen Sie anschließend die Ubuntu-Terminalumgebung und kopieren Sie den entpackten Ordner in Ihr Linux-Dateisystem.
Dieses finden Sie im Dateiexplorer unter „Linux“. Kopieren Sie den Ordner in das Verzeichnis
/home/"username" (wobei username Ihrem bei der Installation gewählten Benutzernamen entspricht).

Öffnen Sie nun das Ubuntu-Terminal und navigieren Sie in den Ordner:

```bash
cd /home/username/MD_Praktikum/
```

In diesem Ordner befinden sich folgende Dateien und Unterordner:

- RNA_hairpin.pdb
- amber14sb_OL15.ff
- mdp/
- Praktikum_MD.ipynb

Die erste Datei enthält die Struktur des RNA-Hairpins
im <a href="https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html" target="_blank">PDB-Format</a>.
Das Kraftfeld (amber14sb_OL15.ff) definiert alle Parameter zur Berechnung der potentiellen Energie des Systems, also
Bindungen, Winkel,
Torsionen und nicht-gebundene Wechselwirkungen.
Die .mdp-Dateien enthalten die Simulationsparameter (Integrationsschritt, Temperatur- und Druckkopplung, Cutoff-Werte
usw.) für die folgenden Phasen der MD-Simulation. Das JupyterNotebook enthält Funktionen für die Auswertung der
Ergebnisse.

Starten Sie anschließend GROMACS, indem Sie die Umgebungsvariablen laden:

```bash
source /usr/local/gromacs/bin/GMXRC 
```

Damit steht der Befehl gmx systemweit zur Verfügung.

```{Note}
Führen Sie diesen Befehl immer aus, wenn Sie das Ubuntu-Terminal neu starten
```

## 5.1 pdb2gmx

Mit dem Befehl pdb2gmx wird aus der PDB-Strukturdatei eine vollständige Moleküldefinition erstellt.
Das Programm generiert dabei drei Ausgabedateien:

1. Topologiedatei (RNA_hairpin.top) – enthält alle atomaren Parameter (Atomtypen, Ladungen, Bindungen, Winkel usw.)
   gemäß
   Kraftfelddefinition.

2. Positionsrestriktionsdatei (RNA_hairpin.itp) – listet Atome, deren Positionen während bestimmter Simulationsphasen (
   z. B.
   Äquilibrierung) fixiert oder harmonisch eingeschränkt werden.

3. Strukturdatei (RNA_hairpin.gro) – GROMACS-kompatibles Koordinatenformat, konsistent mit der Kraftfeldbeschreibung.

Die Positionsrestriktionen sind insbesondere in der Äquilibrierung nützlich, damit sich das Lösungsmittel an die RNA
anpasst und nicht umgekehrt.

Erstellen Sie zunächst ein Unterverzeichnis und führen Sie den folgenden Befehl aus:

```bash
mkdir em
gmx pdb2gmx -f RNA_hairpin.pdb -o em/RNA_hairpin.gro -p RNA_hairpin.top -i em/RNA_hairpin.itp
```

- f: Eingabe (Struktur des Moleküls)
- o: Ausgabe GROMACS kompatible Strukturdatei
- p: Ausgabe Topologie Datei
- i: Ausgabe positionsrestriktion Datei

Nach dem Start von pdb2gmx werden Sie aufgefordert, ein Kraftfeld auszuwählen.

Für dieses Praktikum wird das Amber14SB-OL15-Kraftfeld verwendet, das speziell für RNA- und Protein-Simulationen
optimiert ist.
Wählen Sie die entsprechende Option (`1`) und bestätigen Sie mit `Enter`.

Im nächsten Schritt wählen Sie das Wassermodell.
Da RNA experimentell in wässriger Lösung untersucht wird, wird auch in der Simulation ein realistisches Wassermodell
benötigt.. Für dieses Praktikum nutzen wir das
Wassermodell 1 (TIP3P). Geben Sie dafür `1` ein drücken Sie `Enter`.

## 5.2 Wasserbox definieren und füllen

Im nächsten Schritt wird eine simulationsfähige Box definiert, die anschließend mit Wassermolekülen gefüllt wird.
Die Box legt das Volumen fest, in dem sich das RNA-Molekül und das Lösungsmittel bewegen dürfen. Da
MD-Simulationen unter periodischen Randbedingungen durchgeführt werden, wiederholt sich die Box in alle
Raumrichtungen, dadurch entsteht der Eindruck einer unendlichen Lösung.

#### Boxdefinition mit `editconf`

Zunächst muss die Form und Größe der Box festgelegt werden.
Während eine kubische Box die einfachste Variante ist, wird in der Praxis häufig eine rhombische Dodekaederbox
verwendet.
Diese Form spart Volumen, da sie das Molekül enger umschließt und somit weniger Wassermoleküle benötigt, während der
minimale Abstand zu den Boxgrenzen erhalten bleibt.

Die Box wird mit dem Befehl editconf erstellt:

```bash
gmx editconf -f em/RNA_hairpin.gro -o em/RNA_hairpin.gro -c -d 1.0 -bt dodecahedron

```

Parameter:

- f: Eingabedatei (Struktur des Moleküls)
- o: Ausgabedatei mit definierter Box
- c: Zentriert das Molekül im Zentrum der Box
- d 1.0: legt den minimalen Abstand zwischen Molekül und Boxrand auf 1,0 nm fest
- bt dodecahedron: legt den Boxtyp fest (rhombischer Dodekaeder)

```{admonition} Aufgabe
Notieren Sie sich das Boxvolumen und die Boxabmessungen, die nach Ausführung des Befehls im Terminal ausgegeben werden.
```

#### Wasser zum System hinzufügen

Nachdem die Box definiert wurde, kann sie nun mit Wasser gefüllt werden.
GROMACS verwendet dafür voräquilibrierte Wasserkoordinaten (Standarddatei spc216.gro), die durch Kopieren und
Zuschneiden in das Boxvolumen eingesetzt werden.
Überlappende Wassermoleküle werden automatisch entfernt.

Der Befehl lautet:

```bash
gmx solvate -cp em/RNA_hairpin.gro -cs spc216.gro -o em/RNA_hairpin.gro -p RNA_hairpin.top
```

Parameter:

- cp: Eingabedatei mit der definierten Box (aus editconf)
- cs: Koordinatendatei des Lösungsmittels (Standard: spc216.gro)
- o: Ausgabedatei mit solvatisierter Struktur
- p: aktualisiert die Topologie (topol.top) um die Anzahl der hinzugefügten Wassermoleküle

Hinweis:
Obwohl die Datei spc216.gro SPC-Wasser enthält, kann sie auch für TIP3P-Simulationen verwendet werden.
Die Datei liefert lediglich die Koordinaten für das Auffüllen der Box, die tatsächlichen physikalischen Parameter werden
durch das im Kraftfeld definierte Wassermodell (hier TIP3P) bestimmt.

```{admonition} Aufgabe
Nach diesem Prozess enthält die Ausgabe das RNA-Molekül umgeben von Wassermolekülen. Schauen Sie sich die Box an, indem Sie die erstellte `.gro` Datei in PyMOL laden.
```

## 5.3 Hinzufügen von Ionen

In diesem Schritt wird das System elektrisch neutralisiert und Ionen hinzugefügt.
Das geschieht mit dem GROMACS-Modul genion, welches ausgewählte Wassermoleküle durch Ionen ersetzt.

Damit genion korrekt arbeiten kann, benötigt es eine Run-Input-Datei (.tpr), die alle Informationen zum aktuellen
System
enthält. Diese Datei wird mit dem Tool grompp erzeugt (GROMACS Preprocessor).

grompp kombiniert die Eingabedateien .mdp, .gro und .top zu einer einzigen Binärdatei (.tpr), die das System
vollständig beschreibt.
Diese Datei kann anschließend von allen GROMACS-Modulen, z.B. genion oder mdrun, eingelesen werden.

Führen Sie den folgenden Befehl aus:

```bash
gmx grompp -f ./mdp/em.mdp -c em/RNA_hairpin.gro -p RNA_hairpin.top -o em/RNA_hairpin.tpr -po em/RNA_hairpin.mdp -maxwarn 2
```

Parameter:

- f: Eingabedatei mit den Simulationsparametern (hier: ions.mdp)
- c: aktuelle Koordinatendatei des solvatisierten Systems
- p: Topologiedatei, die Moleküldefinitionen und Kraftfeldparameter enthält
- o: Name der erzeugten Run-Input-Datei (ions.tpr)

Die Datei ions.mdp befindet sich im Ordner mdp/ und enthält die grundlegenden Kontrollparameter für diese
Vorbereitungsschritte.
Obwohl sie für eine „Simulation“ erstellt ist, wird in diesem Fall keine eigentliche MD-Rechnung ausgeführt, grompp
nutzt sie nur, um die vollständige Systembeschreibung auf atomarer Ebene zusammenzustellen.

#### Ionen hinzufügen mit `genion`

Mit der nun vorhandenen run-input-Datei kann das Ionisieren des Systems durchgeführt werden.
genion ersetzt dabei einzelne Wassermoleküle durch Ionen, die in der Topologie ergänzt werden.
Ionen werden benötigt, um die Nettoladung des Systems auszugleichen und damit die elektrostatische Berechnung (z. B. mit
dem Particle Mesh Ewald-Verfahren) physikalisch korrekt bleibt.
Bei Bedarf kann auch eine definierte Salzkonzentration (z.B. 0,1 M KCl) hinzugefügt werden.

Führen Sie den folgenden Befehl aus:

```bash
gmx genion -s em/RNA_hairpin.tpr -o em/RNA_hairpin.gro -p RNA_hairpin.top -pname K -nname Cl -neutral
```

Parameter:

- s: Eingabedatei (ions.tpr), erzeugt durch grompp
- o: Ausgabedatei mit den neuen Ionenkonfigurationen
- p: Topologiedatei, die automatisch um die Anzahl der Ionen ergänzt wird
- pname / nname: Symbolische Bezeichnungen für die positiven bzw. negativen Ionen
- neutral: weist GROMACS an, nur so viele Ionen hinzuzufügen, dass das Gesamtsystem elektrisch neutral wird

Während der Ausführung von genion werden Sie aufgefordert, eine Atomgruppe auszuwählen, deren Moleküle durch Ionen
ersetzt werden sollen.
Wählen Sie hier die Gruppe `SOL` aus und bestätigen Sie mit `Enter`.

```{admonition} Aufgabe
Berechnen Sie die Ionenkonzentration von Kalium in Ihrem System.
```

## 5.4 Energieminimierung

Nach dem Aufbau und der Neutralisierung des Systems kann die Molekulardynamik-Simulation noch nicht direkt gestartet
werden.  
Zunächst müssen sterische Kollisionen oder ungünstige Geometrien im System beseitigt werden, die während der
vorangegangenen Schritte entstanden sein können.  
Dies geschieht durch eine Energie-Minimierung (EM), bei der das System in einen energetisch
günstigen Startzustand gebracht wird.

---

##### Vorbereitung mit `grompp`

Wie in den vorherigen Schritten wird zunächst eine Run-Input-Datei (.tpr) erzeugt, die alle Informationen zur
Minimierung enthält.  
Dazu werden die Parameterdatei, die aktuelle Struktur und die Topologie kombiniert:

```bash
gmx grompp -f ./mdp/em.mdp -c em/RNA_hairpin.gro -p RNA_hairpin.top -o em/RNA_hairpin.tpr
```

Parameter:

- f: Eingabedatei mit Minimierungsparametern (minim.mdp)
- c: Startkoordinaten (ionisiertes System)
- p: Topologiedatei
- o: Ausgabedatei für die Run-Input-Datei (em.tpr)

Die Datei em.mdp enthält Einstellungen für die Energie-Minimierung.
In der Datei gibt es die Einstellung `integrator = steep`. Diese Einstellung minimiert die Energie des Moleküls durch
Änderung der Positionen der Atome. Dabei wird das Verfahren des steilsten Abstiegs verwendet.
Die Minimierung wird so lange fortgesetzt, bis der maximale Kraftwert unter einem definierten Schwellenwert (z. B. 1000
kJ mol⁻¹ nm⁻¹) liegt oder die maximale Schrittzahl erreicht ist.

##### Durchführung mit mdrun

Nach erfolgreicher Vorbereitung kann die Energie-Minimierung gestartet werden:

```bash
gmx mdrun -v -s em/RNA_hairpin.tpr -c em/RNA_hairpin.gro -o em/RNA_hairpin.trr -e em/RNA_hairpin.edr -g em/RNA_hairpin.log -po em/RNA_hairpin.mdp
```

Parameter:

- v: „verbose“: zeigt den Fortschritt der Minimierung in der Konsole an
- s: Eingabedatei (.tpr), erzeugt durch grompp
- c: Ausgabedatei mit der minimierten Struktur
- o: Trajektorie-Datei
- e: Binärdatei mit Energiewerten
- g: Logdatei mit Protokoll des Minimierungslaufs

Während der Minimierung werden keine dynamischen Bewegungen simuliert, sondern die atomaren Positionen entlang des
Gradienten des Potentialfelds angepasst, bis die potentielle Energie minimiert ist.

Nach erfolgreicher Minimierung ist das System nun spannungsfrei und bereit für die Äquilibrierungsphase (NVT/NPT).

## 5.5 Temperatur Äquilibrierung

In der NVT-Äquilibrierung wird das System langsam auf die Zieltemperatur gebracht, während das Volumen konstant
bleibt. Positionsrestriktionen halten die RNA weitgehend fixiert, sodass sich das Lösungsmittel an die Moleküloberfläche
anpassen kann. Dieser Schritt dient dazu, die Temperatur stabil zu etablieren, bevor der Druck angeglichen wird.

Zunächst muss GROMACS alle Eingabedateien (Parameter, Struktur und Topologie) zu einer run-input-Datei zusammenfassen:

```bash
mkdir nvt
gmx grompp -f ./mdp/nvt.mdp -c em/RNA_hairpin.gro -r em/RNA_hairpin.gro -p RNA_hairpin.top -o nvt/RNA_hairpin.tpr -po nvt/RNA_hairpin.mdp -maxwarn 2

```

- f: Parameterdatei im .mdp-Format mit den Einstellungen für die NVT-Simulation.
- c: Eingabestruktur aus der Energie-Minimierung.
- r: Referenzstruktur für Positionsrestriktionen.
- p: Topologiedatei mit Molekül- und Kraftfeldinformationen.
- o: Ausgabedatei (.tpr), die von mdrun eingelesen wird.
- po: schreibt die verwendeten Parameter in eine neue .mdp-Datei.
- maxwarn 2: erlaubt bis zu zwei Warnungen, ohne den Vorgang abzubrechen.

Führen Sie nun die NVT-Äquilibrierung mit folgendem Befehl aus:

```bash
gmx mdrun -v -s nvt/RNA_hairpin.tpr -c nvt/RNA_hairpin.gro -x nvt/RNA_hairpin.xtc -cpo nvt/RNA_hairpin.cpt -e nvt/RNA_hairpin.edr -g nvt/RNA_hairpin.log
```

- v: („verbose“) zeigt den Fortschritt der Simulation direkt im Terminal an.
- s: Eingabedatei im .tpr-Format, die alle Simulationsparameter und das System beschreibt (erstellt mit grompp).
- c: Ausgabedatei für die Endstruktur nach Abschluss der Simulation (.gro).
- x: Ausgabedatei der Trajektorie (.xtc), die die Atompositionen über die Zeit enthält.
- cpo: Checkpoint-Datei (.cpt), in der der aktuelle Zustand der Simulation gespeichert wird. Sie ermöglicht es, den Lauf
  zu einem späteren Zeitpunkt fortzusetzen.
- e: Energiedatei (.edr), die während der Simulation berechnete physikalische Größen wie Temperatur, Druck und Energie
  enthält.
- g: Logdatei (.log) mit einem detaillierten Protokoll des Simulationsverlaufs, inklusive Start- und Endbedingungen,
  Step-Informationen und Warnungen.

## 5.6 Druck Equillibrierung

Nachdem im vorherigen Schritt die Temperatur erfolgreich stabilisiert wurde, muss nun auch der ruck des Systems an
die gewünschten Bedingungen angepasst werden.  
Dieser Schritt erfolgt unter einem NPT-Ensemble, bei dem Teilchenzahl (N), Druck (P) und Temperatur (T)
konstant gehalten werden.

Während der NPT-Äquilibrierung wird das Volumen der Simulationsbox flexibel angepasst, bis die Dichte des Systems (
typischerweise etwa 1000 kg m⁻³ für Wasser) den Zielwert erreicht.  
Die RNA bleibt weiterhin durch Positionsrestriktionen leicht fixiert, sodass sich vor allem Lösungsmittel und Ionen an
die Moleküloberfläche anpassen.

---

Die Erstellung der Run-Input-Datei für die NPT-Phase erfolgt analog zur NVT-Äquilibrierung, jedoch mit geänderten
Parametern aus der Datei `npt.mdp`.  
In dieser Datei sind die Druckkopplung und weitere relevante Parameter definiert.  
Lemkul (2024) empfiehlt den Parrinello-Rahman-Barostat, da er bei biologischen Systemen eine realistische
Druckschwankung ermöglicht.

```bash
mkdir npt
gmx grompp -f ./mdp/npt.mdp -c nvt/RNA_hairpin.gro -r nvt/RNA_hairpin.gro -p RNA_hairpin.top -o npt/RNA_hairpin.tpr -po npt/RNA_hairpin.mdp -maxwarn 2
```

Starten Sie die Druck-Äquilibrierung mit:

```bash
gmx mdrun -v -s npt/RNA_hairpin.tpr -c npt/RNA_hairpin.gro -x npt/RNA_hairpin.xtc -cpo npt/RNA_hairpin.cpt -e npt/RNA_hairpin.edr -g npt/RNA_hairpin.log
```

Während dieses Schritts passt sich das Volumen des Systems an, um den gewünschten Druck (meist 1 bar) zu erreichen.  
Temperatur und Druck werden durch die in der `.mdp`-Datei definierten Thermostat- und Barostat-Kopplungen kontrolliert.

## 5.7 MD Simulaiton

Analog zu den vorherigen Äquilibrierungsschritten wird nun die eigentliche MD-Simulation (Produktionslauf) vorbereitet
und gestartet.  
Hierbei wird keine Positionsrestriktion mehr angewendet, alle Atome können sich frei bewegen, und das System entwickelt
sich entsprechend der in der Kraftfelddefinition vorgegebenen Wechselwirkungen.

#### Anpassung der Simulationszeit

Öffnen Sie nun die Datei mdp/md0.mdp in einem Texteditor.
Legen Sie den Parameter für die Simulationszeit (nsteps) selbst fest und bestimmen Sie, wie lange Ihre Simulation laufen soll.
Beachten Sie, dass das die Simulationslänge von den Zeitschritten abhängt, die auf 2fs eingestellt sind.
Starten Sie anschließend grompp und führen Sie die Simulation aus:

```bash
mkdir md0
gmx grompp -f ./mdp/md0.mdp -c npt/RNA_hairpin.gro -r npt/RNA_hairpin.gro -p RNA_hairpin.top -o md0/RNA_hairpin.tpr -po md0/RNA_hairpin.mdp -maxwarn 2
```

```bash
 gmx mdrun -v -s md0/RNA_hairpin.tpr -c md0/RNA_hairpin.gro -x md0/RNA_hairpin.xtc -cpo md0/RNA_hairpin.cpt -e md0/RNA_hairpin.edr -g md0/RNA_hairpin.log 
```

Beobachten Sie, wie lange die Berechnung auf Ihrem Rechner dauern würde, um ein Gefühl für die benötigte Rechenzeit zu
bekommen. Brechen Sie dann die Simulation ab.

```{admonition} Aufgabe
Im Anschluss stellen Sie die Simulationszeit auf eine Mikrosekunde (1 µs) ein und wiederholen den grompp-Schritt.
Laden Sie die dabei erzeugten Dateien anschließend auf Opal hoch.
Die Simulation wird dort auf einem dafür vorgesehenen Rechner ausgeführt.
Die Ergebnisse erhalten Sie nach Abschluss der Simulation zur weiteren Analyse und Auswertung.

```

