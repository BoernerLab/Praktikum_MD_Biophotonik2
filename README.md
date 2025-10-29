# Praktikum Biophotonik II – Molekulardynamik (MD)

Dieses Repository enthält die Materialien zum **Praktikum Biophotonik II** an der Hochschule Mittweida.  
Im Rahmen des Versuchs wird eine **Molekulardynamik-Simulation eines RNA-Hairpins** mit **GROMACS** durchgeführt und anschließend ausgewertet.

---

## Ziel des Praktikums

- Durchführung einer vollständigen MD-Simulation eines RNA-Hairpins mit GROMACS  
- Verständnis der einzelnen Phasen:  
  1. **Energieminimierung**  
  2. **Äquilibrierung (NVT/NpT)**  
  3. **Produktionslauf**
- Analyse der Trajektorie und Bewertung struktureller Veränderungen mithilfe von **PyMOL** und **Python (MDAnalysis)**

---

## Inhalt des Ordners `MD_Praktikum`

| Datei / Ordner | Beschreibung |
|----------------|---------------|
| `RNA_hairpin.pdb` | Strukturdatei des RNA-Hairpins |
| `amber14sb_OL15.ff/` | Kraftfeldparameter für RNA und Proteine |
| `mdp/` | Parameterdateien für Energieminimierung, NVT-, NpT- und Produktionsläufe |
| `Praktikum_MD.ipynb` | Jupyter Notebook zur Analyse der Simulationsergebnisse |

Autoren: Felix Erichson & Richard Börner