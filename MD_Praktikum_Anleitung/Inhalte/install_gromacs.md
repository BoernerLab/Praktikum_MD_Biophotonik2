# Installationsanleitung: GROMACS

## 1. Vorbereitung: Windows Subsystem for Linux (WSL) und Ubuntu installieren

1. **WSL aktivieren**  
   - Öffnen Sie die Windows-Kommandozeile (`cmd`) oder PowerShell und geben Sie ein:  
     ```bash
     wsl --install
     ```  
   - Bei Windows 11 ist WSL bereits vorinstalliert. Falls nicht, folgen Sie den Anweisungen zur Installation.  
   - Nach der Installation muss der Rechner neu gestartet werden.

2. **Windows-Feature aktivieren**  
   - Suchen Sie in der Windows-Suchleiste nach **„Windows-Features aktivieren oder deaktivieren“**.  
   - Aktivieren Sie die Option **„Windows-Subsystem für Linux“**.  
   - Starten Sie den Rechner erneut.

3. **Ubuntu installieren**  
   - Öffnen Sie den **Microsoft Store** und installieren Sie **Ubuntu** (empfohlen: Ubuntu 22.04 LTS).  
   - Starten Sie Ubuntu und vergeben Sie einen **Benutzernamen** und ein **Passwort**.  
   - Merken Sie sich beides, da Sie es bei Installationen benötigen.

---

## 2. GROMACS installieren

1. **Ubuntu aktualisieren**  
   ```bash
   sudo apt update && sudo apt upgrade -y
   ```

2. **Erforderliche Pakete installieren**  
   ```bash
   sudo apt install build-essential cmake wget -y
   ```

3. **GROMACS herunterladen**  
   ```bash
   wget https://ftp.gromacs.org/gromacs/gromacs-2025.3.tar.gz
   ```

   Mit  
   ```bash
   ls
   ```  
   können Sie prüfen, ob die Datei im aktuellen Verzeichnis liegt.

4. **Archiv entpacken**  
   ```bash
   tar xfz gromacs-2025.3.tar.gz
   cd gromacs-2025.3
   ```

5. **Build-Ordner anlegen und betreten**  
   ```bash
   mkdir build
   cd build
   ```

6. **Konfiguration mit CMake**  
   ```bash
   cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
   ```

7. **Kompilieren**  
   ```bash
   make -j$(nproc)
   ```

8. **Optional: Installation prüfen**  
   ```bash
   make check
   ```

9. **Installation durchführen**  
   ```bash
   sudo make install
   ```

10. **GROMACS-Befehle aktivieren**  
    ```bash
    source /usr/local/gromacs/bin/GMXRC
    ```

---

## 3. Test der Installation

Prüfen Sie die Installation mit:  
```bash
gmx --version
```
Wenn die Versionsnummer von **GROMACS 2025.3** erscheint, ist die Installation erfolgreich abgeschlossen.

---

## 4. Troubleshooting (Häufige Fehler und Lösungen)

### Problem 1: `wsl: command not found`
- Ursache: WSL ist nicht installiert oder nicht aktiviert.  
- Lösung:  
  - Öffnen Sie PowerShell als Administrator.  
  - Führen Sie erneut aus:  
    ```bash
    wsl --install
    ```
  - Starten Sie den Rechner neu.

---

### Problem 2: `cmake: command not found`
- Ursache: CMake wurde nicht installiert.  
- Lösung:  
  ```bash
  sudo apt install cmake -y
  ```

---

### Problem 3: `make: command not found`
- Ursache: Entwicklungswerkzeuge fehlen.  
- Lösung:  
  ```bash
  sudo apt install build-essential -y
  ```

---

### Problem 4: `wget: command not found`
- Ursache: Das Download-Programm wget fehlt.  
- Lösung:  
  ```bash
  sudo apt install wget -y
  ```

---

### Problem 5: `gmx: command not found` nach Installation
- Ursache: GROMACS-Umgebungsvariablen wurden nicht geladen.  
- Lösung:  
  - Laden Sie GROMACS manuell:  
    ```bash
    source /usr/local/gromacs/bin/GMXRC
    ```
  - Um dies dauerhaft zu aktivieren, fügen Sie die Zeile in Ihre `~/.bashrc` ein:  
    ```bash
    echo "source /usr/local/gromacs/bin/GMXRC" >> ~/.bashrc
    ```

---

### Problem 6: Speicher oder Berechtigungsfehler bei `make`
- Ursache: Kompilierung benötigt mehr Speicher oder Adminrechte.  
- Lösung:  
  - Kompilieren mit nur einem Thread:  
    ```bash
    make -j1
    ```
  - Stellen Sie sicher, dass Sie genügend Speicherplatz haben (ca. 3–5 GB).

---

### Problem 7: `tar: command not found`
- Ursache: Das Archivierungsprogramm fehlt.  
- Lösung:  
  ```bash
  sudo apt install tar -y
  ```

---

✅ Mit diesen Hinweisen sollten die häufigsten Installationsprobleme behoben werden können.
