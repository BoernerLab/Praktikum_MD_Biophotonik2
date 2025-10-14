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

4. **Archiv entpacken mit tar**  
   ```bash
   tar xfz gromacs-2025.3.tar.gz
   cd gromacs-2025.3
   ```
   Mit dem Befehl `cd` wechseln Sie in das entpackte Verzeichnis. 

5. **Build-Ordner anlegen und betreten**  
   ```bash
   mkdir build
   cd build
   ```

6. **Konfiguration der Installtion mit CMake**  
   ```bash
   cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON
   ```

7. **Kompilieren (Dieser Schritt kann einige Zeit in Anspruch nehmen)**  
   ```bash
   make -j$(nproc)
   ```

8. **Optional: Installation prüfen (Dieser Schritt kann einige Zeit in Anspruch nehmen)**  
   ```bash
   make check
   ```

9. **Installation durchführen**  
   ```bash
   sudo make install
   ```

10. **GROMACS-Befehle aktivieren**

    Dieser Befehl muss immer bei Neustart der Ubuntu Umgebung unter Windows durchgeführt werden.
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
