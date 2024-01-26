# **PolyGEN**

Design automation of polycistronic tRNA-based genes containing custom RNAs for assembly in Type IIs restriction enzyme-driven Golden Gate experiments. The backbone of PolyGEN is based on
[iBioCAD](https://ibiocad.igb.illinois.edu/) by HamediRad et al. (2019), for which the code is openly available [here](https://github.com/scottweis1/iBioCAD).

The code takes as input an array of custom RNAs and will compute the finished PTG together with the necessary oligomers to produce all parts from a plasmid containing a gRNA-tRNA template. Currently, the produced PTGs can include sgRNAs, pegRNAs, crRNAs and other small RNAs. By default, PolyGEN will use the following parameters, which can be varied

- primer melting temperature between 55 and 65 °C if possible
- Digestion by BsaI
- 'tgcc' and 'gttt' as restriction overlaps with the plasmid
- no additional restriction sites in the plasmid

To calculate the primer melting temperatures, PolyGEN uses the same method and parameters as Benchling: [SantaLucia (1998)](https://www.pnas.org/content/95/4/1460).

_____________

**Setup Linux/Mac**

First, install docker by running apt-get install docker

In the terminal, run through the following pipeline

- Clone this repository via `git clone https://git.hhu.de/urquizag/polygen`
- Navigate into the cloned repo `cd polygen`
- execute `docker-compose up` (requires docker desktop to be active)

In the browser, open localhost:5000

**Setup Windows**

Activate Windows Subsystem for Linux (WSL2) by 

- open Control Panel
- open Turn Windows features on or off
- check features Virtual Machine Platform and Windows Subsystem for Linux
- confirm with OK
- restart the pc
- download the Update Setup from [here](https://wslstorestorage.blob.core.windows.net/wslblob/wsl_update_x64.msi)
- Execute the wsl_update_x64.msi file
- open a command prompt
- run wsl --set-default-version 2

In the BIOS, enable Visualization tools. Next, install [Docker](https://docs.docker.com/desktop/windows/install/) and clone this git repository. With active docker, open a command prompt and navigate into the cloned repository using dir \<location\>. Execute docker-compose up.

In the browser, open localhost:5000
  
____________
  
**Common problems**
  
If the install fails due to issues with docker-snap:
  
- sudo rm -rf /etc/docker
- sudo snap refresh
