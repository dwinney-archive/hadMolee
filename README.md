# hadMolee

## **Had**ronic **Mol**ecules in **e+e-** collisions

Amplitude analysis toolkit for exploring lines-shaped of vector charmonium-like mesons in electron-positron annihilation data. 

## Installation

Installation comes in two parts: the core library and header-only implementations of physics amplitudes. After compiling the core library, amplitudes and analysis is done using ROOT's cling interpreter through such as those contained in the `/scripts` directory. 

The core library is built by cloning the directory and running cmake: 
```
mkdir build && cd build
cmake ..
cmake --build . --target install
```
This will install the linkable library in the `lib/` directory (e.g. `/lib/libHADMOLEE.so` on Linux) and an executable `/bin/hadMolee` will be installed. The executable is simply a shortcut to load the library and physics header-files into a ROOT interactive session and run a specified script with cling. 