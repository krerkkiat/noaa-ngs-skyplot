# National Geodetic Survey's Skyplot

---
**NOTE**

Original source code can be obtained from [C++ code for creating a Skyplot of GPS satellites](https://geodesy.noaa.gov/gps-toolbox/skyplot.htm)

---

## How to use skyplot

1. Create a working directory (ie, `c:\work`) and update `skyplot.inp` with information related to the site. An example `skyplot.inp` is shown below as

```plain
1995  9  2  18 30 0
1995  9  2  20 30 0 
brdc2450.95n
GRAZ 47.0671 15.4935 538.3
10
```

Note 1:  Make sure the orbit file (either broadcast or SP3)
covers the time span requested in lines 1 and 2 of skyplot.inp.

Note 2:  It may be necessary to supply the entire path
for the orbit file (eg. `c:\work\brdc2450.95n`)

2. Execute skyplot.exe at the DOS prompt. For example:

```console
c:\work> c:\bin\skyplot.exe
```

3. Examine `skyplot.log` to ensure that `skyplot.inp` was read correctly.  `skyplot.log` should report 

```plain
Normal Termination
```

4. Assuming GMT (Generic Mapping Tools) exists and is properly loaded on the computer, issue the command `skyplot.bat`.  For example:

```console
c:\work\skyplot.bat
```

Note:  If ERROR messages related to `psxy`, `pstext`,
or `psvelo` appear, then it is likely that GMT is
improperly loaded on the computer.  See the GMT home 
page ([http://gmt.soest.hawaii.edu/](http://gmt.soest.hawaii.edu/)) for assistance.

5. Use ghostview or your favorite postscript viewer to
   examine and print the file 'skyplot.ps'.

## Modifications

- Apply `clang-format`.
- Linux focus (for now).
- Dirty hack to support RINEX file version 2.10.
- Update [GMT](https://docs.generic-mapping-tools.org/latest/index.html)'s commands.
