Nastavi� VC okolje na 64 bitni verziji builda:

c:\Dev\GIS\gdal-2.0.1>"c:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\vcvarsall.bat" amd64

Po�ene� build v gdal-2.0.1 folderju:

nmake /f makefile.vc MSVC_VER=1400 GDAL_HOME="C:\gdal2\bld"

Isti nmake lahko po�ene� najprej tudi le v svojem driver folderju, da vidi� da dela prevajanje.
Sam celoten build traja minuto dve.
In potem �e dve do tri install

Nato �e installira� na lokacijo, ki je nastavljena v build konfiguraciji,
spremenljivka BINDIR in DATADIR:

nmake /f makefile.vc install

Namesto install lahko tudi devinstall, �e rabi� knji�nice za razvijanje.

Install je sedaj nastavljen na pot:
C:\gdal2\bld

kjer so tudi bin -i in data. Na to pot nastavi path, da ti delajo gdal tooli z pravo verzijo.

Primer poganjanja gdalwarp za sestavljanje kosov DMV5 za en Tehni�ni List TN:

c:\gdal2\bld\bin\gdalwarp C26\*.xyz radovljica.tiff

Poganjaj v cygwin, ker zna expandad imena filetov, ki jih poda� z * (asterisk)
gdalwarp, �e ho�e� dodajat najprej naredi s prvim in zadnjim quadrantom da dobi� geographic
extent.

Kako v qgisu naredi� mape nevarnosti plazov:
- Najprej z gdalwarp naredi� tiff z elevation data
- Ta tiff uvozi� kot raster layer
- Na tem tifu naredi� raster operacijo za slope (z factor ostane 1.0)
- Na slope naredi� raster operacijo za relief - poimenuje� terain.
- Terain layerju da� transparentnost 65%
- Lahko naredi� Raster - Extraction - Contour in naredi� �e plastnice (na 10 m)
- Uporabi� terain layer in pod njim OpenStreetMap layer.
- Lahko doda� �e Google Satellite layer za la�jo orientacijo
