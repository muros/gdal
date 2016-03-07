Nastaviš VC okolje na 64 bitni verziji builda:

c:\Dev\GIS\gdal-2.0.1>"c:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\vcvarsall.bat" amd64

Poženeš build v gdal-2.0.1 folderju:

nmake /f makefile.vc MSVC_VER=1400 GDAL_HOME="C:\gdal2\bld"

Isti nmake lahko poženeš najprej tudi le v svojem driver folderju, da vidiš da dela prevajanje.
Sam celoten build traja minuto dve.
In potem še dve do tri install

Nato še installiraš na lokacijo, ki je nastavljena v build konfiguraciji,
spremenljivka BINDIR in DATADIR:

nmake /f makefile.vc install

Namesto install lahko tudi devinstall, èe rabiš knjižnice za razvijanje.

Install je sedaj nastavljen na pot:
C:\gdal2\bld

kjer so tudi bin -i in data. Na to pot nastavi path, da ti delajo gdal tooli z pravo verzijo.

Primer poganjanja gdalwarp za sestavljanje kosov DMV5 za en Tehnièni List TN:

c:\gdal2\bld\bin\gdalwarp C26\*.xyz radovljica.tiff

Poganjaj v cygwin, ker zna expandad imena filetov, ki jih podaš z * (asterisk)
gdalwarp, èe hoèeš dodajat najprej naredi s prvim in zadnjim quadrantom da dobiš geographic
extent.

Kako v qgisu narediš mape nevarnosti plazov:
- Najprej z gdalwarp narediš tiff z elevation data
- Ta tiff uvoziš kot raster layer
- Na tem tifu narediš raster operacijo za slope (z factor ostane 1.0)
- Na slope narediš raster operacijo za relief - poimenuješ terain.
- Terain layerju daš transparentnost 65%
- Lahko narediš Raster - Extraction - Contour in narediš še plastnice (na 10 m)
- Uporabiš terain layer in pod njim OpenStreetMap layer.
- Lahko dodaš še Google Satellite layer za lažjo orientacijo
