/******************************************************************************
 * $Id:$
 *
 * Project:  DMV Reader
 * Purpose:  All code for Slovenian DEM Reader
 * Author:   Uros Mesaric Kunst, upumesar@gmail.com
 *
 ******************************************************************************
 * Copyright (c) 2016, Uros Mesaric Kunst <upumesar@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ****************************************************************************/

#include "gdal_pam.h"
#include "cpl_conv.h"

CPL_CVSID("$Id:$");

CPL_C_START
void GDALRegister_DMV(void);
CPL_C_END

struct xyz
{
	double x;
	double y;
	double z;
};

struct dmvcell
{
	double originX;
	double originY;
	int widthY;
	int heightX;
};

typedef struct xyz XYZ;
typedef struct dmvcell DmvCell;

static int readLine(char line[], FILE* fp);
static void parseXyzLine(char *line, struct xyz* data);
static DmvCell getCellDimension(const char * fileName);
static XYZ getCellSize(int subNum);
static XYZ getCellOrign(char letter, int num, int subNum);

/************************************************************************/
/*                            DMVGetField()                             */
/************************************************************************/

static int DMVGetField(char *pszField, int nWidth)
{
	char szWork[32];

	CPLAssert(nWidth < (int) sizeof(szWork));

	strncpy(szWork, pszField, nWidth);
	szWork[nWidth] = '\0';

	return atoi(szWork);
}

/************************************************************************/
/* ==================================================================== */
/*                          JDEMDataset                                 */
/* ==================================================================== */
/************************************************************************/

class DMVRasterBand;

class DMVDataset: public GDALPamDataset

{
	friend class DMVRasterBand;

	VSILFILE *fp;
	GByte abyHeader[1012];

	int nRasterXOrigin;
	int nRasterYOrigin;

	double minX;
	double maxX;
	double minY;
	double maxY;

public:
	DMVDataset();
	~DMVDataset();

	static GDALDataset *Open(GDALOpenInfo *);
	static int Identify(GDALOpenInfo *);

	CPLErr GetGeoTransform(double * padfTransform);
	const char *GetProjectionRef();

	int GetRasterXOrigin();
	int GetRasterYOrigin();

	double GetMinX();
	double GetMaxX();
	double GetMinY();
	double GetMaxY();
};

/************************************************************************/
/* ==================================================================== */
/*                            DMVRasterBand                             */
/* ==================================================================== */
/************************************************************************/

class DMVRasterBand: public GDALPamRasterBand
{
	friend class DMVDataset;

	int nRecordSize;
	int bBufferAllocFailed;
	double **band;

public:

	DMVRasterBand(DMVDataset *, int);
	~DMVRasterBand();
	double calcPixAvg(int, int, int, int);

	virtual CPLErr IReadBlock(int, int, void *);
};

/*
 * Calculate average of pixels around location bandY, bandX.
 */
double DMVRasterBand::calcPixAvg(int bandY, int bandX, int maxy, int maxx)

{
	// 1 2 3
	// 4 X 5
	// 6 7 8
	double pixAvg = 0;
	double pixSum = 0;
	int pixCnt = 0;
	int x, y;

	// 1
	y = bandY - 1;
	x = bandX - 1;
	if ((y > 0) && (y < maxy) &&
		(x > 0) && (x < maxx)) {
		if (band[y][x] > 0) {
			pixSum += band[y][x];
			pixCnt++;
		}
	}
	// 2
	y = bandY - 1;
	x = bandX;
	if ((y > 0) && (y < maxy) &&
		(x > 0) && (x < maxx)) {
		if (band[y][x] > 0) {
			pixSum += band[y][x];
			pixCnt++;
		}
	}
	// 3
	y = bandY - 1;
	x = bandX + 1;
	if ((y > 0) && (y < maxy) &&
		(x > 0) && (x < maxx)) {
		if (band[y][x] > 0) {
			pixSum += band[y][x];
			pixCnt++;
		}
	}
	// 4
	y = bandY;
	x = bandX - 1;
	if ((y > 0) && (y < maxy) &&
		(x > 0) && (x < maxx)) {
		if (band[y][x] > 0) {
			pixSum += band[y][x];
			pixCnt++;
		}
	}
	// 5
	y = bandY;
	x = bandX + 1;
	if ((y > 0) && (y < maxy) &&
		(x > 0) && (x < maxx)) {
		if (band[y][x] > 0) {
			pixSum += band[y][x];
			pixCnt++;
		}
	}
	// 6
	y = bandY + 1;
	x = bandX - 1;
	if ((y > 0) && (y < maxy) &&
		(x > 0) && (x < maxx)) {
		if (band[y][x] > 0) {
			pixSum += band[y][x];
			pixCnt++;
		}
	}
	// 7
	y = bandY + 1;
	x = bandX;
	if ((y > 0) && (y < maxy) &&
		(x > 0) && (x < maxx)) {
		if (band[y][x] > 0) {
			pixSum += band[y][x];
			pixCnt++;
		}
	}
	// 8
	y = bandY + 1;
	x = bandX + 1;
	if ((y > 0) && (y < maxy) &&
		(x > 0) && (x < maxx)) {
		if (band[y][x] > 0) {
			pixSum += band[y][x];
			pixCnt++;
		}
	}

	pixAvg = pixSum / pixCnt;

	return pixAvg;
}


/************************************************************************/
/*                           DMVRasterBand()                            */
/************************************************************************/

DMVRasterBand::DMVRasterBand(DMVDataset *poDS, int nBand)

{
	this->poDS = poDS;
	this->nBand = nBand;

	eDataType = GDT_UInt32;

	nBlockXSize = poDS->GetRasterXSize();
	nBlockYSize = 1;

	/*
	 * All data is read in constructor, bacause it is not linear in file.
	 */
	char * line = NULL;
	char lline[80];
	int len = 0;
	XYZ * xyzs = NULL;
	struct xyz data;
	double maxx = 0;
	double minx = DBL_MAX;
	double maxy = 0;
	double miny = DBL_MAX;
	long lineCnt = 0;

	xyzs = (XYZ *) malloc(100 * sizeof *xyzs);

	// malloc band data and init to NO data value
	band = (double**) malloc(poDS->GetRasterYSize() * sizeof(double*));
	for (int i = 0; i < poDS->GetRasterYSize(); i++)
	{
		band[i] = (double*) malloc(poDS->GetRasterXSize() * sizeof(double));
		for (int j = 0; j < poDS->GetRasterXSize(); j++)
		{
			band[i][j] = (double)0;
		}
	}

	line = (char *) malloc(sizeof(char) * (80 + 1));

	VSIFSeekL( poDS->fp, 0, SEEK_SET );
	len = readLine(lline, poDS->fp);
	while ((len != EOF) && (len > 5))
	{
		lineCnt++;
		parseXyzLine(lline, &data);
		int bandY = (int) ((data.x - poDS->GetRasterYOrigin()) / 5);
		int bandX = (int) ((data.y - poDS->GetRasterXOrigin()) / 5);
		// TODO What to do with overshoots, skip or move down a pixel?
		if (bandY >= poDS->GetRasterYSize()) {
			printf("Line %d Band Y overshoot %d\n", lineCnt, bandY);
			bandY = poDS->GetRasterYSize() - 1;
		}
		if (bandX >= poDS->GetRasterXSize()) {
			printf("Line %d Band X overshoot %d\n", lineCnt, bandX);
			bandX = poDS->GetRasterXSize() - 1;
		}
		// There is duplicate erroneous data. Why is there such data?
		//printf("%ld : %ld, %ld = %.2lf\n", lineCnt, bandY, bandX, data.z);
		double pixAvg = calcPixAvg(bandY, bandX, poDS->GetRasterYSize(), poDS->GetRasterXSize());
		if (pixAvg == 0) {
			band[bandY][bandX] = data.z;
		} else {
			double delta = abs(data.z - pixAvg);
			if (delta > 30) {
				//printf("Delta %ld : %ld, %ld = %.2lf\n", lineCnt, bandY, bandX, delta);
			} else {
				band[bandY][bandX] = data.z;
			}
		}

		lline[0] = '\n';
		if (data.x > maxx)
		{
			maxx = data.x;
		}
		if (data.x < minx)
		{
			minx = data.x;
		}
		if (data.y > maxy)
		{
			maxy = data.y;
		}
		if (data.y < miny)
		{
			miny = data.y;
		}
		len = readLine(lline, poDS->fp);
	}
	printf("Minmax based on read data with x/y in source format (x ordinata, y abcisa):\n");
	printf("minx: %.2lf, maxx: %.2lf\n", minx, maxx);
	printf("miny: %.2lf, maxy: %.2lf\n", miny, maxy);
	printf("Read lines of data: %ld\n\n", lineCnt);

	poDS->minX = poDS->GetRasterXOrigin();
	poDS->maxX = poDS->GetRasterXOrigin() + (5 * poDS->GetRasterXSize());
	poDS->minY = poDS->GetRasterYOrigin();
	poDS->maxY = poDS->GetRasterYOrigin() + (5 * poDS->GetRasterYSize());

	printf("Minmax based on origin and size x/y in destination format (x abcisa, y ordinata):\n");
	printf("MINX: %.2lf, MAXX: %.2lf\n", poDS->minX, poDS->maxX);
	printf("MINY: %.2lf, MAXY: %.2lf\n\n", poDS->minY, poDS->maxY);

	if (line)
		free(line);
	if (xyzs)
		free(xyzs);
}

/************************************************************************/
/*                          ~DMVRasterBand()                            */
/************************************************************************/

DMVRasterBand::~DMVRasterBand()

{
	// Free band data
	for (int i = 0; i < this->poDS->GetRasterYSize(); i++)
	{
		free(band[i]);
	}
	free(band);
}

/************************************************************************/
/*                             IReadBlock()                             */
/************************************************************************/

CPLErr DMVRasterBand::IReadBlock( CPL_UNUSED int nBlockXOff,
		int nBlockYOff,
		void * pImage )

{
	DMVDataset *poGDS = (DMVDataset *) poDS;
	int i;

	//printf("\n*** IReadBlock %d, %d ***\n", nBlockXOff, nBlockYOff);

	for (i = 0; i < nBlockXSize; i++)
	{
		//printf("%.2lf, %ld\n", band[nBlockYOff][i], (int)(band[nBlockYOff][i]));
		// ta dela z gdaltranslate
		//((float *) pImage)[i] = (float) band[poGDS->GetRasterYSize()-nBlockYOff-1][i];
		// ta dela z gdalwarp
		int y = poGDS->GetRasterYSize() - nBlockYOff - 1;
		((int *) pImage)[i] = (int)(band[y][i]);
		//printf("\n*** IReadBlock band[%d, %d] = %d ***\n", y, i, (int)(band[nBlockYOff][i]));
	}

	return CE_None;
}

/************************************************************************/
/* ==================================================================== */
/*                             DMVDataset                               */
/* ==================================================================== */
/************************************************************************/

/************************************************************************/
/*                            DMVDataset()                             */
/************************************************************************/

DMVDataset::DMVDataset() :
		fp(NULL)
{
	fp = NULL;
}

/************************************************************************/
/*                            ~DMVDataset()                             */
/************************************************************************/

DMVDataset::~DMVDataset()

{
	FlushCache();
	if (fp != NULL)
		VSIFCloseL(fp);
}

/************************************************************************/
/*                          GetGeoTransform()                           */
/************************************************************************/

CPLErr DMVDataset::GetGeoTransform(double * padfTransform)

{

	padfTransform[0] = GetMinX(); // LLLong
	padfTransform[1] = 5; // pixel width
	padfTransform[2] = 0.0; // ?
	padfTransform[3] = GetMaxY(); // URLat
	padfTransform[4] = 0.0; // ?
	padfTransform[5] = -5; // pixel height


	return CE_None;
}

/************************************************************************/
/*                          GetProjectionRef()                          */
/************************************************************************/

const char *DMVDataset::GetProjectionRef()

{
	// http://epsg.io/3912	EPSG:3912	MGI 1901 / Slovene National Grid	D48/GK
	return ("PROJCS[\"MGI 1901 / Slovene National Grid\",GEOGCS[\"MGI 1901\",DATUM[\"MGI_1901\",SPHEROID[\"Bessel 1841\",6377397.155,299.1528128,AUTHORITY[\"EPSG\",\"7004\"]],TOWGS84[682,-203,480,0,0,0,0],AUTHORITY[\"EPSG\",\"1031\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",	\"8901\"]],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"3906\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",15],PARAMETER[\"scale_factor\",0.9999],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",-5000000],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Y\",EAST],AXIS[\"X\",NORTH],AUTHORITY[\"EPSG\",\"3912\"]]");
}

int DMVDataset::GetRasterXOrigin()

{
	return this->nRasterXOrigin;
}

int DMVDataset::GetRasterYOrigin()

{
	return this->nRasterYOrigin;
}

double DMVDataset::GetMinX()

{
	return this->minX;
}

double DMVDataset::GetMaxX()

{
	return this->maxX;
}

double DMVDataset::GetMinY()

{
	return this->minY;
}

double DMVDataset::GetMaxY()

{
	return this->maxY;
}

/************************************************************************/
/*                              Identify()                              */
/************************************************************************/

int DMVDataset::Identify(GDALOpenInfo *poOpenInfo)

{
	/* -------------------------------------------------------------------- */
	/*      First check file naming to confirm to:                          */
	/*      VT<letter><two digit num><two digit num>.xyz                    */
	/*                                                                      */
	/*      Check that first three rows have what it appears as correct     */
	/*      data of three double values separated by whitespace.            */
	/*      Relatively weak test.                                           */
	/* -------------------------------------------------------------------- */
	if (poOpenInfo->nHeaderBytes < 100)
		return FALSE;

	const char *fileName = CPLGetBasename(poOpenInfo->pszFilename);
	const char *fileExt = CPLGetExtension(poOpenInfo->pszFilename);

	/* Check file extension and naming. */
	if (!EQUALN(fileExt, "xyz", 3))
	{
		return FALSE;
	}
	if (!EQUALN(fileName, "VT", 2))
	{
		return FALSE;
	}
	char letter[1];
	strncpy(letter, fileName + 2, 1);
	if (!((letter[0] > 64) && (letter[0] < 77)))
	{
		return FALSE;
	}
	char numStr[2];
	strncpy(numStr, fileName + 3, 2);
	int checkNum = atoi(numStr);
	if (checkNum == 0)
	{
		return FALSE;
	}
	strncpy(numStr, fileName + 5, 2);
	checkNum = atoi(numStr);
	if (checkNum == 0)
	{
		return FALSE;
	}

	/* TODO Add first three line data check. */

	return TRUE;
}

/************************************************************************/
/*                                Open()                                */
/************************************************************************/

GDALDataset *DMVDataset::Open(GDALOpenInfo * poOpenInfo)

{
	/* -------------------------------------------------------------------- */
	/*      Confirm that the file name and first rows of data are           */
	/*      compatible with a DMV dataset.                                  */
	/* -------------------------------------------------------------------- */
	if (!Identify(poOpenInfo))
		return NULL;

	/* -------------------------------------------------------------------- */
	/*      Confirm the requested access is supported.                      */
	/* -------------------------------------------------------------------- */
	if (poOpenInfo->eAccess == GA_Update)
	{
		CPLError(CE_Failure, CPLE_NotSupported,
				"The DMV driver does not support update access to existing"
						" datasets.\n");
		return NULL;
	}

	/* Check that the file pointer from GDALOpenInfo* is available */
	if (poOpenInfo->fpL == NULL)
	{
		return NULL;
	}

	/* -------------------------------------------------------------------- */
	/*      Create a corresponding GDALDataset.                             */
	/* -------------------------------------------------------------------- */
	DMVDataset *poDS;

	poDS = new DMVDataset();

	/* Borrow the file pointer from GDALOpenInfo* */
	poDS->fp = poOpenInfo->fpL;
	poOpenInfo->fpL = NULL;

	/* Raster size is based on file name. */
	DmvCell cellDim;
	cellDim = getCellDimension(CPLGetBasename(poOpenInfo->pszFilename));
	poDS->nRasterYSize = cellDim.heightX;
	poDS->nRasterXSize = cellDim.widthY;
	poDS->nRasterYOrigin = cellDim.originX;
	poDS->nRasterXOrigin = cellDim.originY;

	if (poDS->nRasterXSize <= 0 || poDS->nRasterYSize <= 0)
	{
		CPLError(CE_Failure, CPLE_AppDefined, "Invalid dimensions : %d x %d",
				poDS->nRasterXSize, poDS->nRasterYSize);
		delete poDS;
		return NULL;
	}

	// Branje fileta
	//VSIFReadL(poDS->abyHeader, 1, 1012, poDS->fp);

	/* -------------------------------------------------------------------- */
	/*      Create band information objects.                                */
	/* -------------------------------------------------------------------- */
	poDS->SetBand(1, new DMVRasterBand(poDS, 1));

	/* -------------------------------------------------------------------- */
	/*      Initialize any PAM information.                                 */
	/* -------------------------------------------------------------------- */
	poDS->SetDescription(poOpenInfo->pszFilename);
	poDS->TryLoadXML();

	/* -------------------------------------------------------------------- */
	/*      Check for overviews.                                            */
	/* -------------------------------------------------------------------- */
	poDS->oOvManager.Initialize(poDS, poOpenInfo->pszFilename);

	return (poDS);
}

/************************************************************************/
/*                          GDALRegister_DMV()                          */
/************************************************************************/

void GDALRegister_DMV()

{
	GDALDriver *poDriver;

	if (GDALGetDriverByName("DMV") == NULL)
	{
		poDriver = new GDALDriver();

		poDriver->SetDescription("DMV");
		poDriver->SetMetadataItem(GDAL_DCAP_RASTER, "YES");
		poDriver->SetMetadataItem(GDAL_DMD_LONGNAME, "Slovenian DEM (.xyz)");
		poDriver->SetMetadataItem(GDAL_DMD_HELPTOPIC, "frmt_various.html#DMV");
		poDriver->SetMetadataItem(GDAL_DMD_EXTENSION, "xyz");
		poDriver->SetMetadataItem(GDAL_DCAP_VIRTUALIO, "YES");

		poDriver->pfnOpen = DMVDataset::Open;
		poDriver->pfnIdentify = DMVDataset::Identify;

		GetGDALDriverManager()->RegisterDriver(poDriver);
	}
}

/*
 * Read line of data till new line character.
 *
 * @param line buffer for data read
 * @param fp file handle
 * @return number of bytes read or EOF if end of file
 */
static int readLine(char line[], FILE* fp)

{
	int ch;
	char chr[1];
	int i = 0;

	VSIFReadL( chr, 1, 1, fp );
	while (chr[0] != EOF)
	{
		if (chr[0] == '\n')
		{
			// line read
			line[i++] = '\0';
			return i;
		}
		else if (chr[0] == '\r')
		{
			// skip
		}
		else
		{
			// char in line
			line[i++] = chr[0];
		}
		VSIFReadL( chr, 1, 1, fp );
	}

	return EOF;
}

/**
 * Parse line of data and return XYZ structure.
 *
 * @param line line of data from file.
 * @param data structure to be filed wiht data from line.
 */
static void parseXyzLine(char *line, struct xyz* data)
{
	int i = 0;
	int j = 0;
	char xStr[80];
	char yStr[80];
	char zStr[80];
	char ch;

	memset(yStr, 0, 80);
	while ((ch = line[i++]) != ' ')
	{
		yStr[j++] = ch;
	}
	j = 0;
	memset(xStr, 0, 80);
	while ((ch = line[i++]) != ' ')
	{
		xStr[j++] = ch;
	}
	j = 0;
	memset(zStr, 0, 80);
	while ((ch = line[i++]) != '\0')
	{
		zStr[j++] = ch;
	}

	data->x = atof(xStr);
	data->y = atof(yStr);
	data->z = atof(zStr);
}

/**
 * Sub Cell sizes are of different sizes depending if they are part of left border
 * upper border or both in case of left upper corner.
 * Those cells are one pixel wider or higer or both.
 *
 * @param subNum sub section num from file name.
 * @return XYZ, only x and y are used for heightX and widthY.
 */
static XYZ getCellSize(int subNum)

{

	// TODO Refactor because it looks like it is allways the same
	// this was based on reversing original data. Spec was bad.
	XYZ size;

	if (subNum == 1)
	{
		size.x = 600;
		size.y = 450;
	}
	else if (subNum > 1 && subNum < 11)
	{
		size.x = 600;
		size.y = 450; // 450
	}
	else if (subNum == 11 || subNum == 21 || subNum == 31 || subNum == 41)
	{
		size.x = 600;
		size.y = 450;
	}
	else
	{
		size.x = 600;
		size.y = 450; // 450
	}

	return size;
}

/**
 * Cell origin is based on TTN5 section name that is part of file name.
 * Inside section subSection is based on 5 x 10 subgrid.
 * Origin of sections is in lower left corner.
 * Sub sections are enumerated from 1 to 50, starting in upper left corner,
 * going 10 across and 5 down.
 *
 * @param letter section letter form file name.
 * @param num section num from file name.
 * @param subNum sub section num from file name.
 * @return XYZ, only x and y are used for x,y coordinates of cell origin.
 */
static XYZ getCellOrign(char letter, int num, int subNum)

{
	printf("Letter: %c, Num: %d, SubNum: %d\n", letter, num, subNum);
	// A - L, 19 - 30
	// 19A, 20A, .., 30A
	// 19B, 20B, .., 30B
	// ..
	// 19L, 20L, .., 30L
	// [][][0] - Y
	// [][][1] - X
	double originXY[12][12][2] =
	{
			// A
			{
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 365000, 85000 },
			{ 365000, 100000 },
			{ 365000, 115000 },
			{ 365000, 130000 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 } },
			// B
			{
			{ 387500, 25000 },
			{ 387500, 40000 },
			{ 387500, 55000 },
			{ 387500, 70000 },
			{ 387500, 85000 },
			{ 387500, 100000 },
			{ 387500, 115000 },
			{ 387500, 130000 },
			{ 387500, 145000 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 } },
			// C
			{
			{ 410000, 25000 },
			{ 410000, 40000 },
			{ 410000, 55000 },
			{ 410000, 70000 },
			{ 410000, 85000 },
			{ 410000, 100000 },
			{ 410000, 115000 },
			{ 410000, 130000 },
			{ 410000, 145000 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 } },
			// D
			{
			{ 432500, 25000 },
			{ 432500, 40000 },
			{ 432500, 55000 },
			{ 432500, 70000 },
			{ 432500, 85000 },
			{ 432500, 100000 },
			{ 432500, 115000 },
			{ 432500, 130000 },
			{ 432500, 145000 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 } },
			// E
			{
			{ 0, 0 },
			{ 455000, 40000 },
			{ 455000, 55000 },
			{ 455000, 70000 },
			{ 455000, 85000 },
			{ 455000, 100000 },
			{ 455000, 115000 },
			{ 455000, 130000 },
			{ 455000, 145000 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 } },
			// F
			{
			{ 477500, 25000 },
			{ 477500, 40000 },
			{ 477500, 55000 },
			{ 477500, 70000 },
			{ 477500, 85000 },
			{ 477500, 100000 },
			{ 477500, 115000 },
			{ 477500, 130000 },
			{ 477500, 145000 },
			{ 477500, 160000 },
			{ 0, 0 },
			{ 0, 0 } },
			// G
			{
			{ 500000, 25000 },
			{ 500000, 40000 },
			{ 500000, 55000 },
			{ 500000, 70000 },
			{ 500000, 85000 },
			{ 500000, 100000 },
			{ 500000, 115000 },
			{ 500000, 130000 },
			{ 500000, 145000 },
			{ 500000, 160000 },
			{ 0, 0 },
			{ 0, 0 } },
			// H
			{
			{ 522500, 25000 },
			{ 522500, 40000 },
			{ 522500, 55000 },
			{ 522500, 70000 },
			{ 522500, 85000 },
			{ 522500, 100000 },
			{ 522500, 115000 },
			{ 522500, 130000 },
			{ 522500, 145000 },
			{ 522500, 160000 },
			{ 0, 0 },
			{ 0, 0 } },
			// I
			{
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 545000, 70000 },
			{ 545000, 85000 },
			{ 545000, 100000 },
			{ 545000, 115000 },
			{ 545000, 130000 },
			{ 545000, 145000 },
			{ 545000, 160000 },
			{ 545000, 175000 },
			{ 0, 0 } },
			// J
			{
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 567500, 115000 },
			{ 567500, 130000 },
			{ 567500, 145000 },
			{ 567500, 160000 },
			{ 567500, 175000 },
			{ 567500, 190000 } },
			// K
			{
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 590000, 130000 },
			{ 590000, 145000 },
			{ 590000, 160000 },
			{ 590000, 175000 },
			{ 590000, 190000 } },
			// L
			{
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 0, 0 },
			{ 612500, 145000 },
			{ 612500, 160000 },
			{ 0, 0 },
			{ 0, 0 } } };

	int y = letter - 65;
	int x = num - 19;
	XYZ origin;

	// Origin based on TTN5
	if ((y > -1 && y < 12) && (x > -1 && x < 12))
	{
		origin.y = originXY[y][x][0]; // y orign
		origin.x = originXY[y][x][1]; // x orign
	}
	else
	{
		// Out of range
		origin.x = -1;
		origin.y = -1;
	}

	// Offset inside TTN5
	if (origin.x > 0 && origin.y > 0)
	{
		int subX = 4 - ((subNum - 1) / 10);
		int subY = (subNum - 1) % 10;

		if (subX > 0)
		{
			origin.x += (5 * 600 * subX);
		}
		if (subY > 0)
		{
			if (subY == 1)
			{
				origin.y += 5 * 450;
			}
			else
			{
				origin.y += (5 * 450) + (5 * 450 * (subY - 1));
			}
		}
	}

	return origin;
}

/**
 * Calculate origin of cell and its matrix size.
 *
 * @param fileName name of file.
 * @return structure with cell dimension.
 */
static DmvCell getCellDimension(const char * fileName)

{
	DmvCell dmvCell;
	XYZ xy;
	XYZ wh;
	char letter;
	char numStr[3];
	char subNumStr[3];
	int num;
	int subNum;

	letter = fileName[2];
	strncpy(numStr, &fileName[3], 2);
	numStr[2] = '\0';
	strncpy(subNumStr, &fileName[5], 2);
	subNumStr[2] = '\0';
	num = atoi(numStr);
	subNum = atoi(subNumStr);

	xy = getCellOrign(letter, num, subNum);
	wh = getCellSize(subNum);
	dmvCell.originX = xy.x;
	dmvCell.originY = xy.y;
	dmvCell.heightX = (int) wh.x;
	dmvCell.widthY = (int) wh.y;

	return dmvCell;
}
