#ifndef	PCREADER_H
#define PCREADER_H

#include <QString>
#include <QFileInfo>

#include "PointCloud.h"

enum PC_TYPE
{
	UNKNOWN_PC,
	ASCII_PC,
	LAS_PC
};

class PCReader
{
public:
	static bool readPCFile(PointCloud* pc, const QString& filename);

protected:

	virtual bool loadPC(PointCloud* pc, const QString& filename) = 0;

	static PC_TYPE getPCType(const QString& filename);

	static PCReader* getPCReader(PC_TYPE type);
};

#endif