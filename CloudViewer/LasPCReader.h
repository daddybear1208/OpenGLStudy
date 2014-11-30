#ifndef	LASPCREADER_H
#define LASPCREADER_H

#include "PCReader.h"

class LasPCReader :public PCReader
{
public:
	bool loadPC(PointCloud* pc, const QString& filename)
	{
		return true;
	}
};

#endif