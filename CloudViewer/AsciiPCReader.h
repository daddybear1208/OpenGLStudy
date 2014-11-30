#ifndef ASCIIPCREADER_H
#define ASCIIPCREADER_H
#include "PCReader.h"

#include <QTextStream>
#include <QStringList>

class AsciiPCReader:public PCReader
{
public:
	bool loadPC(PointCloud* pc, const QString& filename)
	{
		QFile file(filename);
		if (!file.open(QFile::ReadOnly))
		{
			return false;
		}

		QTextStream stream(&file);
		QString currLine;
		do 
		{
			currLine = stream.readLine();
			QStringList points = currLine.split(",", QString::SkipEmptyParts);
			if (points.size() < 3)
				continue;
			Eigen::Vector3f pt(points[0].toFloat(), points[1].toFloat(), points[2].toFloat());
			pc->addPoint(pt);
		} while (!currLine.isEmpty());
		return true;
	}
};


#endif