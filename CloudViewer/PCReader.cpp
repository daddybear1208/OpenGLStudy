#include "AsciiPCReader.h"
#include "LasPCReader.h"

#include <QDebug>

bool PCReader::readPCFile(PointCloud* pc, const QString& filename)
{
	PCReader* reader = getPCReader(getPCType(filename));
	if (!reader) return false;
	return reader->loadPC(pc, filename);
}

PCReader* PCReader::getPCReader(PC_TYPE type)
{
	switch (type)
	{
	case UNKNOWN_PC:
		return 0;
		break;
	case ASCII_PC:
		return new AsciiPCReader();
		break;
	case LAS_PC:
		return new LasPCReader();
		break;
	default:
		return 0;
		break;
	}
}

PC_TYPE PCReader::getPCType(const QString& filename)
{
	PC_TYPE type = UNKNOWN_PC;
	QFileInfo file(filename);

	QString suffix = file.suffix();

	if (suffix == "txt")
	{
		type = ASCII_PC;
	}
	else if (suffix == "las")
	{
		type = LAS_PC;
	}
	return type;
}

