#ifndef cmd_build_sbt_H
#define cmd_build_sbt_H

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>

#include "commands.h"

class BuildSBTCommand: public Command
	{
public:
	BuildSBTCommand(const std::string& name): Command(name) {}
	virtual ~BuildSBTCommand() {}
	virtual void short_description (std::ostream& s);
	virtual void usage (std::ostream& s, const std::string& message="");
	virtual void parse (int _argc, char** _argv);
	virtual int execute (void);

	std::string inTreeFilename;
	std::string outTreeFilename;
	std::uint32_t bfKind;
	std::uint32_t compressor;
	};

#endif // cmd_build_sbt_H
