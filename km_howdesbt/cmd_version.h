#ifndef cmd_version_H
#define cmd_version_H

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>

#include "commands.h"

class VersionCommand: public Command
	{
public:
	static const unsigned int  major    = 2;
	static const unsigned int  minor    = 0;
	static const unsigned int  subMinor = 4;
	static const std::uint32_t date     = 0x20210430;

public:
	VersionCommand(const std::string& name): Command(name) {}
	virtual ~VersionCommand() {}
	virtual void short_description (std::ostream& s);
	virtual void usage (std::ostream& s, const std::string& message="");
	virtual void parse (int _argc, char** _argv);
	virtual int execute (void);
	};

#endif // cmd_version_H
