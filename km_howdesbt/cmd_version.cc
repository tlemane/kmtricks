// cmd_version.cc-- report this program's version

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <iomanip>

#include "utilities.h"

#include "support.h"
#include "commands.h"
#include "cmd_version.h"

using std::string;
using std::cout;
using std::endl;

void VersionCommand::short_description
   (std::ostream& s)
	{
	s << commandName << "-- report this program's version" << endl;
	}

void VersionCommand::usage
   (std::ostream& s,
	const string& message)
	{
	if (!message.empty())
		{
		s << message << endl;
		s << endl;
		}

	short_description(s);
	s << "usage: " << commandName << endl;
	}

void VersionCommand::parse
   (int		_argc,
	char**	_argv)
	{
	if (_argc > 1) chastise("give me no arguments");
	return;
	}


int VersionCommand::execute()
	{
	cout << "version "
	     << major
	     << "." << std::setfill('0') << std::setw(2) << minor
	     << "." << std::setfill('0') << std::setw(2) << subMinor
	     << " " << std::setfill('0') << std::setw(8) << std::hex << date
	     << endl;

	return EXIT_SUCCESS;
	}

