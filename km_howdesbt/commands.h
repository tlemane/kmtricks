#ifndef commands_H
#define commands_H

#include <string>
#include <iostream>
#include <set>
#include "utilities.h"
//----------
//
// Command class--
//	Subclasses should override everything here, except (usually) the destructor
//	and main().
//
//----------

class Command
	{
public:
	Command(const std::string& name) { commandName = name; }
	virtual ~Command() {}
	virtual int main (int argc, char** argv)
		{ parse (argc, argv);  return execute (); }
	virtual void chastise (const std::string& message = "")
		{ usage (std::cerr, message);  std::exit (EXIT_FAILURE); }
	virtual void short_description (std::ostream& s) {}
	virtual void usage (std::ostream& s, const std::string& message = "") {}
	virtual void parse (int _argc, char** _argv) {}
	virtual bool in_debug(std::string keyword) { return contains(debug,keyword); }
	virtual int execute (void) { return 0; }

public:
	std::string commandName;
	std::set<std::string> debug;
	std::vector<std::string> deferredCommands;
	};

#endif // commands_H
