// howdesbt.cc-- work with HowDe sequence bloom trees.

#include <string>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <vector>
#include <list>

#include "utilities.h"
#include "support.h"
#include "commands.h"
// subcommands
#include "cmd_cluster.h"
#include "cmd_build_sbt.h"
#include "cmd_query_km.h"
#include "cmd_version.h"


using std::string;
using std::vector;
using std::list;
using std::pair;
using std::ostream;
using std::cout;
using std::cerr;
using std::endl;

string programName = "howdesbt";

//----------
//
// prototypes--
//
//----------

int main (int argc, char** argv);
//----------
//
// classes--
//
//----------

class MainCommand: public Command
	{
public:
	MainCommand(const string& name): Command(name) {}
	virtual ~MainCommand();
	virtual void short_description (ostream& s);
	virtual void usage (ostream& s, const string& message = "");
	virtual void usage_subcommands (ostream& s);
	virtual void parse (int _argc, char** _argv);
	virtual int execute (void);
	virtual void add_subcommand(Command* subCmd);
	virtual void add_command_alias(const string& name);
	virtual Command* find_subcommand(const string& name);

	vector<Command*> subCommands;
	vector<pair<string,Command*>> commandAliases;

	Command* subCommand = nullptr;
	int      subArgc = 0;
	char**   subArgv = nullptr;
	};

//----------
//
// main program--
//
//----------

int main
   (int		argc,
	char**	argv)
	{
	MainCommand* cmd = new MainCommand(programName);
	string programExe(argv[0]);

	// primary commands

	cmd->add_subcommand (new ClusterCommand      ("cluster"));
	cmd->add_subcommand (new BuildSBTCommand     ("build"));
	cmd->add_subcommand (new QueryCommandKm      ("queryKm"));
	cmd->add_subcommand (new VersionCommand      ("version"));

	// secondary commands


	// perform the user's command; if it was succesful, collect any additional
	// command(s) it would like us to perform

	cmd->parse(argc,argv);
	int successCode = cmd->execute ();

	list<string> toDoList;
	if ((successCode == EXIT_SUCCESS) and (cmd->subCommand != nullptr))
		{
		for (const auto& commandLine : cmd->subCommand->deferredCommands)
			toDoList.emplace_back(commandLine);
		cmd->subCommand->deferredCommands.clear();
		}

	// perform any additional commands; note that this may also generate
	// additional commands, which we'll put in the front of the queue

	while (not toDoList.empty())
		{
		string commandLine = toDoList.front();
		toDoList.pop_front();

		// $$$ move all this tokenization to a subroutine
		std::vector<string> argS = quoted_tokenize(commandLine);
		size_t argC = argS.size();

		char* argV[argC];
		for (size_t ix=0 ; ix<argC; ix++)
			{
			string arg = argS[ix];
			argV[ix] = new char[arg.size()+1];
			strcpy(/*to*/argV[ix],/*from*/arg.c_str());
			}

		cmd->parse(argC,argV);
		int successCode = cmd->execute ();
		for (size_t ix=0 ; ix<argC; ix++)
			delete[] argV[ix];
		if (successCode != EXIT_SUCCESS) break;

		if (cmd->subCommand != nullptr)
			{
			int numCommands = cmd->subCommand->deferredCommands.size();
			for (int cmdIx=numCommands-1 ; cmdIx>=0 ; cmdIx--)
				{
				string commandStr = cmd->subCommand->deferredCommands[cmdIx];
				toDoList.emplace_front(commandStr);
				}
			cmd->subCommand->deferredCommands.clear();
			}
		}

	delete cmd;
	return successCode;
	}

//----------
//
// main command--
//
//----------

MainCommand::~MainCommand()
	{
	for (const auto& subCmd : subCommands)
		delete subCmd;
	}

void MainCommand::short_description
   (ostream& s)
	{
	s << commandName << "-- work with HowDe sequence bloom trees" << endl;
	}

void MainCommand::usage
   (ostream& s,
	const string& message)
	{
	if (!message.empty())
		{
		s << message << endl;
		s << endl;
		}

	short_description(s);
	s << "usage: " << programName << " <command> [arguments]" << endl;
	//    123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
	s << "  <command>           perform a particular command" << endl;
	s << "  --help[=<command>]  get detail about a particular command" << endl;
	s << "  ?                   list available commands with brief descriptions" << endl;
	s << "  ?<command>          same as --help=<command>" << endl;
	s << endl;
	s << "For a list of available commands, do \"" << programName << " ?\"." << endl;
	s << "For more detailed descriptions of the commands, do \"" << programName << " --help\"." << endl;
	}

void MainCommand::usage_subcommands
   (ostream&	s)
	{
	vector<pair<string,string>> commandDescriptions;
	string prefix;
	string suffix;
	size_t maxPrefixLen = 0;

	// pre-scan command table, collecting prefix and suffix of each line

	for (const auto& subCmd : subCommands)
		{
		if (subCmd == nullptr)
			{
			commandDescriptions.emplace_back("","");
			continue;
			}

		std::stringstream ss;
		subCmd->short_description(ss);
		string description = strip_blank_suffix(strip_suffix(ss.str(),"\n"));
		size_t hyphensIx = description.find("--");
		if (hyphensIx == string::npos)
			{ prefix = description;  suffix = ""; }
		else
			{
			prefix = description.substr(0,hyphensIx+2);
			suffix = strip_blank_prefix(description.substr(hyphensIx+2));
			}
		commandDescriptions.emplace_back(prefix,suffix);
		if (prefix.length() > maxPrefixLen) maxPrefixLen = prefix.length();
		}

	// print command table, lining up the two columns

	s << "Primary commands (general form is <command> [arguments]):" << endl;

	for (const auto& prefixAndSuffix : commandDescriptions)
		{
		prefix = prefixAndSuffix.first;
		suffix = prefixAndSuffix.second;
		if (prefix == "")
			{
			s << endl << "Other commands (used less frequently):" << endl;
			continue;
			}
		s << std::setw(maxPrefixLen+1) << std::left << prefix << suffix << endl;
		}
	}

void MainCommand::parse
   (int		_argc,
	char**	_argv)
	{
	int		argc;
	char**	argv;

	// skip command name

	argv = _argv+1;  argc = _argc - 1;

	if (argc <= 0)
		{ usage (cerr);  std::exit (EXIT_FAILURE); }

	//////////
	// scan arguments
	//////////

	for (int argIx=0 ; argIx<argc ; argIx++)
		{
		string arg = argv[argIx];
		string argVal;
		if (arg.empty()) continue;

		string::size_type argValIx = arg.find('=');
		if (argValIx == string::npos) argVal = "";
		                              else argVal = arg.substr(argValIx+1);

		// sub command

		if (!is_prefix_of (arg, "--"))
			{
			Command* subCmd = find_subcommand (arg);
			if (subCmd != nullptr)
				{
				subCommand = subCmd;
				subArgc = argc-argIx;
				subArgv = argv+argIx;
				return;
				}
			}

		// --help and --help=<operator> (and variations starting with '?')

		if ((arg == "?")
		 || (arg == "--?"))
			{
			if (argIx == argc-1)
				{
				usage_subcommands (cerr);
				exit (EXIT_SUCCESS);
				}

			vector<pair<string,string>> commandDescriptions;
			string prefix;
			string suffix;
			size_t maxPrefixLen = 0;

			for (int helpIx=argIx+1 ; helpIx<argc ; helpIx++)
				{
				string subCmdName = argv[helpIx];
				Command* subCmd = find_subcommand (subCmdName);
				if (subCmd == nullptr)
					chastise ("\"" + subCmdName + "\" is not a known command");
				std::stringstream ss;
				subCmd->short_description(ss);
				string description = strip_blank_suffix(strip_suffix(ss.str(),"\n"));
				size_t hyphensIx = description.find("--");
				if (hyphensIx == string::npos)
					{ prefix = description;  suffix = ""; }
				else
					{
					prefix = description.substr(0,hyphensIx+2);
					suffix = strip_blank_prefix(description.substr(hyphensIx+2));
					}
				commandDescriptions.emplace_back(prefix,suffix);
				if (prefix.length() > maxPrefixLen) maxPrefixLen = prefix.length();
				}

			for (const auto& prefixAndSuffix : commandDescriptions)
				{
				prefix = prefixAndSuffix.first;
				suffix = prefixAndSuffix.second;
				cerr << std::setw(maxPrefixLen+1) << std::left << prefix << suffix << endl;
				}

			exit (EXIT_SUCCESS);
			}

		if (is_prefix_of (arg, "?="))
			goto alias_for_help_equals;

		if (is_prefix_of (arg, "?"))
			{ argVal = arg.substr(1);  goto help_for_one_subcommand; }

		if ((is_prefix_of (arg, "--help="))
		 || (is_prefix_of (arg, "--?=")))
			{
		alias_for_help_equals:
			if (argVal == "*") goto help_for_all_subcommands;
			if (argVal == "")  goto help_for_all_subcommands;
		help_for_one_subcommand:
			Command* subCmd = find_subcommand (argVal);
			if (subCmd == nullptr)
				chastise ("\"" + argVal + "\" is not a known command");

			cerr << "=== " << subCmd->commandName << " ===" << endl;
			subCmd->usage (cerr);
			exit (EXIT_SUCCESS);
			}

		if (arg == "--help")
			{
		help_for_all_subcommands:
		    for (const auto& subCmd : subCommands)
		    	{
		    	if (subCmd == nullptr) continue;
				cerr << "=== " << subCmd->commandName << " ===" << endl;
				subCmd->usage (cerr);
				}
			exit (EXIT_SUCCESS);
			}

		// --version

		if ((arg == "--version")
		 || (arg == "--v")
		 || (arg == "--V")
		 || (arg == "-v")
		 || (arg == "-V"))
			{
			VersionCommand cmd("version");
			cmd.execute();
			exit (EXIT_SUCCESS);
			}

		// unrecognized argument

		chastise ("unrecognized argument: \"" + arg + "\"");
		}

	return;
	}

void MainCommand::add_subcommand
   (Command*	subCmd)
	{
	subCommands.emplace_back(subCmd);
	}

void MainCommand::add_command_alias
   (const string& name)
	{
	if (subCommands.empty())
		fatal ("internal error: attempt to add alias \""
		     + name + "\" before any sub commands");

	Command* subCmd = subCommands.back();

	commandAliases.emplace_back(name,subCmd);
	}

Command* MainCommand::find_subcommand
   (const string& name)
	{
	for (const auto& subCmd : subCommands)
		{
	   	if (subCmd == nullptr) continue;
		if (subCmd->commandName == name) return subCmd;
		}

	for (const auto& nameAndCommand : commandAliases)
		{ if (nameAndCommand.first == name) return nameAndCommand.second; }

	return nullptr;
	}

int MainCommand::execute()
	{
	if (subCommand == nullptr)
		return 0; // success
	else
		return subCommand->main (subArgc, subArgv);
	}
