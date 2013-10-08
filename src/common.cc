#include "common.hh"

#include <iostream>

banner::banner()
{
  std::cout << std::endl
  	    << " __           _    ___  ___ _____          __   _____ "		<< std::endl
  	    << "/ _|         | |   |  \\/  ||  ___|        /  | |  _  |"	<< std::endl
  	    << "| |_ __ _ ___| |_  | .  . || |__   __   __`| | | |/' |"		<< std::endl
  	    << "|  _/ _` / __| __| | |\\/| ||  __|  \\ \\ / / | | |  /| |"	<< std::endl
  	    << "| || (_| \\__ \\ |_  | |  | || |___   \\ V / _| |_\\ |_/ /"	<< std::endl
  	    << "|_| \\__,_|___/\\__| \\_|  |_/\\____/    \\_/  \\___(_)___/ "	<< std::endl << std::endl;
}

const banner _msg;
