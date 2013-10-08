#ifndef LOG_HH
#define LOG_HH

#include <iostream>
#include <string>

#include "log4cpp/Category.hh"
#include "log4cpp/Appender.hh"
#include "log4cpp/FileAppender.hh"
#include "log4cpp/OstreamAppender.hh"
#include "log4cpp/Layout.hh"
#include "log4cpp/BasicLayout.hh"
#include "log4cpp/PatternLayout.hh"
#include "log4cpp/SimpleLayout.hh"
#include "log4cpp/Priority.hh"

typedef log4cpp::Priority::PriorityLevel priority_t;

class logger 
{
public:
  static priority_t msg_level;

  static void initialize( priority_t lvl = logger::msg_level, 
			  const std::string& logfile = "" );

  static void set_level( priority_t lvl );

  static log4cpp::Category& log();
  
private:
  logger();
  logger( const logger& ) {}
  logger& operator=( const logger& ) { return (*this); }

  ~logger();

  static log4cpp::Appender* _stream;
  static log4cpp::Appender* _file;

};

typedef log4cpp::Priority msg;

#endif
