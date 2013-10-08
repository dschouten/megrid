#include "logger.hh"

priority_t logger::msg_level = log4cpp::Priority::DEBUG;

log4cpp::Appender* logger::_stream = NULL;
log4cpp::Appender* logger::_file = NULL;

logger::logger() { }

logger::~logger() 
{
  if( _stream ) delete _stream;
  if( _file ) delete _file;
}

void logger::initialize( priority_t lvl, 
			 const std::string& logfile ) 
{  
  log4cpp::Category& root = log4cpp::Category::getRoot();
  root.setPriority( lvl );

  msg_level = lvl;

  if( _stream==0x0 )
  {
    _stream = new log4cpp::OstreamAppender("console", &std::cout);
    _stream->setLayout( new log4cpp::SimpleLayout() );
    root.addAppender( _stream );
  }
  if( _file==0x0 )
  {
    if( logfile != "" )
    {
      _file = new log4cpp::FileAppender("default", logfile.c_str());
      _file->setLayout( new log4cpp::SimpleLayout() );
      root.addAppender( _file );
    }
  }
}

void logger::set_level( priority_t lvl )
{
  if( _stream==0x0 )
  {
    initialize( lvl );
  }
  log4cpp::Category::getRoot().setPriority( lvl );
  msg_level = lvl;
}

log4cpp::Category& logger::log()
{
  if( _stream==0x0 )
  {
    initialize( msg_level );
  }
  return log4cpp::Category::getRoot();
}
