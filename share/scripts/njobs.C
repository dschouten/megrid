void njobs( const char* F, const char* T, unsigned N )
{
  Int_t nJobs, nEventsPerJob, nEventsTotal;
  TFile* fp = TFile::Open( F, "read" );
  if( fp != NULL && !(fp->IsZombie()) )
  {
    TTree* tp = (TTree*)fp->Get( T );
    if( tp != NULL )
    {
      nEventsTotal = tp->GetEntries();
    }
  }
  nJobs = 1;
  while( nEventsTotal / nJobs > N )
  {
    nJobs += 1;
  }
  std::cout << nJobs+1 << std::endl;
}
